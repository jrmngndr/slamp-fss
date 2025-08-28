use crate::{field_arithmetic::*, config::*, error::*, linalg_util::*, prg::*, util::*};
use ndarray::{Array1, Array2};
use std::collections::{HashMap, HashSet};

pub type SecretKey = (
    u128,           // x_shared
    u128,           // t_shared
    Vec<u128>,      // all_wi0
    Vec<u128>,      // all_wi1
    Vec<Vec<u128>>, // all_di
);

pub type SecretKeyG = (SecretKey, Vec<u128>); // (sk, g)

pub fn gen_key(points: &[(u128, u128)]) -> Result<(SecretKeyG, SecretKeyG)> {
    let (a, b): (Vec<_>, Vec<_>) = points.iter().cloned().unzip();

    let mut attempts = 1;
    loop {
        match gen_rand(&a) {
            Ok((state, sk1, sk2)) => {
                let (_, z1_shared, z2_shared, _z) = state;
                match de_rand(&z1_shared, &z2_shared, &a, &b) {
                    Ok(g) => {
                        let sk1_g = (sk1, g.clone());
                        let sk2_g = (sk2, g);
                        println!("Gen key attempts: {}", attempts);
                        return Ok((sk1_g, sk2_g));
                    }
                    Err(e) => {
                        println!("Failed to generate key: {:?}", e);
                        attempts += 1;
                        if attempts > 100 {
                            return Err(SlampFssError::KeyGenFailed { attempts });
                        }
                        continue;
                    }
                }
            }
            Err(e) => {
                println!("Failed in gen_rand: {:?}", e);
                attempts += 1;
                if attempts > 100 {
                    return Err(SlampFssError::KeyGenFailed { attempts });
                }
                continue;
            }
        }
    }
}

type StateDict = HashMap<u128, (u128, u128, u128, u128, u128, u128)>; //rb, (x_rb1, x_rb2, x_rb, t_rb1, t_rb2, t_rb)
type State = (Vec<u128>, StateDict);

type FinalState = (
    Vec<u128>, //Rn
    HashMap<u128, u128>, //z1_shared
    HashMap<u128, u128>, //z2_shared
    HashMap<u128, u128>, //z_values
);

enum ProcessLayerResult {
    Intermediate(State, Vec<u128>),
    Final(FinalState, Vec<u128>),
}

fn gen_rand(a: &[u128]) -> Result<(FinalState, SecretKey, SecretKey)> {
    let x_e = random_from_field(V);
    let xe_1_shared = random_from_field(V);
    let xe_2_shared = gf_add(x_e, xe_1_shared, V);

    let t_e1_shared = random_from_field(V);
    let mut t_e2_shared = random_from_field(V);
    while t_e1_shared == t_e2_shared {
        t_e2_shared = random_from_field(V);
    }
    let t_e = gf_add(t_e1_shared, t_e2_shared, V);

    let initial_data = (
        xe_1_shared.clone(),
        xe_2_shared.clone(),
        x_e.clone(),
        t_e1_shared.clone(),
        t_e2_shared.clone(),
        t_e,
    );
    let mut state = (vec![0u128], {
        let mut map = HashMap::new();
        map.insert(0u128, initial_data);
        map
    });

    let mut all_wi0 = Vec::new();
    let mut all_wi1 = Vec::new();
    let mut all_di = Vec::new();

    for i in 0..N {
        let (wi0, wi1, layer_result) = process_layer(i, a, state)?;
        all_wi0.push(wi0);
        all_wi1.push(wi1);

        match layer_result {
            ProcessLayerResult::Final((rn, z1_shared, z2_shared, z_values), d_solution) => {
                all_di.push(d_solution);
                let sk1 = (
                    xe_1_shared,
                    t_e1_shared,
                    all_wi0.clone(),
                    all_wi1.clone(),
                    all_di.clone(),
                );
                let sk2 = (xe_2_shared, t_e2_shared, all_wi0, all_wi1, all_di);
                return Ok(((rn, z1_shared, z2_shared, z_values), sk1, sk2));
            }
            ProcessLayerResult::Intermediate(new_state, d_solution) => {
                all_di.push(d_solution);
                state = new_state;
            }
        }
    }

    unreachable!()
}

fn de_rand(
    z1_shared: &HashMap<u128, u128>,
    z2_shared: &HashMap<u128, u128>,
    a: &[u128],
    b: &[u128],
) -> Result<Vec<u128>> {
    let t = a.len();
    let mut z = Array2::zeros((t, V_BAR));

    for (i, &ai) in a.iter().enumerate() {
        let ai = ai >> (128-N);
        let z1_bits = to_bit_vector(z1_shared[&ai], V_BAR);
        let z2_bits = to_bit_vector(z2_shared[&ai], V_BAR);

        for j in 0..V_BAR {
            z[[i, j]] = (z1_bits[j] ^ z2_bits[j]) as u128;
        }
    }

    let z_rank = gf2_matrix_rank(&z);

    let mut augmented = Array2::zeros((t, V_BAR + 1));
    for i in 0..t {
        for j in 0..V_BAR {
            augmented[[i, j]] = z[[i, j]];
        }
        augmented[[i, V_BAR]] = b[i];
    }

    let augmented_rank = gf2_matrix_rank(&augmented);
    if z_rank != augmented_rank {
        return Err(SlampFssError::RankMismatch {
            expected: z_rank,
            actual: augmented_rank,
        });
    }

    let b_u128 = Array1::from_vec(b.iter().map(|&x| from_u128(x, V_BAR)).collect());
    solve_binary_system(&z, &b_u128, V_BAR).ok_or(SlampFssError::NoSolution)
}

fn process_layer(
    i: usize,
    a: &[u128],
    state: State,
) -> Result<(u128, u128, ProcessLayerResult)> {
    let (r_i_minus_1, state_dict) = state;
    let field_size = if i < N - 1 { K } else { V_BAR };

    let mut wi0 = random_from_field(field_size);
    while wi0 == 0 {
        wi0 = random_from_field(field_size);
    }

    let mut wi1 = random_from_field(field_size);
    while wi0 == wi1 || wi1 == 0 {
        wi1 = random_from_field(field_size);
    }

    let (r_i_active, r_i_inactive, r_i_minus_1_double) = compute_sets(&r_i_minus_1, a, i);
    let (a_matrix, b_vector) = add_constraints(
        &r_i_inactive,
        &r_i_minus_1_double,
        field_size,
        &state_dict,
        wi0,
        wi1,
    )?;

    let binary_rank = gf2_matrix_rank(&a_matrix);
    if binary_rank != a_matrix.nrows() {
        println!("Rank mismatch: expected {}, actual {}", a_matrix.nrows(), binary_rank);
        println!("A matrix: {}", a_matrix );
        return Err(SlampFssError::NotFullRank);
    }

    let d_solution =
        solve_binary_system(&a_matrix, &b_vector, field_size).ok_or(SlampFssError::NoSolution)?;

    if i == N - 1 {
        let mut z_dict = HashMap::new();

        for &rb in &r_i_active {
            let r_parent = rb >> 1;
            let b_bit = rb & 1;

            if let Some((x_r1, x_r2, _, t_r1, t_r2, _)) = state_dict.get(&r_parent) {
                let x_r1_vec = to_bit_vector(*x_r1, field_size);
                let x_r2_vec = to_bit_vector(*x_r2, field_size);

                let inner_prod1 = inner_product(&x_r1_vec, &d_solution, field_size);
                let inner_prod2 = inner_product(&x_r2_vec, &d_solution, field_size);

                let w_val = if b_bit == 0 { wi0 } else { wi1 };

                let z_rb1 = gf_add(
                    inner_prod1,
                    gf_multiply(*t_r1, w_val, field_size),
                    field_size,
                );
                let z_rb2 = gf_add(
                    inner_prod2,
                    gf_multiply(*t_r2, w_val, field_size),
                    field_size,
                );
                let z_rb = gf_add(z_rb1, z_rb2, field_size);

                z_dict.insert(rb, (z_rb1, z_rb2, z_rb));
            }
        }

        let rn: Vec<u128> = r_i_active.into_iter().collect();
        let z1_shared: HashMap<u128, u128> =
            rn.iter().map(|&r| (r, z_dict[&r].0)).collect();
        let z2_shared: HashMap<u128, u128> =
            rn.iter().map(|&r| (r, z_dict[&r].1)).collect();
        let z_values: HashMap<u128, u128> =
            rn.iter().map(|&r| (r, z_dict[&r].2)).collect();

        Ok((
            wi0,
            wi1,
            ProcessLayerResult::Final((rn, z1_shared, z2_shared, z_values), d_solution),
        ))
    } else {
        let mut new_state_dict = HashMap::new();

        for &rb in &r_i_active {
            let r_parent = rb >> 1;
            let b_bit = rb & 1;

            if let Some((x_r1, x_r2, _, t_r1, t_r2, _)) = state_dict.get(&r_parent) {
                let x_r1_vec = to_bit_vector(*x_r1, field_size);
                let x_r2_vec = to_bit_vector(*x_r2, field_size);

                let inner_prod1 = inner_product(&x_r1_vec, &d_solution, field_size);
                let inner_prod2 = inner_product(&x_r2_vec, &d_solution, field_size);

                let w_val = if b_bit == 0 { wi0 } else { wi1 };

                let prg_input1 = gf_add(
                    inner_prod1,
                    gf_multiply(*t_r1, w_val, field_size),
                    field_size,
                );
                let prg_input2 = gf_add(
                    inner_prod2,
                    gf_multiply(*t_r2, w_val, field_size),
                    field_size,
                );

                let (x_rb1, t_rb1, x_rb2, t_rb2) = if i == N - 2 {
                    let (x1, t1) = prg(prg_input1, V_BAR, V_BAR);
                    let (x2, t2) = prg(prg_input2, V_BAR, V_BAR);
                    (x1, t1, x2, t2)
                } else {
                    let (x1, t1) = prg_intermediate(prg_input1);
                    let (x2, t2) = prg_intermediate(prg_input2);
                    (x1, t1, x2, t2)
                };

                let x_rb = gf_add(x_rb1, x_rb2, field_size);
                let t_rb = gf_add(t_rb1, t_rb2, field_size);

                if t_rb == 0 {
                    return Err(SlampFssError::InvalidState);
                }

                new_state_dict.insert(rb, (x_rb1, x_rb2, x_rb, t_rb1, t_rb2, t_rb));
            }
        }

        let new_state = (r_i_active.into_iter().collect(), new_state_dict);
        Ok((
            wi0,
            wi1,
            ProcessLayerResult::Intermediate(new_state, d_solution),
        ))
    }
}

//adds constraints to the linear system
fn add_constraints(
    r_i_inactive: &HashSet<u128>,
    r_i_minus_1_double: &HashSet<u128>,
    field_size: usize,
    state_dict: &StateDict,
    wi0: u128,
    wi1: u128,
) -> Result<(Array2<u128>, Array1<u128>)> {
    let num_equations = r_i_inactive.len() + r_i_minus_1_double.len();
    let mut a = Array2::zeros((num_equations, field_size));
    let mut b = vec![0u128; num_equations];
    let mut equation_idx = 0;

    for &rb in r_i_inactive {
        let r_parent = rb >> 1;
        let b_bit = rb & 1;

        if let Some((_, _, x_r, _, _, t_r)) = state_dict.get(&r_parent) {
            let x_r_vec = to_bit_vector(*x_r, field_size);
            for j in 0..field_size.min(x_r_vec.len()) {
                a[[equation_idx, j]] = x_r_vec[j] as u128;
            }

            let w_val = if b_bit == 0 { wi0 } else { wi1 };
            b[equation_idx] = gf_multiply(*t_r, w_val, field_size);
        }
        equation_idx += 1;
    }

    for &r in r_i_minus_1_double {
        if let Some((_, _, x_r, _, _, t_r)) = state_dict.get(&r) {
            let mut wr2 = random_from_field(field_size);
            while wr2 == wi1 || wr2 == wi0 || wr2 == 0 {
                wr2 = random_from_field(field_size);
            }

            let x_r_vec = to_bit_vector(*x_r, field_size);
            for j in 0..field_size.min(x_r_vec.len()) {
                a[[equation_idx, j]] = x_r_vec[j] as u128;
            }
            b[equation_idx] = gf_multiply(*t_r, wr2, field_size);
        }
        equation_idx += 1;
    }
    let b_array = Array1::from_vec(b);
    Ok((a, b_array))
}


//compute active and inactive nodes in process layer
fn compute_sets(
    r_i_minus_1: &[u128],
    a: &[u128],
    i: usize,
) -> (HashSet<u128>, HashSet<u128>, HashSet<u128>) {
    let mut r_i_active = HashSet::new();
    let mut r_i_inactive = HashSet::new();
    let mut r_i_minus_1_double = HashSet::new();

    for &r in r_i_minus_1 {
        for b in [0, 1] {
            let r_concat_b = (r << 1) | b;
            let mut prefix_match = false;

            for &aj in a {
                if starts_with_prefix_msb(aj, r_concat_b, i + 1) {
                    prefix_match = true;
                    break;
                }
            }

            if prefix_match {
                r_i_active.insert(r_concat_b);
            } else {
                r_i_inactive.insert(r_concat_b);
            }
        }

        let r_0 = r << 1;
        let r_1 = (r << 1) | 1;

        if r_i_active.contains(&r_0) && r_i_active.contains(&r_1) {
            r_i_minus_1_double.insert(r);
        }
    }

    (r_i_active, r_i_inactive, r_i_minus_1_double)
}