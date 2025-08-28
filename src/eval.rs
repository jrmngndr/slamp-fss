use crate::{config::*, field_arithmetic::*, prg::*, util::*, key_gen::*};

pub fn eval_function(x: u128, sk_g: &SecretKeyG) -> u128 {
    let (sk_prime, g) = sk_g;
    let z_share = eval_rand(x, sk_prime);
    de_rand_eval(z_share, g)
}

pub fn eval_rand(x: u128, sk_prime: &SecretKey) -> u128 {
    let (x_shared, t_shared, all_wi0, all_wi1, all_di) = sk_prime;

    let x_bits = extract_input_bits(x);
    let mut x_shared_current = *x_shared;
    let mut t_shared_current = *t_shared;

    let field_size = K;

    //process intermediate layers
    for i in 0..(N - 2) {
        let x_i = x_bits[i];
        let wi0_i = all_wi0[i];
        let wi1_i = all_wi1[i];
        let d_vec = &all_di[i];

        let x_vec = to_bit_vector(x_shared_current, field_size);
        let inner_prod = inner_product(&x_vec, d_vec, field_size);
        let w_i_xi = if x_i == 0 { wi0_i } else { wi1_i };
        let prg_input = gf_add(inner_prod, gf_multiply(t_shared_current, w_i_xi, field_size), field_size);

        let (new_x, new_t) = prg_intermediate(prg_input);
        x_shared_current = new_x;
        t_shared_current = new_t;
    }

    // process transition layer
    let i = N - 2;
    let x_i = x_bits[i];
    let wi0_i = all_wi0[i];
    let wi1_i = all_wi1[i];
    let d_vec = &all_di[i];

    let x_vec = to_bit_vector(x_shared_current, field_size);
    let inner_prod = inner_product(&x_vec, d_vec, field_size);
    let w_i_xi = if x_i == 0 { wi0_i } else { wi1_i };
    let prg_input = gf_add(inner_prod, gf_multiply(t_shared_current, w_i_xi, field_size), field_size);

    let (x_shared_current_final, t_shared_current_final) = prg(prg_input, V_BAR, V_BAR);

    //process final layer
    let i = N - 1;
    let x_i = x_bits[i];
    let wi0_i = all_wi0[i];
    let wi1_i = all_wi1[i];
    let d_vec = &all_di[i];

    let x_vec = to_bit_vector(x_shared_current_final, V_BAR);
    let inner_prod = inner_product(&x_vec, d_vec, V_BAR);
    let w_i_xi = if x_i == 0 { wi0_i } else { wi1_i };

    gf_add(inner_prod, gf_multiply(t_shared_current_final, w_i_xi, V_BAR), V_BAR)
}

// derandomize the share
pub fn de_rand_eval(z_share: u128, g: &[u128]) -> u128 {
    let z_bits = to_bit_vector(z_share, V_BAR);
    let mut y = 0u128;
    let min_len = g.len().min(z_bits.len());
    for i in 0..min_len {
        if z_bits[i] == 1 {
            y = gf_add(y, g[i], V_BAR);
        }
    }
    y
}

pub fn eval_two_party(x: u128, sk1_g: &SecretKeyG, sk2_g: &SecretKeyG) -> u128 {
    let (sk1_prime, g1) = sk1_g;
    let (sk2_prime, _g2) = sk2_g;
    let z1_share = eval_rand(x, sk1_prime);
    let z2_share = eval_rand(x, sk2_prime);
    let z_combined = gf_add(z1_share, z2_share, V_BAR);
    de_rand_eval(z_combined, g1)
}


//the tests are (mostly) generated with claude ai
//right now N is hardcoded in the config and should be changed for full eval tests
#[cfg(test)]
mod tests {
    use std::collections::{HashMap, HashSet};
    use std::time::Instant;
    use super::*;

    #[test]
    fn test_slamp_fss_64bit_domain() {
        const NUM_POINTS: usize = 32;
        const DOMAIN_BITS: usize = 64; // 2^64 domain
        const SHIFT_AMOUNT: usize = 128 - DOMAIN_BITS; // Shift to MSB position

        println!("Testing SLAMP-FSS with {} points in 2^{} domain", NUM_POINTS, DOMAIN_BITS);

        let mut points = Vec::new();
        let mut used_inputs = HashSet::new();

        // Generate 32 unique random points in the 2^64 domain
        while points.len() < NUM_POINTS {
            // Generate random 64-bit input, shift to MSB position
            let input_64bit = rand::random::<u64>() as u128;
            let input_128bit = input_64bit << SHIFT_AMOUNT;

            if used_inputs.contains(&input_128bit) {
                continue;
            }
            used_inputs.insert(input_128bit);

            // Random output value (can be full 128-bit)
            let output = rand::random::<u128>() >> SHIFT_AMOUNT;
            points.push((input_128bit, output));
        }

        points.sort_by_key(|&(a, _)| a);

        println!("Sample inputs (64-bit): {:016x}, {:016x}, {:016x}",
                 points[0].0 >> SHIFT_AMOUNT,
                 points[1].0 >> SHIFT_AMOUNT,
                 points[2].0 >> SHIFT_AMOUNT);

        // Measure key generation time
        let start_time = Instant::now();
        let result = gen_key(&points);
        let gen_time = start_time.elapsed();

        match result {
            Ok((sk1_g, sk2_g)) => {
                println!("Key generation successful in {:?}", gen_time);

                // Test correctness on interpolation points
                println!("Testing correctness on interpolation points...");
                let eval_start = Instant::now();
                let mut all_correct = true;

                for (i, &(input, expected_output)) in points.iter().enumerate() {
                    let result = eval_two_party(input>>SHIFT_AMOUNT, &sk1_g, &sk2_g);
                    let correct = result == expected_output;

                    if !correct {
                        all_correct = false;
                        println!("MISMATCH at point {}: f({:016x}) = {:x}, expected {:x}",
                                 i, input >> SHIFT_AMOUNT, result, expected_output);
                    } else if i < 5 { // Show first 5 results
                        println!("CORRECT: f({:016x}) = {:x}",
                                 input >> SHIFT_AMOUNT, result);
                    }
                }

                let eval_time = eval_start.elapsed();
                println!("{} evaluations completed in {:?}", NUM_POINTS, eval_time);
                println!("Average evaluation time: {:?}", eval_time / NUM_POINTS as u32);

                assert!(all_correct, "Not all interpolation points were correct");

                const NUM_POINTS_ZERO: usize = 1000;
                let mut points = Vec::new();

                while points.len() < NUM_POINTS_ZERO {
                    // Generate random 64-bit input, shift to MSB position
                    let input_64bit = rand::random::<u64>() as u128;
                    let input_128bit = input_64bit << SHIFT_AMOUNT;

                    if used_inputs.contains(&input_128bit) {
                        continue;
                    }
                    used_inputs.insert(input_128bit);

                    // Random output value (can be full 128-bit)
                    points.push(input_128bit);
                }

                points.sort_by_key(|&a| a);
                // Test on random non-interpolation points
                println!("Testing on non-interpolation points...");
                let eval_start = Instant::now();
                for i in 0..NUM_POINTS_ZERO {
                    // let random_input_64 = rand::random::<u64>() as u128;
                    // let random_input_128 = random_input_64 << SHIFT_AMOUNT;

                    // // Skip if this happens to be one of our interpolation points
                    // if used_inputs.contains(&random_input_128) {
                    //     continue;
                    // }

                    let result = eval_two_party(points[i]>>SHIFT_AMOUNT, &sk1_g, &sk2_g);
                    if i < 3 { // Show first 3 results
                        println!("Non-interpolation: f({:016x}) = {:x}",
                                 points[i], result);
                    }
                }
                let eval_time = eval_start.elapsed();
                println!("{} non encrypted evaluations completed in {:?}", NUM_POINTS_ZERO, eval_time);
                println!("Average evaluation time: {:?}", eval_time / NUM_POINTS_ZERO as u32);
                print!("evals per second {}", NUM_POINTS_ZERO as f64 / eval_time.as_secs_f64());

                println!("✓ All tests passed!");
            }
            Err(e) => {
                panic!("Key generation failed: {:?}", e);
            }
        }
    }

        #[test]
        fn benchmark_keygen() {
            const NUM_ITERATIONS: usize = 1000;
            const NUM_POINTS: usize = 32;
            const DOMAIN_BITS: usize = 64;
            const SHIFT_AMOUNT: usize = 128 - DOMAIN_BITS;

            println!("Benchmarking key generation with {} iterations", NUM_ITERATIONS);

            let mut total_time = std::time::Duration::ZERO;
            let mut successful_generations = 0;
            let mut failed_generations = 0;

            for iteration in 0..NUM_ITERATIONS {
                // Generate fresh points for each iteration
                let mut points = Vec::new();
                let mut used_inputs = HashSet::new();

                while points.len() < NUM_POINTS {
                    let input_64bit = rand::random::<u64>() as u128;
                    let input_128bit = input_64bit << SHIFT_AMOUNT;

                    if used_inputs.contains(&input_128bit) {
                        continue;
                    }
                    used_inputs.insert(input_128bit);

                    let output = rand::random::<u128>() >> SHIFT_AMOUNT;
                    points.push((input_128bit, output));
                }

                points.sort_by_key(|&(a, _)| a);

                // Measure key generation time
                let start_time = Instant::now();
                let result = gen_key(&points);
                let gen_time = start_time.elapsed();

                match result {
                    Ok(_) => {
                        total_time += gen_time;
                        successful_generations += 1;
                    }
                    Err(e) => {
                        failed_generations += 1;
                        println!("Key generation failed at iteration {}: {:?}", iteration, e);
                    }
                }

                if (iteration + 1) % 100 == 0 {
                    println!("Completed {} iterations...", iteration + 1);
                }
            }

            println!("\n=== KEY GENERATION BENCHMARK RESULTS ===");
            println!("Total iterations: {}", NUM_ITERATIONS);
            println!("Successful generations: {}", successful_generations);
            println!("Failed generations: {}", failed_generations);
            println!("Success rate: {:.2}%", (successful_generations as f64 / NUM_ITERATIONS as f64) * 100.0);

            if successful_generations > 0 {
                let average_time = total_time / successful_generations as u32;
                println!("Total time for successful generations: {:?}", total_time);
                println!("Average key generation time: {:?}", average_time);
                println!("Key generations per second: {:.2}", 1.0 / average_time.as_secs_f64());
            }
        }

        #[test]
        fn benchmark_eval() {
            const NUM_ITERATIONS: usize = 1000;
            const NUM_POINTS: usize = 32;
            const DOMAIN_BITS: usize = 64;
            const SHIFT_AMOUNT: usize = 128 - DOMAIN_BITS;

            println!("Benchmarking evaluation with {} iterations", NUM_ITERATIONS);

            // First, generate a single key pair to use for all evaluations
            let mut points = Vec::new();
            let mut used_inputs = HashSet::new();

            while points.len() < NUM_POINTS {
                let input_64bit = rand::random::<u64>() as u128;
                let input_128bit = input_64bit << SHIFT_AMOUNT;

                if used_inputs.contains(&input_128bit) {
                    continue;
                }
                used_inputs.insert(input_128bit);

                let output = rand::random::<u128>() >> SHIFT_AMOUNT;
                points.push((input_128bit, output));
            }

            points.sort_by_key(|&(a, _)| a);

            println!("Generating keys for benchmark...");
            let (sk1_g, sk2_g) = match gen_key(&points) {
                Ok(keys) => keys,
                Err(e) => {
                    panic!("Failed to generate keys for benchmark: {:?}", e);
                }
            };
            println!("Keys generated successfully!");

            // Generate random inputs for evaluation (avoiding interpolation points)
            let interpolation_inputs: HashSet<u128> = points.iter().map(|(input, _)| *input).collect();
            let mut eval_inputs = Vec::new();

            while eval_inputs.len() < NUM_ITERATIONS {
                let input_64bit = rand::random::<u64>() as u128;
                let input_128bit = input_64bit << SHIFT_AMOUNT;

                // Skip if this is an interpolation point
                if interpolation_inputs.contains(&input_128bit) {
                    continue;
                }

                eval_inputs.push(input_128bit);
            }

            println!("Starting evaluation benchmark...");
            let start_time = Instant::now();

            // Perform evaluations
            for (i, &input) in eval_inputs.iter().enumerate() {
                let _result = eval_two_party(input >> SHIFT_AMOUNT, &sk1_g, &sk2_g);

                if (i + 1) % 100 == 0 {
                    println!("Completed {} evaluations...", i + 1);
                }
            }

            let total_time = start_time.elapsed();

            println!("\n=== EVALUATION BENCHMARK RESULTS ===");
            println!("Total evaluations: {}", NUM_ITERATIONS);
            println!("Total time: {:?}", total_time);
            println!("Average evaluation time: {:?}", total_time / NUM_ITERATIONS as u32);
            println!("Evaluations per second: {:.2}", NUM_ITERATIONS as f64 / total_time.as_secs_f64());
        }

        #[test]
        fn benchmark_single_party_eval() {
            const NUM_ITERATIONS: usize = 1000;
            const NUM_POINTS: usize = 32;
            const DOMAIN_BITS: usize = 64;
            const SHIFT_AMOUNT: usize = 128 - DOMAIN_BITS;

            println!("Benchmarking single-party evaluation with {} iterations", NUM_ITERATIONS);

            // Generate key pair
            let mut points = Vec::new();
            let mut used_inputs = HashSet::new();

            while points.len() < NUM_POINTS {
                let input_64bit = rand::random::<u64>() as u128;
                let input_128bit = input_64bit << SHIFT_AMOUNT;

                if used_inputs.contains(&input_128bit) {
                    continue;
                }
                used_inputs.insert(input_128bit);

                let output = rand::random::<u128>() >> SHIFT_AMOUNT;
                points.push((input_128bit, output));
            }

            points.sort_by_key(|&(a, _)| a);

            println!("Generating keys for single-party benchmark...");
            let (sk1_g, _sk2_g) = match gen_key(&points) {
                Ok(keys) => keys,
                Err(e) => {
                    panic!("Failed to generate keys for benchmark: {:?}", e);
                }
            };
            println!("Keys generated successfully!");

            // Generate random inputs for evaluation
            let interpolation_inputs: HashSet<u128> = points.iter().map(|(input, _)| *input).collect();
            let mut eval_inputs = Vec::new();

            while eval_inputs.len() < NUM_ITERATIONS {
                let input_64bit = rand::random::<u64>() as u128;
                let input_128bit = input_64bit << SHIFT_AMOUNT;

                if interpolation_inputs.contains(&input_128bit) {
                    continue;
                }

                eval_inputs.push(input_128bit);
            }

            println!("Starting single-party evaluation benchmark...");
            let start_time = Instant::now();

            // Perform single-party evaluations
            for (i, &input) in eval_inputs.iter().enumerate() {
                let _result = eval_function(input >> SHIFT_AMOUNT, &sk1_g);

                if (i + 1) % 100 == 0 {
                    println!("Completed {} single-party evaluations...", i + 1);
                }
            }

            let total_time = start_time.elapsed();

            println!("\n=== SINGLE-PARTY EVALUATION BENCHMARK RESULTS ===");
            println!("Total evaluations: {}", NUM_ITERATIONS);
            println!("Total time: {:?}", total_time);
            println!("Average evaluation time: {:?}", total_time / NUM_ITERATIONS as u32);
            println!("Evaluations per second: {:.2}", NUM_ITERATIONS as f64 / total_time.as_secs_f64());
        }

    
    #[test]
    fn test_full_domain_evaluation() {
        const NUM_POINTS: usize = 16;
        const DOMAIN_BITS: usize = 20;       
        const SHIFT_AMOUNT: usize = 128 - DOMAIN_BITS;
        const DOMAIN_SIZE: usize = 1 << DOMAIN_BITS;

        println!("Testing full domain evaluation with {} points in 2^{} domain ({} total evaluations)",
                 NUM_POINTS, DOMAIN_BITS, DOMAIN_SIZE);

        let mut points = Vec::new();
        let mut used_inputs = HashSet::new();

        // Generate random interpolation points
        while points.len() < NUM_POINTS {
            let input_bits = rand::random::<u32>() & ((1u32 << DOMAIN_BITS) - 1);
            let input_128bit = (input_bits as u128) << SHIFT_AMOUNT;

            if used_inputs.contains(&input_128bit) {
                continue;
            }
            used_inputs.insert(input_128bit);

            let output = rand::random::<u128>() & ((1u128 << 64) - 1); // 64-bit output
            points.push((input_128bit, output));
        }

        points.sort_by_key(|&(a, _)| a);

        println!("Sample interpolation points:");
        for i in 0..std::cmp::min(5, points.len()) {
            println!("  f({:05x}) = {:016x}",
                     points[i].0 >> SHIFT_AMOUNT, points[i].1);
        }

        // Generate keys
        println!("Generating keys...");
        let start_time = Instant::now();
        let (sk1_g, sk2_g) = match gen_key(&points) {
            Ok(keys) => keys,
            Err(e) => {
                panic!("Key generation failed: {:?}", e);
            }
        };
        let keygen_time = start_time.elapsed();
        println!("Key generation completed in {:?}", keygen_time);

        // Verify correctness on interpolation points first
        println!("Verifying interpolation points...");
        for &(input, expected_output) in &points {
            let result = eval_two_party(input >> SHIFT_AMOUNT, &sk1_g, &sk2_g);
            assert_eq!(result, expected_output,
                       "Interpolation point mismatch: f({:05x}) = {:016x}, expected {:016x}",
                       input >> SHIFT_AMOUNT, result, expected_output);
        }
        println!("✓ All interpolation points correct");

        // Create interpolation point lookup for quick checking
        let interpolation_map: HashMap<u128, u128> = points.iter()
            .map(|&(input, output)| (input >> SHIFT_AMOUNT, output))
            .collect();

        // Full domain evaluation
        println!("Starting full domain evaluation...");
        let eval_start = Instant::now();

        let mut interpolation_hits = 0;
        let mut non_interpolation_evaluations = 0;
        let mut sample_outputs = Vec::new();

        // Evaluate all possible inputs in the domain
        for i in 0..DOMAIN_SIZE {
            let input = i as u128;
            let result = eval_two_party(input, &sk1_g, &sk2_g);

            // Check if this is an interpolation point
            if let Some(&expected) = interpolation_map.get(&input) {
                assert_eq!(result, expected,
                           "Mismatch at interpolation point: f({:05x}) = {:016x}, expected {:016x}",
                           input, result, expected);
                interpolation_hits += 1;
            } else {
                non_interpolation_evaluations += 1;
            }

            // Collect some sample outputs for inspection
            if sample_outputs.len() < 10 && !interpolation_map.contains_key(&input) {
                sample_outputs.push((input, result));
            }

            // Progress indicator
            if (i + 1) % 50000 == 0 {
                let elapsed = eval_start.elapsed();
                let rate = (i + 1) as f64 / elapsed.as_secs_f64();
                println!("  Completed {}/{} evaluations ({:.1}%), rate: {:.0} eval/sec",
                         i + 1, DOMAIN_SIZE, 100.0 * (i + 1) as f64 / DOMAIN_SIZE as f64, rate);
            }
        }

        let total_eval_time = eval_start.elapsed();

        println!("\n=== FULL DOMAIN EVALUATION RESULTS ===");
        println!("Domain size: 2^{} = {} points", DOMAIN_BITS, DOMAIN_SIZE);
        println!("Interpolation points: {}", NUM_POINTS);
        println!("Interpolation hits during evaluation: {}", interpolation_hits);
        println!("Non-interpolation evaluations: {}", non_interpolation_evaluations);
        println!("Total evaluation time: {:?}", total_eval_time);
        println!("Average evaluation time: {:?}", total_eval_time / DOMAIN_SIZE as u32);
        println!("Evaluations per second: {:.0}", DOMAIN_SIZE as f64 / total_eval_time.as_secs_f64());
        println!("Throughput: {:.2} MB/s (assuming 16 bytes per eval)",
                 (DOMAIN_SIZE * 16) as f64 / (1024.0 * 1024.0) / total_eval_time.as_secs_f64());

        println!("\nSample non-interpolation outputs:");
        for &(input, output) in &sample_outputs {
            println!("  f({:05x}) = {:016x}", input, output);
        }

        // Verify we hit all interpolation points
        assert_eq!(interpolation_hits, NUM_POINTS,
                   "Should have hit all {} interpolation points, but only hit {}",
                   NUM_POINTS, interpolation_hits);

        assert_eq!(interpolation_hits + non_interpolation_evaluations, DOMAIN_SIZE,
                   "Total evaluations should equal domain size");

        println!("✓ Full domain evaluation completed successfully!");
    }

    #[test]
    fn test_full_domain_evaluation_small() {
        const NUM_POINTS: usize = 8;         // Fewer interpolation points
        const DOMAIN_BITS: usize = 16;       // 2^16 = 65,536 domain size
        const SHIFT_AMOUNT: usize = 128 - DOMAIN_BITS;
        const DOMAIN_SIZE: usize = 1 << DOMAIN_BITS;

        println!("Testing small full domain evaluation with {} points in 2^{} domain ({} total evaluations)",
                 NUM_POINTS, DOMAIN_BITS, DOMAIN_SIZE);

        // Generate interpolation points
        let mut points = Vec::new();
        let mut used_inputs = HashSet::new();

        while points.len() < NUM_POINTS {
            let input_bits = rand::random::<u16>() as u128;
            let input_128bit = input_bits << SHIFT_AMOUNT;

            if used_inputs.contains(&input_128bit) {
                continue;
            }
            used_inputs.insert(input_128bit);

            let output = rand::random::<u128>() & ((1u128 << 32) - 1); // 32-bit output
            points.push((input_128bit, output));
        }

        points.sort_by_key(|&(a, _)| a);

        // Generate keys
        println!("Generating keys...");
        let (sk1_g, sk2_g) = match gen_key(&points) {
            Ok(keys) => keys,
            Err(e) => panic!("Key generation failed: {:?}", e),
        };
        println!("✓ Keys generated");

        // Build interpolation map
        let interpolation_map: HashMap<u128, u128> = points.iter()
            .map(|&(input, output)| (input >> SHIFT_AMOUNT, output))
            .collect();

        // Full domain evaluation with detailed analysis
        println!("Performing full domain evaluation...");
        let eval_start = Instant::now();

        let mut zero_outputs = 0;
        let mut nonzero_outputs = 0;
        let mut max_output = 0u128;
        let mut output_histogram = HashMap::new();

        for i in 0u128..DOMAIN_SIZE as u128 {
            let result = eval_two_party(i, &sk1_g, &sk2_g);

            // Verify interpolation points
            if let Some(&expected) = interpolation_map.get(&i) {
                assert_eq!(result, expected);
            }

            // Collect statistics
            if result == 0 {
                zero_outputs += 1;
            } else {
                nonzero_outputs += 1;
                max_output = max_output.max(result);
            }

            // Count output frequency (first 16 bits only for histogram)
            let output_key = (result & 0xFFFF) as u16;
            *output_histogram.entry(output_key).or_insert(0) += 1;
        }

        let total_eval_time = eval_start.elapsed();

        println!("\n=== SMALL DOMAIN EVALUATION ANALYSIS ===");
        println!("Total evaluations: {}", DOMAIN_SIZE);
        println!("Zero outputs: {} ({:.2}%)", zero_outputs, 100.0 * zero_outputs as f64 / DOMAIN_SIZE as f64);
        println!("Non-zero outputs: {} ({:.2}%)", nonzero_outputs, 100.0 * nonzero_outputs as f64 / DOMAIN_SIZE as f64);
        println!("Maximum output value: {:016x}", max_output);
        println!("Evaluation time: {:?}", total_eval_time);
        println!("Evaluations per second: {:.0}", DOMAIN_SIZE as f64 / total_eval_time.as_secs_f64());

        println!("Most common output values (low 16 bits):");
        let mut histogram_vec: Vec<_> = output_histogram.iter().collect();
        histogram_vec.sort_by(|a, b| b.1.cmp(a.1));
        for (value, count) in histogram_vec.iter().take(10) {
            println!("  {:04x}: {} times ({:.3}%)", value, count, 100.0 * **count as f64 / DOMAIN_SIZE as f64);
        }

        println!("✓ Small domain evaluation completed successfully!");
    }

}