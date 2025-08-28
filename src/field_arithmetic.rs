use crate::config::*;


pub fn gf_multiply(a: u128, b: u128, field_size: usize) -> u128 {
    if a == 0 || b == 0 {
        return 0;
    }

    let mut result = 0u128;
    let mut a_shifted = a;
    let mut b_remaining = b;
    let irreducible = irreducible_for_field(field_size);

    while b_remaining != 0 {
        if (b_remaining & 1) == 1 {
            result ^= a_shifted;
        }

        a_shifted <<= 1;

        if field_size<128 && (a_shifted >> field_size) != 0 {
            a_shifted ^= irreducible;
        }

        b_remaining >>= 1;
    }

    mask_to_field(result, field_size)
}

pub fn gf_add(a: u128, b: u128, field_size: usize) -> u128 {
    mask_to_field(a ^ b, field_size)
}

pub fn inner_product(x_bits: &[u8], d_vec: &[u128], field_size: usize) -> u128 {
    let mut result = 0u128;

    let min_len = x_bits.len().min(d_vec.len());
    for i in 0..min_len {
        if x_bits[i] != 0 {
            result = gf_add(result, d_vec[i], field_size);
        }
    }
    result
}