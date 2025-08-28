use crate::config::*;

pub fn starts_with_prefix_msb(aj: u128, prefix: u128, prefix_length: usize) -> bool {
    if prefix_length == 0 {
        return true;
    }
    if prefix_length > N {
        return false;
    }

    let aj_prefix = aj >> (128 - prefix_length);
    aj_prefix == prefix
}

pub fn extract_input_bits(x: u128) -> Vec<u8> {
    let mut bits = vec![0u8; N];
    for i in 0..N {
        bits[i] = ((x >> (N - 1 - i)) & 1) as u8;
    }
    bits
}

pub fn to_bit_vector(value: u128, length: usize) -> Vec<u8> { //LSB
    (0..length).map(|i| {
        if i < 128 && (value & (1u128 << i)) != 0 { 1 } else { 0 }
    }).collect()
}