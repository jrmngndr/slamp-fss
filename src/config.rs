pub const K: usize = 64;
pub const V: usize = 64;
pub const V_BAR: usize = 128; //final output bits
pub const N: usize = 64; //input bits

pub const IRREDUCIBLE_K: u128 = (1u128 << 64) | (1u128 << 4) | (1u128 << 3) | (1u128 << 1) | 1u128;
pub const IRREDUCIBLE_V_BAR: u128 = (1u128 << 7) | (1u128 << 2) | (1u128 << 1) | 1u128;

pub fn field_mask(field_size: usize) -> u128 {
    if field_size == 128 {
        u128::MAX
    } else {
        (1u128 << field_size) - 1
    }
}

pub fn irreducible_for_field(field_size: usize) -> u128 {
    match field_size {
        K => IRREDUCIBLE_K,
        V_BAR => IRREDUCIBLE_V_BAR,
        _ => panic!("Unsupported field size: {}", field_size),
    }
}

pub fn mask_to_field(value: u128, field_size: usize) -> u128 {
    value & field_mask(field_size)
}

pub fn from_u128(value: u128, field_size: usize) -> u128 {
    mask_to_field(value, field_size)
}

pub fn random_from_field(field_size: usize) -> u128 {
    let value = rand::random::<u128>();
    mask_to_field(value, field_size)
}
