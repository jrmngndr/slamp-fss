use aes::{Aes128, cipher::{BlockEncrypt, KeyInit}};
use aes::cipher::generic_array::GenericArray;
use crate::{config::*};

//PRG call for final layers
pub fn prg(input: u128, output1_field_size: usize, output2_field_size: usize) -> (u128, u128) {

    let key_bytes = input.to_le_bytes();
    let key = GenericArray::from_slice(&key_bytes);
    let cipher = Aes128::new(key);

    let mut output_bytes = Vec::with_capacity(32);

    for i in 0u128..2 {
        let mut block = GenericArray::from(i.to_le_bytes());
        cipher.encrypt_block(&mut block);
        output_bytes.extend_from_slice(&block);
    }

    let output1_bytes = &output_bytes[0..16];
    let output2_bytes = &output_bytes[16..32];

    let output1_raw = u128::from_le_bytes(output1_bytes.try_into().unwrap());
    let output2_raw = u128::from_le_bytes(output2_bytes.try_into().unwrap());

    let output1 = mask_to_field(output1_raw, output1_field_size);
    let output2 = mask_to_field(output2_raw, output2_field_size);

    (output1, output2)
}

//PRG call for intermediate layers
pub fn prg_intermediate(input: u128) -> (u128, u128) {
    let key_bytes = input.to_le_bytes();
    let key = GenericArray::from_slice(&key_bytes);
    let cipher = Aes128::new(key);

    let mut output_bytes = Vec::with_capacity(16);
    let mut block = GenericArray::from([0u8; 16]);
    cipher.encrypt_block(&mut block);
    output_bytes.extend_from_slice(&block);

    let output1_raw = u128::from_le_bytes(output_bytes.try_into().unwrap());

    let output1 = mask_to_field(output1_raw >> 64, V);
    let output2 = mask_to_field(output1_raw, V);

    (output1, output2)
}