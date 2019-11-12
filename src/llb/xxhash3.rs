/*
   xxHash - Extremely Fast Hash algorithm
   Rust implementation of XXH3_64bits
   Copyright (C) 2019-present, Yann Collet.

   BSD 2-Clause License (http://www.opensource.org/licenses/bsd-license.php)

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

       Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
       Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following disclaimer
   in the documentation and/or other materials provided with the
   distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   input, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   You can contact the author at :
   - xxHash source repository : https://github.com/Cyan4973/xxHash


 This is a relatively simple, but fast XXH3_64bits single run implementation in Rust, 
 written by easyaspi314.

 This doesn't use SIMD code, although it does use some unsafe hacks to avoid
 annoying bounds checking. It is still very performant.

 There is a rust namespaced version of this, as well as with the cfg 

 It is highly recommended to run this in release mode. In debug mode, this
 runs about 400x slower. */


use std::convert::TryInto;

/// Minimum size for the secret
pub const XXH3_SECRET_SIZE_MIN: usize = 136;
/// Default secret size
pub const XXH_SECRET_DEFAULT_SIZE: usize = 192;

/// Just a typedef
pub type XXHash64Hash = u64;

// Primes
const PRIME32_1: u32 = 0b10011110001101110111100110110001u32;
const PRIME32_2: u32 = 0b10000101111010111100101001110111u32;
const PRIME32_3: u32 = 0b11000010101100101010111000111101u32;

const PRIME64_1: u64 = 0b1001111000110111011110011011000110000101111010111100101010000111u64;
const PRIME64_2: u64 = 0b1100001010110010101011100011110100100111110101001110101101001111u64;
const PRIME64_3: u64 = 0b0001011001010110011001111011000110011110001101110111100111111001u64;
const PRIME64_4: u64 = 0b1000010111101011110010100111011111000010101100101010111001100011u64;
const PRIME64_5: u64 = 0b0010011111010100111010110010111100010110010101100110011111000101u64;

// Various constants
const ACC_NB: usize = 64 / std::mem::size_of::<u64>();
const XXH_SECRET_CONSUME_RATE: usize = 8;
const STRIPE_LEN: usize = 64;

// do not align on 8, so that secret is different from scrambler
const XXH_SECRET_LASTACC_START: usize = 7;
const XXH_SECRET_MERGEACCS_START: usize = 11;

const XXH3_MIDSIZE_STARTOFFSET: usize = 3;
const XXH3_MIDSIZE_LASTOFFSET: usize = 9;

///
/// The default secret key. This is used for the default seed, and is also used to mix
/// with a numeric seed.
///
/// Minimum XXH3_SECRET_SIZE_MIN
///
#[rustfmt::skip]
const DEFAULT_SECRET: &[u8; XXH_SECRET_DEFAULT_SIZE] = &[
    0xb8, 0xfe, 0x6c, 0x39, 0x23, 0xa4, 0x4b, 0xbe, 0x7c, 0x01, 0x81, 0x2c, 0xf7, 0x21, 0xad, 0x1c,
    0xde, 0xd4, 0x6d, 0xe9, 0x83, 0x90, 0x97, 0xdb, 0x72, 0x40, 0xa4, 0xa4, 0xb7, 0xb3, 0x67, 0x1f,
    0xcb, 0x79, 0xe6, 0x4e, 0xcc, 0xc0, 0xe5, 0x78, 0x82, 0x5a, 0xd0, 0x7d, 0xcc, 0xff, 0x72, 0x21,
    0xb8, 0x08, 0x46, 0x74, 0xf7, 0x43, 0x24, 0x8e, 0xe0, 0x35, 0x90, 0xe6, 0x81, 0x3a, 0x26, 0x4c,
    0x3c, 0x28, 0x52, 0xbb, 0x91, 0xc3, 0x00, 0xcb, 0x88, 0xd0, 0x65, 0x8b, 0x1b, 0x53, 0x2e, 0xa3,
    0x71, 0x64, 0x48, 0x97, 0xa2, 0x0d, 0xf9, 0x4e, 0x38, 0x19, 0xef, 0x46, 0xa9, 0xde, 0xac, 0xd8,
    0xa8, 0xfa, 0x76, 0x3f, 0xe3, 0x9c, 0x34, 0x3f, 0xf9, 0xdc, 0xbb, 0xc7, 0xc7, 0x0b, 0x4f, 0x1d,
    0x8a, 0x51, 0xe0, 0x4b, 0xcd, 0xb4, 0x59, 0x31, 0xc8, 0x9f, 0x7e, 0xc9, 0xd9, 0x78, 0x73, 0x64,

    0xea, 0xc5, 0xac, 0x83, 0x34, 0xd3, 0xeb, 0xc3, 0xc5, 0x81, 0xa0, 0xff, 0xfa, 0x13, 0x63, 0xeb,
    0x17, 0x0d, 0xdd, 0x51, 0xb7, 0xf0, 0xda, 0x49, 0xd3, 0x16, 0x55, 0x26, 0x29, 0xd4, 0x68, 0x9e,
    0x2b, 0x16, 0xbe, 0x58, 0x7d, 0x47, 0xa1, 0xfc, 0x8f, 0xf8, 0xb8, 0xd1, 0x7a, 0xd0, 0x31, 0xce,
    0x45, 0xcb, 0x3a, 0x8f, 0x95, 0x16, 0x04, 0x28, 0xaf, 0xd7, 0xfb, 0xca, 0xbb, 0x4b, 0x40, 0x7e
];

//////////////////////////////////////////////////////////////////////////
///                             Utilities                              ///
//////////////////////////////////////////////////////////////////////////

///
/// Returns a subslice of a u8 slice.
///
/// This is a debug feature which uses bounds checking. In release mode, we
/// turn bounds checking off.
///
#[inline(always)]
#[cfg(debug_assertions)]
fn read<I>(slice: &[u8], index: I) -> &I::Output
where
    I: std::slice::SliceIndex<[u8]>,
{
    &slice[index]
}

///
/// Returns a subslice of a u8 slice.
///
/// Because we know that none of this will ever go out of bounds,
/// we can safely disable bounds checking here.
///
#[inline(always)]
#[cfg(not(debug_assertions))]
fn read<I>(slice: &[u8], index: I) -> &I::Output
where
    I: std::slice::SliceIndex<[u8]>,
{
    unsafe { slice.get_unchecked(index) }
}

///
/// Portably reads a 32-bit little endian integer from p.
///
#[inline(always)]
fn read32(input: &[u8], offset: usize) -> u32 {
    let p = read(input, offset..offset + 4);
    u32::from_le_bytes(p.try_into().unwrap())
}

///
/// Portably reads a 64-bit little endian integer from p.
///
#[inline(always)]
fn read64(input: &[u8], offset: usize) -> u64 {
    let p = read(input, offset..offset + 8);
    u64::from_le_bytes(p.try_into().unwrap())
}

//////////////////////////////////////////////////////////////////////////
///                         Hashing subroutines                        ///
//////////////////////////////////////////////////////////////////////////

///
/// Does a 64-bit to 128-bit multiply, then xor's the low bits of the product with
/// the high bits for a 64-bit result.
///
/// This uses the native u128 type that Rust provides and is enabled for 64-bit targets.
///
#[inline]
#[cfg(all(any(target_pointer_width = "64", feature = "u128_arith"), not(feature = "no_u128_arith")))]
fn folded_mult_128(lhs: u64, rhs: u64) -> u64 {
    let product = lhs as u128 * rhs as u128;
    (product as u64) ^ ((product >> 64) as u64)
}

///
/// Does a 64-bit to 128-bit multiply, then xor's the low bits of the product with
/// the high bits for a 64-bit result.
///
/// This is the manual version which only uses native arithmetic.
///
/// Despite the fact that Rust now supports u128 arithmetic on all targets, we still
/// want to use the manual routine for 32-bit because multiplication results in a call
/// to compiler runtime.
///
#[rustfmt::skip]
#[cfg(all(any(not(target_pointer_width = "64"), feature = "no_u128_arith"), not(feature = "u128_arith")))]
fn folded_mult_128(lhs: u64, rhs: u64) -> u64 {
    let lo_lo = (lhs & 0xFFFFFFFF)   * (rhs & 0xFFFFFFFF);
    let hi_lo = (lhs >> 32)          * (rhs & 0xFFFFFFFF) + (lo_lo >> 32);
    let lo_hi = (lhs & 0xFFFFFFFF)   * (rhs >> 32)        + (hi_lo & 0xFFFFFFFF);
    let hi_64 = (lhs >> 32)          * (rhs >> 32)        + (lo_hi >> 32) + (hi_lo >> 32);
    let lo_64 = (lo_lo & 0xFFFFFFFF) | (lo_hi << 32);
    hi_64 ^ lo_64
}

///
/// Performs a final mix on the hash to improve distribution.
///
fn avalanche(mut hash: u64) -> XXHash64Hash {
    hash ^= hash >> 37;
    hash = hash.wrapping_mul(PRIME64_3);
    hash ^= hash >> 32;
    hash
}

///
/// Mixes 16 bytes
///
#[inline(always)]
fn mix_16(input: &[u8], input_offset: usize, key: &[u8], key_offset: usize, seed: u64) -> u64 {
    folded_mult_128(
        read64(input, input_offset + 0) ^ read64(key, key_offset + 0).wrapping_add(seed),
        read64(input, input_offset + 8) ^ read64(key, key_offset + 8).wrapping_sub(seed),
    )
}

///
/// Combines two accumulators with two keys.
///
#[inline]
fn mix_accs(acc: &[u64; 8], acc_offset: usize, key: &[u8], key_offset: usize) -> u64 {
    folded_mult_128(
        acc[acc_offset + 0] ^ read64(key, key_offset),
        acc[acc_offset + 1] ^ read64(key, key_offset + 8),
    )
}

///
/// Combines 8 accumulators with keys into 1 finalized 64-bit hash.
///
#[inline]
fn merge_accs(acc: &[u64; 8], key: &[u8], key_offset: usize, start: u64) -> XXHash64Hash {
    let mut result64 = start; // last bytes
    result64 = result64.wrapping_add(mix_accs(acc, 0, key, key_offset + 0));
    result64 = result64.wrapping_add(mix_accs(acc, 2, key, key_offset + 16));
    result64 = result64.wrapping_add(mix_accs(acc, 4, key, key_offset + 32));
    result64 = result64.wrapping_add(mix_accs(acc, 6, key, key_offset + 48));
    avalanche(result64)
}

//////////////////////////////////////////////////////////////////////////
///                             Short keys                             ///
//////////////////////////////////////////////////////////////////////////

///
/// Hashes short keys from 1 to 3 bytes.
///
#[inline]
fn hash_64_1_to_3(input: &[u8], key: &[u8], seed: XXHash64Hash) -> XXHash64Hash {
    let byte1 = input[0] as u32;
    let byte2 = if input.len() > 1 { input[1] } else { input[0] } as u32;
    let byte3 = input[input.len() - 1] as u32;
    let combined = byte1 | (byte2 << 8) | (byte3 << 16) | ((input.len() as u32) << 24);
    let keyed = (combined as u64) ^ ((read32(key, 0) as u64).wrapping_add(seed));
    let mixed = keyed.wrapping_mul(PRIME64_1);
    avalanche(mixed)
}

///
/// Hashes short keys from 4 to 8 bytes.
///
#[inline]
fn hash_64_4_to_8(input: &[u8], key: &[u8], seed: XXHash64Hash) -> XXHash64Hash {
    let input_lo = read32(input, 0);
    let input_hi = read32(input, input.len() - 4);
    let input64 = (input_lo as u64).wrapping_add((input_hi as u64) << 32);
    let keyed = input64 ^ ((read64(key, 0) as u64).wrapping_add(seed));
    let mixed = (input.len() as u64)
        .wrapping_add((keyed ^ (keyed >> 51)).wrapping_mul(PRIME32_1 as u64));
    let mixed = (mixed ^ (mixed >> 47)).wrapping_mul(PRIME64_2);
    avalanche(mixed)
}

///
/// Hashes short keys from 9 to 16 bytes.
///
#[inline]
fn hash_64_9_to_16(input: &[u8], key: &[u8], seed: XXHash64Hash) -> XXHash64Hash {
    let input1 = read64(input, 0) ^ (read64(key, 0).wrapping_add(seed));
    let input2 = read64(input, input.len() - 8) ^ (read64(key, 8).wrapping_sub(seed));
    let acc = (input.len() as u64)
        .wrapping_add(input1)
        .wrapping_add(input2)
        .wrapping_add(folded_mult_128(input1, input2));
    avalanche(acc)
}

///
/// Hashes short keys that are less than or equal to 16 bytes.
///
#[inline]
fn hash_64_0_to_16(input: &[u8], key: &[u8], seed: XXHash64Hash) -> XXHash64Hash {
    if input.len() > 8 {
        hash_64_9_to_16(input, key, seed)
    } else if input.len() >= 4 {
        hash_64_4_to_8(input, key, seed)
    } else if input.len() != 0 {
        hash_64_1_to_3(input, key, seed)
    } else {
        0
    }
}

//////////////////////////////////////////////////////////////////////////
///                           Medium keys                              ///
//////////////////////////////////////////////////////////////////////////

///
/// Hashes data that is from 17 to 128 bytes long.
///
#[inline]
fn hash_64_17_to_128(input: &[u8], secret: &[u8], seed: XXHash64Hash) -> XXHash64Hash {
    let mut acc = (input.len() as u64).wrapping_mul(PRIME64_1);

    if input.len() > 96 {
        acc = acc.wrapping_add(mix_16(input, 48, secret, 24 * 4, seed));
        acc = acc.wrapping_add(mix_16(input, input.len() - 64, secret, 28 * 4, seed));
    }
    if input.len() > 64 {
        acc = acc.wrapping_add(mix_16(input, 32, secret, 16 * 4, seed));
        acc = acc.wrapping_add(mix_16(input, input.len() - 48, secret, 20 * 4, seed));
    }
    if input.len() > 32 {
        acc = acc.wrapping_add(mix_16(input, 16, secret, 8 * 4, seed));
        acc = acc.wrapping_add(mix_16(input, input.len() - 32, secret, 12 * 4, seed));
    }
    if input.len() > 16 {
        acc = acc.wrapping_add(mix_16(input, 0, secret, 0 * 4, seed));
        acc = acc.wrapping_add(mix_16(input, input.len() - 16, secret, 4 * 4, seed));
    }
    avalanche(acc)
}

///
/// Hashes data that is 129 to 240 bytes long.
///
#[inline]
fn hash_64_129_to_240(input: &[u8], secret: &[u8], seed: XXHash64Hash) -> XXHash64Hash {
    let mut acc = (input.len() as u64).wrapping_mul(PRIME64_1);
    let nb_rounds = input.len() / 16;

    for i in 0..8 {
        acc = acc.wrapping_add(mix_16(input, 16 * i, secret, 16 * i, seed));
    }

    acc = avalanche(acc);

    for i in 8..nb_rounds {
        acc = acc.wrapping_add(mix_16(
            input,
            16 * i,
            secret,
            (16 * (i - 8)) + XXH3_MIDSIZE_STARTOFFSET,
            seed,
        ));
    }
    acc = acc.wrapping_add(mix_16(
        input,
        input.len() - 16,
        secret,
        XXH3_SECRET_SIZE_MIN - XXH3_MIDSIZE_LASTOFFSET,
        seed,
    ));
    avalanche(acc)
}

//////////////////////////////////////////////////////////////////////////
///                            Large keys                              ///
//////////////////////////////////////////////////////////////////////////

///
/// This is the main loop for large data blocks (>240 bytes).
///
/// This is usually written in SIMD code.
///
/// This is a modified version of the long multiply method used in UMAC and FARSH.
/// The changes are that data and key are xored instead of added, and the original data
/// is directly added to the mix after the multiply to prevent multiply-by-zero issues.
///
#[inline(always)]
fn accumulate_512(
    acc: &mut [u64; ACC_NB],
    input: &[u8],
    input_offset: usize,
    key: &[u8],
    key_offset: usize,
) {
    for i in 0..ACC_NB {
        let data_val = read64(input, input_offset + (8 * i));
        let key_val = read64(key, key_offset + (8 * i));
        let data_key = data_val ^ key_val;
        acc[i] = acc[i].wrapping_add((data_key & 0xFFFFFFFFu64) * (data_key >> 32));
        acc[i] = acc[i].wrapping_add(data_val);
    }
}

///
/// Scrambles each lane. This is usually written in SIMD code, as it is usually part of the main loop.
///
#[inline(always)]
fn scramble_acc(acc: &mut [u64; ACC_NB], key: &[u8], key_offset: usize) {
    for i in 0..ACC_NB {
        let key_val = read64(key, key_offset + (8 * i));
        acc[i] ^= acc[i] >> 47;
        acc[i] ^= key_val;
        acc[i] = acc[i].wrapping_mul(PRIME32_1 as u64);
    }
}

///
/// Processes a full block.
///
#[inline]
fn accumulate(
    acc: &mut [u64; 8],
    input: &[u8],
    input_offset: usize,
    secret: &[u8],
    nb_stripes: usize,
) {
    for n in 0..nb_stripes {
        accumulate_512(
            acc,
            input,
            input_offset + (n * STRIPE_LEN),
            secret,
            n * XXH_SECRET_CONSUME_RATE,
        );
    }
}
///
/// Controls the long hash function.
///
fn hash_long_internal_loop(acc: &mut [u64; 8], input: &[u8], secret: &[u8]) {
    let nb_rounds = (secret.len() - STRIPE_LEN) / XXH_SECRET_CONSUME_RATE;
    let block_len = STRIPE_LEN * nb_rounds;
    let nb_blocks = input.len() / block_len;
    let nb_stripes = (input.len() - (block_len * nb_blocks)) / STRIPE_LEN;
    for n in 0..nb_blocks {
        accumulate(acc, input, n * block_len, secret, nb_rounds);
        scramble_acc(acc, secret, secret.len() - STRIPE_LEN);
    }
    accumulate(acc, input, nb_blocks * block_len, secret, nb_stripes);
    if input.len() % STRIPE_LEN != 0 {
        accumulate_512(
            acc,
            input,
            input.len() - STRIPE_LEN,
            secret,
            secret.len() - STRIPE_LEN - XXH_SECRET_LASTACC_START,
        );
    };
}

///
/// Hashes a long block of data.
///
fn hash_long_internal(input: &[u8], secret: &[u8]) -> XXHash64Hash {
    let mut acc = [
        PRIME32_3 as u64,
        PRIME64_1,
        PRIME64_2,
        PRIME64_3,
        PRIME64_4,
        PRIME32_2 as u64,
        PRIME64_5,
        PRIME32_1 as u64,
    ];
    hash_long_internal_loop(&mut acc, input, secret);
    merge_accs(
        &mut acc,
        secret,
        XXH_SECRET_MERGEACCS_START,
        (input.len() as u64).wrapping_mul(PRIME64_1),
    )
}

fn hash_64_long_default_secret(input: &[u8]) -> XXHash64Hash {
    hash_long_internal(input, DEFAULT_SECRET)
}

//////////////////////////////////////////////////////////////////////////
///                         Public functions                           ///
//////////////////////////////////////////////////////////////////////////

////
/// The 64 bit non-seeded hash function.
/// input: The data to hash.
/// returns: The 64-bit calculated hash value.
///
pub extern fn hash_64(input: &[u8]) -> XXHash64Hash {
    if input.len() <= 16 {
        hash_64_0_to_16(input, DEFAULT_SECRET, 0)
    } else if input.len() <= 128 {
        hash_64_17_to_128(input, DEFAULT_SECRET, 0)
    } else if input.len() <= 240 {
        hash_64_129_to_240(input, DEFAULT_SECRET, 0)
    } else {
        hash_64_long_default_secret(input)
    }
}

#[cfg(feature = "c_wrapper")]
use std::ffi::c_void;

///
/// C wrappers. These are ABI-compatible with the official implementation.
///
#[no_mangle]
#[cfg(feature = "c_wrapper")]
pub unsafe extern "C" fn XXH3_64bits(data: *const c_void, length: usize) -> XXHash64Hash {
    if data.is_null() {
        let tmp: &[u8; 0] = &[];
        hash_64(tmp)
    } else {
        hash_64(std::slice::from_raw_parts(data as *const u8, length))
    }
}

#[no_mangle]
#[cfg(feature = "c_wrapper")]
pub unsafe extern "C" fn XXH3_64bits_withSeed(
    data: *const c_void,
    length: usize,
    seed: XXHash64Hash,
) -> XXHash64Hash {
    if data.is_null() {
        let tmp: &[u8; 0] = &[];
        hash_64_with_seed(tmp, seed)
    } else {
        hash_64_with_seed(std::slice::from_raw_parts(data as *const u8, length), seed)
    }
}

#[no_mangle]
#[cfg(feature = "c_wrapper")]
pub unsafe extern "C" fn XXH3_64bits_withSecret(
    data: *const c_void,
    length: usize,
    secret: *const c_void,
    secret_size: usize,
) -> XXHash64Hash {
    assert!(!secret.is_null());
    assert!(secret_size >= XXH3_SECRET_SIZE_MIN);
    if data.is_null() {
        let tmp: &[u8; 0] = &[];
        hash_64_with_secret(
            tmp,
            std::slice::from_raw_parts(secret as *const u8, secret_size),
        )
    } else {
        hash_64_with_secret(
            std::slice::from_raw_parts(data as *const u8, length),
            std::slice::from_raw_parts(secret as *const u8, secret_size),
        )
    }
}
