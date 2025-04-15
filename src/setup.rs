use bicycl::{CL_HSMqk, PublicKey, RandGen, Mpz};
use curv::{
    arithmetic::Converter,
    BigInt,
};
use sha2::{Digest, Sha256};

use crate::utils::{Zq, G};

const CL_DELTA_K_NBITS: usize = 1827;
const NIZK_KST: usize = 40;

#[derive(Clone)]
pub struct Crs {
    pub cl: CL_HSMqk,
    pub pk: PublicKey,
}

impl Crs {
    pub fn new() -> Self {
        let mut rng = RandGen::new();
        rng.set_seed(&Mpz::from(&Zq::random()));

        let cl = CL_HSMqk::with_rand_gen(
            &Mpz::from_bytes(&Zq::group_order().to_bytes()),
            1,
            CL_DELTA_K_NBITS,
            &mut rng,
            &(Mpz::from_bytes(&(BigInt::from(1) << 40).to_bytes())),
            false,
        );

        let sk = cl.secret_key_gen(&mut rng);
        let pk = cl.public_key_gen(&sk);

        Crs { cl, pk }
    }
}

pub struct SignatureHashFunction {
    prefix: Vec<u8>,
}

impl SignatureHashFunction {
    pub fn new(prefix: &[u8]) -> Self {
        Self {
            prefix: prefix.to_vec(),
        }
    }

    pub fn hash(&self, msg: &[u8]) -> Zq {
        let mut hasher = Sha256::new();
        hasher.update(&self.prefix);
        hasher.update(msg);
        let hash_result = hasher.finalize();
        
        // 将哈希结果转换为Zq类型
        let hash_bigint = BigInt::from_bytes(&hash_result);
        Zq::from(hash_bigint % Zq::group_order())
    }
}

pub struct HashFunction {
    prefix: Vec<u8>,
}

impl HashFunction {
    pub fn new(prefix: &[u8]) -> Self {
        Self {
            prefix: prefix.to_vec(),
        }
    }

    pub fn hash(&self, msg: &[u8]) -> Zq {
        let mut hasher = Sha256::new();
        hasher.update(&self.prefix);
        hasher.update(msg);
        let hash_result = hasher.finalize();
        
        // 将哈希结果转换为Zq类型
        let hash_bigint = BigInt::from_bytes(&hash_result);
        Zq::from(hash_bigint % Zq::group_order())
    }
}

pub struct PublicParams {
    pub n: u8,
    pub t: u8,
    pub g: G,                        // generator
    pub q: Zq,                       // order
    pub hsig: SignatureHashFunction, // signature hash function
    pub h1: HashFunction,            // hash function H1
    pub h2: HashFunction,            // hash function H2
    pub crs: Crs,                    // CRS
    pub cl_delta_k_nbits: usize,     // CL_DELTA_K_NBITS
    pub nizk_kst: usize,             // NIZK_KST
}

/// Setup function: generate threshold ECDSA signature scheme public parameters
pub fn setup(_security_param: u8, n: u8, t: u8) -> PublicParams {    
    // use secp256k1 elliptic curve parameters
    // create a scalar multiple of the Zq base point, as the base point g
    let scalar_one = Zq::from(1u64);
    let g = G::generator() * &scalar_one;
    let q = Zq::from(Zq::group_order());
    
    let hsig = SignatureHashFunction::new(b"THRESHOLD_ECDSA_SIGNATURE");

    let h1 = HashFunction::new(b"THRESHOLD_ECDSA_H1");
    let h2 = HashFunction::new(b"THRESHOLD_ECDSA_H2");
    
    let crs = Crs::new();
    
    PublicParams {
        n,
        t,
        g,
        q,
        hsig,
        h1,
        h2,
        crs,
        cl_delta_k_nbits: CL_DELTA_K_NBITS,
        nizk_kst: NIZK_KST,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bicycl::Mpz;

    #[test]
    fn test_setup() {
        let pp = setup(128, 3, 2);
        
        // verify the generator is not zero
        assert!(!pp.g.is_zero());
        
        // verify the order is the curve order
        assert_eq!( pp.q, Zq::from(Zq::group_order()));
        
        // test hash function
        let msg1 = b"test message 1";
        let msg2 = b"test message 2";
        
        // test H1 and H2
        let hash1_1 = pp.h1.hash(msg1);
        let hash1_2 = pp.h1.hash(msg1);
        let hash2_1 = pp.h1.hash(msg2);
        
        // same message hash value should be the same
        assert_eq!(hash1_1, hash1_2);
        // different message hash value should be different
        assert_ne!(hash1_1, hash2_1);
        
        // verify H1 and H2 are different hash functions
        let h1_result = pp.h1.hash(msg1);
        let h2_result = pp.h2.hash(msg1);
        assert_ne!(h1_result, h2_result);
        
        // test signature hash function
        let sig_hash1 = pp.hsig.hash(msg1);
        let sig_hash2 = pp.hsig.hash(msg1);
        assert_eq!(sig_hash1, sig_hash2);
        
        // verify CRS is generated successfully
        let test_mpz = Mpz::from(123456u64);
        let _ = pp.crs.pk.exp(&pp.crs.cl, &test_mpz);
        
        println!("Setup test completed successfully!");
    }
}

