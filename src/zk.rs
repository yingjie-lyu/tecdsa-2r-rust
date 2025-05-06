use bicycl::{CipherText, ClearText, Mpz, RandGen, QFI};
use crate::setup::*;
use crate::utils::*;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use curv::{
    arithmetic::Converter,
    BigInt,
};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct CLProof {
    c_tiled: CipherText,
    V_tiled: G,
    sr: Mpz,
    sv: Zq,
}

impl CLProof {
    pub fn to_bytes(&self) -> Vec<u8> {
        let c_tiled_bytes = self.c_tiled.to_bytes();
        let V_tiled_bytes = self.V_tiled.to_bytes(true).to_vec();
        let sr_bytes = self.sr.to_bytes();
        let sv_bytes = self.sv.to_bytes().to_vec();
        [c_tiled_bytes, V_tiled_bytes, sr_bytes, sv_bytes].concat()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct PedProof {
    c_tiled: QFI,
    V_tiled: G,
    sr: Mpz,
    sv: Mpz,
}

impl PedProof {
    pub fn to_bytes(&self) -> Vec<u8> {
        let c_tiled_bytes = self.c_tiled.to_bytes();
        let V_tiled_bytes = self.V_tiled.to_bytes(true).to_vec();
        let sr_bytes = self.sr.to_bytes();
        let sv_bytes = self.sv.to_bytes();
        [c_tiled_bytes, V_tiled_bytes, sr_bytes, sv_bytes].concat()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NIZKAoK {
    r_bound: Mpz,
    v_bound: Mpz,
}

impl NIZKAoK {
    pub fn new (pp: &PublicParams) -> Self {
        let r_bound = Self::get_r_bound(&pp);
        let v_bound = Self::get_v_bound(&pp);
        Self { r_bound, v_bound }
    }

    pub fn prove_cl(&self, pp: &PublicParams, rng: &mut RandGen, c: &CipherText, V: &G, r: &Mpz, v: &Zq) -> CLProof {
        // 1)
        let r_tiled = rng.random_mpz(&self.r_bound);
        let v_tiled = Zq::random();

        // 2)
        let clear_text = ClearText::with_mpz(&pp.crs.cl, &Mpz::from(&v_tiled));
        let c_tiled = pp.crs.cl.encrypt_with_r(&pp.crs.pk, &clear_text, &r_tiled);
        let V_tiled = &pp.g * &v_tiled;

        // 3) hash(g0, g1, f, g, c(c0,c1), V, c_tiled(c_tiled0,c_tiled1), V_tiled)
        let h0_bytes = pp.crs.cl.h().to_bytes();
        let h1_bytes = pp.crs.pk.to_bytes();
        let f_bytes = pp.crs.cl.power_of_f(&Mpz::from(1u64)).to_bytes();
        let g_bytes = pp.g.to_bytes(false).to_vec();
        let c_bytes = c.to_bytes();
        let V_bytes = V.to_bytes(false).to_vec();
        let c_tiled_bytes = c_tiled.to_bytes();
        let V_tiled_bytes = V_tiled.to_bytes(false).to_vec();
        let e = Self::hash(&[h0_bytes, h1_bytes, f_bytes, g_bytes, c_bytes, V_bytes, c_tiled_bytes, V_tiled_bytes].concat());

        // 4) sr = r_tiled + e * r, sv = v_tiled + e * v
        let sr = r_tiled + Mpz::from(&e) * r; // TODO: 好像是整数
        let sv = v_tiled + &e * v;

        CLProof { c_tiled, V_tiled, sr, sv }
    }

    pub fn verify_cl(&self, pp: &PublicParams, c: &CipherText, V: &G, proof: &CLProof) -> bool {
        let h0_bytes = pp.crs.cl.h().to_bytes();
        let h1_bytes = pp.crs.pk.to_bytes();
        let f_bytes = pp.crs.cl.power_of_f(&Mpz::from(1u64)).to_bytes();
        let g_bytes = pp.g.to_bytes(false).to_vec();
        let c_bytes = c.to_bytes();
        let V_bytes = V.to_bytes(false).to_vec();
        let c_tiled_bytes = proof.c_tiled.to_bytes();
        let V_tiled_bytes = proof.V_tiled.to_bytes(false).to_vec();

        let e = Self::hash(&[h0_bytes, h1_bytes, f_bytes, g_bytes, c_bytes, V_bytes, c_tiled_bytes, V_tiled_bytes].concat());

        //check r < bound and sv < Zq::group_order()
        if proof.sr.to_bytes() > self.r_bound.to_bytes() || proof.sv.to_bytes().to_vec() > pp.crs.cl.q().to_bytes() {
            return false;
        }

        let clear_text = ClearText::with_mpz(&pp.crs.cl, &Mpz::from(&proof.sv));
        let c_temp: CipherText = pp.crs.cl.encrypt_with_r(&pp.crs.pk, &clear_text, &proof.sr);

        let c0_tiled_mut_c0_pow_e = proof.c_tiled.c1().compose(&pp.crs.cl, &c.c1().exp(&pp.crs.cl, &Mpz::from(&e))) ;

        let c1_tiled_mut_c1_pow_e = proof.c_tiled.c2().compose(&pp.crs.cl, &c.c2().exp(&pp.crs.cl, &Mpz::from(&e)));


        let g_mut_sv = &pp.g * &proof.sv;

        let V_tiled_mut_V_pow_e = &proof.V_tiled + V * &e;


        if c_temp.c1() != c0_tiled_mut_c0_pow_e || c_temp.c2() != c1_tiled_mut_c1_pow_e || g_mut_sv != V_tiled_mut_V_pow_e {
            return false;
        }

        true
    }

    pub fn prove_ped(&self, pp: &PublicParams, rng: &mut RandGen, c: &QFI, V: &G, r: &Mpz, v: &Zq) -> PedProof {
        let v_tiled = rng.random_mpz(&self.v_bound);
        let r_tiled = rng.random_mpz(&self.r_bound);

        let h0_r = pp.crs.cl.power_of_h(&r_tiled);
        let h1_v = pp.crs.pk.exp(&pp.crs.cl, &v_tiled);
        let c_tiled = h0_r.compose(&pp.crs.cl, &h1_v);

        let V_tiled = &pp.g * &Zq::from_bigint(&BigInt::from_bytes(&v_tiled.to_bytes()));

        let h0_bytes = pp.crs.cl.h().to_bytes();
        let h1_bytes = pp.crs.pk.to_bytes();
        let g_bytes = pp.g.to_bytes(false).to_vec();
        let c_bytes = c.to_bytes();
        let V_bytes = V.to_bytes(false).to_vec();
        let c_tiled_bytes = c_tiled.to_bytes();
        let V_tiled_bytes = V_tiled.to_bytes(false).to_vec();
        let e = Self::hash(&[h0_bytes, h1_bytes, g_bytes, c_bytes, V_bytes, c_tiled_bytes, V_tiled_bytes].concat());
        
        let sr = r_tiled + Mpz::from(&e) * r;
        let sv = v_tiled + Mpz::from(&e) * Mpz::from(v);
        
        PedProof { c_tiled, V_tiled, sr, sv }
    }

    pub fn verify_ped(&self, pp: &PublicParams, c: &QFI, V: &G, proof: &PedProof) -> bool {
        let h0_bytes = pp.crs.cl.h().to_bytes();
        let h1_bytes = pp.crs.pk.to_bytes();
        let g_bytes = pp.g.to_bytes(false).to_vec();
        let c_bytes = c.to_bytes();
        let V_bytes = V.to_bytes(false).to_vec();
        let c_tiled_bytes = proof.c_tiled.to_bytes();
        let V_tiled_bytes = proof.V_tiled.to_bytes(false).to_vec();

        let e = Self::hash(&[h0_bytes, h1_bytes, g_bytes, c_bytes, V_bytes, c_tiled_bytes, V_tiled_bytes].concat());

        //check sr < r_bound and sv < v_bound
        if proof.sr.to_bytes() > self.r_bound.to_bytes() || proof.sv.to_bytes() > self.v_bound.to_bytes() {
            return false;
        }

        let h0_sr_mut_h1_sv = pp.crs.cl.power_of_h(&proof.sr).compose(&pp.crs.cl, &pp.crs.pk.exp(&pp.crs.cl, &proof.sv));
        let c_tiled_mut_c_pow_e = proof.c_tiled.compose(&pp.crs.cl, &c.exp(&pp.crs.cl, &Mpz::from(&e)));

        let g_sv = &pp.g * &Zq::from_bigint(&BigInt::from_bytes(&proof.sv.to_bytes()));
        let V_tiled_mut_V_pow_e = &proof.V_tiled + V * &e;

        if h0_sr_mut_h1_sv != c_tiled_mut_c_pow_e || g_sv != V_tiled_mut_V_pow_e {
            return false;
        }

        true
    }
    

    fn hash(msg: &[u8]) -> Zq {
        let mut hasher = Sha256::new();
        hasher.update(b"THRESHOLD_ECDSA_SIGNATURE_ZKP".to_vec());
        hasher.update(msg);
        let hash_result = hasher.finalize();
        
        let hash_bigint = BigInt::from_bytes(&hash_result);
        Zq::from(hash_bigint % Zq::group_order())
    }

    fn get_r_bound(pp: &PublicParams) -> Mpz {
        
        // log(√(pq)) = log(pq)/2
        let log_sqrt_pq = pp.cl_delta_k_nbits as f64 / 2.0;
        // ceil and convert to Mpz type
        let tm = Mpz::from(log_sqrt_pq.ceil() as u64);

        let z = Mpz::from(pp.nizk_kst as u64 * 2) + tm;
        let z_value = z.to_string().parse::<u64>().unwrap_or(0);
        
        let result: Mpz = Mpz::from(2u64).pow(z_value);

        result * &pp.crs.cl.q()
    }

    fn get_v_bound(pp: &PublicParams) -> Mpz {
        // q ^ 2
        let q_square = &pp.crs.cl.q().pow(2 as u64);
        Mpz::from(2u64).pow(pp.nizk_kst as u64) * q_square
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bicycl::{Mpz, RandGen};
    use std::time::Instant;
    use crate::nim::*;
    use crate::setup::setup;

    #[test]
    fn test_nizk() {
        let pp = setup(128, 5, 3);
        let mut nim = Nim::new(pp.crs.clone());
        let mut rnd = RandGen::new();
        let x = Zq::random();
        let y = Zq::random();
        let x_mpz = Mpz::from(&x);
        let y_mpz = Mpz::from(&y);

        let (pe_A, st_A) = nim.encode_A(&x_mpz).unwrap();
        let (pe_B, st_B) = nim.encode_B(&y_mpz).unwrap();

        let now = Instant::now();
        let nizk = NIZKAoK::new(&pp);
        let elapsed = now.elapsed();
        println!("nizk new time: {:?}", elapsed);

        let V = &pp.g * &y;
        let now = Instant::now();
        let proof_cl = nizk.prove_cl(&pp, &mut rnd,  &pe_B, &V, &st_B, &y);
        let elapsed = now.elapsed();
        println!("prove_cl time: {:?}", elapsed);

        let now = Instant::now();
        let result_cl = nizk.verify_cl(&pp, &pe_B, &V, &proof_cl);
        let elapsed = now.elapsed();
        println!("verify_cl time: {:?}", elapsed);
        assert!(result_cl);

        let V = G::generator() * &x;
        let now = Instant::now();
        let proof_ped = nizk.prove_ped(&pp, &mut rnd, &pe_A, &V, &st_A, &x);
        let elapsed = now.elapsed();
        println!("prove_ped time: {:?}", elapsed);

        let now = Instant::now();
        let result_ped = nizk.verify_ped(&pp, &pe_A, &V, &proof_ped);
        let elapsed = now.elapsed();
        println!("verify_ped time: {:?}", elapsed);
        assert!(result_ped);
    }
}