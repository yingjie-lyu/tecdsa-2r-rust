use bicycl::{CipherText, Mpz, RandGen, ClearText, QFI};
use std::result::Result;

use crate::utils::*;
use crate::setup::*;

pub struct Nim {
    crs: Crs,
    rng: RandGen,
}

impl Nim {
    /// create a new NIM instance
    pub fn new(crs: Crs) -> Self {
        let mut rng = RandGen::new();
        rng.set_seed(&Mpz::from(&Zq::random()));
        Nim { crs, rng }
    }

    /// generate a random number
    fn generate_random(&mut self) -> Result<Mpz, &'static str> {
        Ok(self.rng.random_mpz(&self.crs.cl.encrypt_randomness_bound()))
    }

    /// compute h0^r * h1^x
    fn compute_pe(&self, r: &Mpz, x: &Mpz) -> Result<QFI, &'static str> {
        let h0_r = self.crs.cl.power_of_h(r);
        let h1_x = self.crs.pk.exp(&self.crs.cl, x);
        Ok(h0_r.compose(&self.crs.cl, &h1_x))
    }

    /// encode A's input
    pub fn encode_A(&mut self, x: &Mpz) -> Result<(QFI, Mpz), &'static str> {
        let r = self.generate_random()?;
        let pe_A = self.compute_pe(&r, x)?;
        Ok((pe_A, r))
    }

    /// encode B's input
    pub fn encode_B(&mut self, y: &Mpz) -> Result<(CipherText, Mpz), &'static str> {
        let clear_text = ClearText::with_mpz(&self.crs.cl, y);
        let s = self.generate_random()?;
        let pe_B = self.crs.cl.encrypt_with_r(&self.crs.pk, &clear_text, &s);
        Ok((pe_B, s))
    }

    /// decode A's result
    pub fn decode_A(&self, pe_B: &CipherText, st_A: &Mpz, x: &Mpz) -> Result<QFI, &'static str> {
        let c1 = pe_B.c1();
        let c2 = pe_B.c2();

        let c1_exp_st_A = c1.exp(&self.crs.cl, st_A);
        let c2_exp_x = c2.exp(&self.crs.cl, x);
        let zA = c1_exp_st_A.compose(&self.crs.cl, &c2_exp_x);

        // compute label z
        let mut label_zA = zA.clone();
        label_zA.to_maximal_order(&self.crs.cl.q(), &self.crs.cl.discriminant_K(), true);
        label_zA.lift(&self.crs.cl.q());

        let label_zA_invert = label_zA.exp(&self.crs.cl, &Mpz::from(-1 as i64));
        let alpha = zA.compose(&self.crs.cl, &label_zA_invert);
        Ok(alpha)
    }

    /// decode B's result
    pub fn decode_B(&self, pe_A: &QFI, st_B: &Mpz) -> Result<QFI, &'static str> {
        let zB = pe_A.exp(&self.crs.cl, st_B);

        // compute label z
        let mut label_zB = zB.clone();
        label_zB.to_maximal_order(&self.crs.cl.q(), &self.crs.cl.discriminant_K(), true);
        label_zB.lift(&self.crs.cl.q());

        let zB_invert = zB.exp(&self.crs.cl, &Mpz::from(-1 as i64));
        let beta = label_zB.compose(&self.crs.cl, &zB_invert);
        Ok(beta)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use curv::{
        arithmetic::Converter,
        BigInt,
    };

    #[test]
    fn test_nim_basic() {
        let crs = Crs::new();
        let mut nim = Nim::new(crs);
        let x = Mpz::from(Zq::random());
        let y = Mpz::from(Zq::random());

        let (pe_A, st_A) = nim.encode_A(&x).unwrap();
        let (pe_B, st_B) = nim.encode_B(&y).unwrap();

        let zA = nim.decode_A(&pe_B, &st_A, &x).unwrap();
        let zB = nim.decode_B(&pe_A, &st_B).unwrap();
        println!("zA: {:?}", zA);
        println!("zB: {:?}", zB);

        let xy = x * y;

        let alpha_mpz = nim.crs.cl.dlog_in_F(&zA);
        let beta_mpz = nim.crs.cl.dlog_in_F(&zB);

        let alpha_Zq = Zq::from_bigint(&BigInt::from_bytes(&alpha_mpz.to_bytes()));
        let beta_Zq = Zq::from_bigint(&BigInt::from_bytes(&beta_mpz.to_bytes()));
        println!("alpha_mpz: {:?}", alpha_mpz);
        println!("beta_mpz: {:?}", beta_mpz);

        let alpha_add_beta_Zq = alpha_Zq + beta_Zq;
        println!("alpha_add_beta_Zq: {:?}", alpha_add_beta_Zq);

        let xy_Zq = Zq::from_bigint(&BigInt::from_bytes(&xy.to_bytes()));
        println!("xy_Zq: {:?}", xy_Zq);

        assert!(xy_Zq == alpha_add_beta_Zq);
    }

    #[test]
    fn test_nim_zero() {
        let crs = Crs::new();
        let mut nim = Nim::new(crs);
        let x = Mpz::from(0 as i64);
        let y = Mpz::from(9 as i64);

        let (pe_A, st_A) = nim.encode_A(&x).unwrap();
        let (pe_B, st_B) = nim.encode_B(&y).unwrap();

        let zA = nim.decode_A(&pe_B, &st_A, &x).unwrap();
        let zB = nim.decode_B(&pe_A, &st_B).unwrap();

        let xy = x * y;
        
        let alpha_mpz = nim.crs.cl.dlog_in_F(&zA);
        let beta_mpz = nim.crs.cl.dlog_in_F(&zB);

        let alpha_Zq = Zq::from_bigint(&BigInt::from_bytes(&alpha_mpz.to_bytes()));
        let beta_Zq = Zq::from_bigint(&BigInt::from_bytes(&beta_mpz.to_bytes()));
        println!("alpha_mpz: {:?}", alpha_mpz);
        println!("beta_mpz: {:?}", beta_mpz);

        let alpha_add_beta_Zq = alpha_Zq + beta_Zq;
        println!("alpha_add_beta_Zq: {:?}", alpha_add_beta_Zq);

        let xy_Zq = Zq::from_bigint(&BigInt::from_bytes(&xy.to_bytes()));
        println!("xy_Zq: {:?}", xy_Zq);

        assert!(xy_Zq == alpha_add_beta_Zq);
    }
}
