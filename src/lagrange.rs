use std::collections::BTreeMap;

use crate::utils::*;

pub fn compute_lagrange_coefficients(parties: &[Id]) -> Result<BTreeMap<Id, Zq>, Box<dyn std::error::Error>> {
    if parties.is_empty() {
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Empty parties list"
        )));
    }

    if parties.len() < 2 {
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Need at least 2 parties for interpolation"
        )));
    }

    let mut lagrange_coefficients = BTreeMap::new();
    
    for &i in parties {
        let mut num = Zq::from(1u64);
        let mut den = Zq::from(1u64);
        
        for &j in parties {
            if i != j {
                let j_zq = Zq::from(j as u64);
                let diff = Zq::from(j as i32 - i as i32);
                
                num = num * j_zq;
                den = den * diff;
            }
        }
        
        let den_inv = match den.invert() {
            Some(inv) => inv,
            None => {
                return Err(Box::new(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to compute modular inverse for party {}", i)
                )));
            }
        };

        let lambda = num * den_inv;
        
        lagrange_coefficients.insert(i, lambda);
    }
    
    Ok(lagrange_coefficients)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::*;
    use crate::setup::setup;
    use curv::BigInt;
    use curv::arithmetic::Converter;

    #[test]
    fn test_lagrange_coefficients() {
        let n = 3;
        let t = 2;
        
        let pp = setup(128, n, t);
        
        println!("generate secret shares...");
        let (_, x, share_map, _) = keygen_step_1(&pp);
        
        println!("original secret x: {:?}", x.to_bigint().to_hex());
        for (id, share) in &share_map {
            println!("share {}: {:?}", id, share.to_bigint().to_hex());
        }
        
        println!("test lagrange coefficients calculation...");
        let parties = vec![1, 2, 3];
        let coeffs = compute_lagrange_coefficients(&parties).unwrap();
        
        for (id, coeff) in &coeffs {
            let coeff_zq = Zq::from_bigint(&BigInt::from_bytes(&coeff.to_bytes()));
            println!("lagrange coefficient {}: {:?}", id, coeff_zq.to_bigint().to_hex());
        }
        
        assert_eq!(coeffs.len(), 3, "should have 3 lagrange coefficients");
        
        for (_, coeff) in &coeffs {
            assert!(!coeff.to_bytes().is_empty(), "lagrange coefficient should not be 0");
        }
        
        let mut sum = Zq::from(0u64);
        for (_id, coeff) in &coeffs {
            let coeff_zq = Zq::from_bigint(&BigInt::from_bytes(&coeff.to_bytes()));
            sum = sum + coeff_zq;
        }
        assert_eq!(sum, Zq::from(1u64), "sum of lagrange coefficients should be 1");
        
        println!("verify secret reconstruction...");
        let mut reconstructed = Zq::from(0u64);
        for (id, coeff) in &coeffs {
            let share = share_map.get(id).unwrap();
            let coeff_zq = Zq::from_bigint(&BigInt::from_bytes(&coeff.to_bytes()));
            println!("ID {}: share={:?}, coefficient={:?}", id, share.to_bigint().to_hex(), coeff_zq.to_bigint().to_hex());
            reconstructed = reconstructed + coeff_zq * share;
        }
        
        println!("reconstructed secret: {:?}", reconstructed.to_bigint().to_hex());
        println!("original secret: {:?}", x.to_bigint().to_hex());
        assert_eq!(reconstructed, x, "reconstructed secret should be equal to original secret");
        
        println!("\ntest using 2 parties...");
        let parties = vec![1, 2];
        let coeffs = compute_lagrange_coefficients(&parties).unwrap();
        
        assert_eq!(coeffs.len(), 2, "should have 2 lagrange coefficients");
        
        let mut sum = Zq::from(0u64);
        for (_id, coeff) in &coeffs {
            let coeff_zq = Zq::from_bigint(&BigInt::from_bytes(&coeff.to_bytes()));
            sum = sum + coeff_zq;
        }
        assert_eq!(sum, Zq::from(1u64), "sum of lagrange coefficients should be 1");
        
        let mut reconstructed = Zq::from(0u64);
        for (id, coeff) in &coeffs {
            let share = share_map.get(id).unwrap();
            let coeff_zq = Zq::from_bigint(&BigInt::from_bytes(&coeff.to_bytes()));
            reconstructed = reconstructed + coeff_zq * share;
        }
        
        assert_eq!(reconstructed, x, "when using 2 parties, the reconstructed secret should be equal to the original secret");
        
        println!("\ntest error case...");
        let parties = vec![];
        let result = compute_lagrange_coefficients(&parties);
        assert!(result.is_err(), "empty parties list should return error");
        
        let parties = vec![1];
        let result = compute_lagrange_coefficients(&parties);
        assert!(result.is_err(), "only 1 party should return error");
    }
}
