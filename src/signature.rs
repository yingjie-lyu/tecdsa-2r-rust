use std::collections::BTreeMap;

use crate::setup::*;
use crate::utils::*;
use crate::presign::*;

pub struct Signature {
    pub r: Zq,
    pub s: Zq
}

pub fn partial_sign(
    pp: &PublicParams, 
    public_key: &G,
    msg: &[u8], 
    K: &G, 
    presignM: BTreeMap<Id, PreSignM>,
    pw: Zq, 
    pu: Zq, 
    gamma_i: Zq
) -> (Zq, Zq, Zq) {
    let m = pp.hsig.hash(msg);

    let public_key_bytes = public_key.to_bytes(false).to_vec();

    let presignM_bytes = presignM.values().map(|m| m.to_bytes()).collect::<Vec<_>>();
    let presignM_bytes_concat = presignM_bytes.concat();
    // h1.hash(public_key||msg||presignM)
    let z = pp.h1.hash(&[public_key_bytes, msg.to_vec(), presignM_bytes_concat].concat());
    let z_bytes: Vec<u8> = z.to_bytes().to_vec();
    let y = pp.h2.hash(&z_bytes);

    // R = Kz+gy
    let R = K * &z + &pp.g * &y;
    let r =Zq::from(R.x_coord().unwrap());

    // 计算最终的表达式
    let w_i = &m * &gamma_i + r.clone() * pw;

    let u_i = y * &gamma_i + z * pu;

    (w_i, u_i, r)
}

pub fn parallel_partial_sign(
    output: Vec<(G, Zq, Zq, Zq, BTreeMap<Id, PreSignM>)>,
    pp: &PublicParams,
    pk: &G,
    msg: &[u8],
) -> (Vec<Zq>, Vec<Zq>, Vec<Zq>) {
    let mut results = Vec::new();
    for (K, pw, pu, gamma_i, presign_m_list) in output {
        results.push(partial_sign(pp, pk, msg, &K, presign_m_list, pw, pu, gamma_i));
    }

    let mut w_i_list = Vec::with_capacity(results.len());
    let mut u_i_list = Vec::with_capacity(results.len());
    let mut r_list = Vec::with_capacity(results.len());

    for (w_i, u_i, r) in results {
        w_i_list.push(w_i);
        u_i_list.push(u_i);
        r_list.push(r);
    }

    (w_i_list, u_i_list, r_list)
}

pub fn combine_sign(
    w_i_list: Vec<Zq>,
    u_i_list: Vec<Zq>,
    r_list: Vec<Zq>,
) -> Signature {
    let w = w_i_list.iter().sum::<Zq>();
    let u = u_i_list.iter().sum::<Zq>();
    let r = r_list[0].clone();

    let s = &w * &u.invert().unwrap();        // s = w * u^-1
    Signature { r, s }
}

pub fn verify_sign(
    pp: &PublicParams,
    pk: &G,
    msg: &[u8],
    signature: Signature,
) -> bool {
    let m = pp.hsig.hash(msg);
    let s_invert = signature.s.invert().unwrap();
    let u1 = &m * &s_invert * &pp.g;
    let u2 = &signature.r * &s_invert * pk;
    let r_verify = u1 + u2;
    signature.r == Zq::from(r_verify.x_coord().unwrap())
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use crate::setup::*;
    use crate::keygen::*;
    use crate::presign::*;
    use crate::signature::*;


    const N: u8 = 7;
    const T: u8 = 5;
    const MSG: &[u8] = "Hello, world!".as_bytes();

    #[tokio::test]
    async fn test_sign_verify() {
        println!("N: {}, T: {}", N, T);
        let now = Instant::now();
        let pp = setup(128, N, T);
        let elapsed = now.elapsed();
        println!("setup time: {:?}", elapsed);


        let now = Instant::now();
        let (_x, pk, _pub_shares_map, aux_map, share_map, st_x_i_map) = keygen_simulation(&pp).await.unwrap();
        let elapsed = now.elapsed();
        println!("keygen time: {:?}", elapsed);

        let participant = vec![1, 2, 3, 4, 5, 6, 7];
        let now = Instant::now();
        let output = presign_simulation(&participant, &pp, &share_map, &st_x_i_map, &aux_map).await;
        let elapsed = now.elapsed();
        println!("presign time: {:?}", elapsed);


        let now = Instant::now();
        let (w_i_list, u_i_list, r_list) = parallel_partial_sign(output, &pp, &pk, &MSG);
        let elapsed = now.elapsed();
        println!("signing time: {:?}", elapsed);

        let now = Instant::now();
        let signature = combine_sign(w_i_list, u_i_list, r_list);
        let elapsed = now.elapsed();
        println!("combine time: {:?}", elapsed);

        // verify
        let now = Instant::now();
        let verify = verify_sign(&pp, &pk, &MSG, signature);
        let elapsed = now.elapsed();
        println!("verify time: {:?}", elapsed);

        println!("verify: {:?}", verify);
    }
}