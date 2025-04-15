use bicycl::{Mpz, RandGen, SecretKey};
use futures::SinkExt;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use rayon::prelude::*;

use round_based::rounds_router::{simple_store::{RoundInput, RoundInputError}, CompleteRoundError, RoundsRouter};
use round_based::{Delivery, Mpc, MpcParty, Outgoing, ProtocolMessage};
use thiserror::Error;
use round_based::simulation::Simulation;

use crate::utils::*;
use crate::spdz::*;
use crate::nim::*;
use crate::setup::*;
use crate::zk::*;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AuxInfoMsg {
    pub commitment: bicycl::CipherText,
    pub proof: CLProof,
}

#[derive(Debug, Error)]
pub enum Error<RE, SE> {
    #[error("error in sending message in round 1")]
    Round1Send(#[source] SE),
    #[error("error in receiving message in round 1")]
    Round1Recv(#[source] CompleteRoundError<RoundInputError, RE>),
    #[error("error in sending message in round 2")]
    Round2Send(#[source] SE),
    #[error("error in receiving message in round 2")]
    Round2Recv(#[source] CompleteRoundError<RoundInputError, RE>),
    #[error("error in sending message in round 3")]
    Round3Send(#[source] SE),
    #[error("error in receiving message in round 3")]
    Round3Recv(#[source] CompleteRoundError<RoundInputError, RE>),
    #[error("error in step 2")]
    Step2Error(String),
}

#[derive(Clone, Debug, PartialEq, ProtocolMessage, Serialize, Deserialize)]
pub enum DkgMsg {
    Pvss(PvssMsg),
    PowOpen(OpenPowerMsg),
    AuxInfo(AuxInfoMsg),
}

impl From<PvssMsg> for DkgMsg {
    fn from(msg: PvssMsg) -> Self {
        DkgMsg::Pvss(msg)
    }
}

impl From<OpenPowerMsg> for DkgMsg {
    fn from(msg: OpenPowerMsg) -> Self {
        DkgMsg::PowOpen(msg)
    }
}

impl From<AuxInfoMsg> for DkgMsg {
    fn from(msg: AuxInfoMsg) -> Self {
        DkgMsg::AuxInfo(msg)
    }
}


pub async fn dkg<M>(
    party: M,
    my_id: Id, // in the range 1..=n, subtract one before use
    pp: &PubParams,
    h: &G,
    my_cl_sk: &SecretKey,
    lazy_verification: bool,
    crs: &Crs,
) -> Result<(G, BTreeMap<Id, G>, Zq, BTreeMap<Id, bicycl::CipherText>, Mpz), Error<M::ReceiveError, M::SendError>>
where
    M: Mpc<ProtocolMessage = DkgMsg>,
{
    let mut rng = RandGen::new();
    rng.set_seed(&Mpz::from(&Zq::random()));

    // boilerplate
    let MpcParty { delivery, .. } = party.into_party();
    let (incoming, mut outgoing) = delivery.split();

    let i = (my_id - 1) as u16;
    let n = pp.n as u16;

    let mut rounds = RoundsRouter::<DkgMsg>::builder();
    let round1 = rounds.add_round(RoundInput::<PvssMsg>::broadcast(i, n));
    let round2 = rounds.add_round(RoundInput::<OpenPowerMsg>::broadcast(i, n));
    let round3 = rounds.add_round(RoundInput::<AuxInfoMsg>::broadcast(i, n));
    let mut rounds = rounds.listen(incoming);

    // Round 1 interaction
    let pvss_msg = PvssMsg::random(pp, &mut rng, h);
    outgoing
        .send(Outgoing::broadcast(DkgMsg::Pvss(pvss_msg.clone())))
        .await
        .map_err(Error::Round1Send)?;

    let pvss_messages = rounds.complete(round1).await.map_err(Error::Round1Recv)?;

    // Round 1 processing
    let mut pvss_dealings = BTreeMap::new();
    pvss_dealings.insert(my_id, pvss_msg.dealing);

    pvss_messages
        .into_iter_indexed()
        .map(|(inner_id, _, msg)| ((inner_id + 1) as Id, msg))
        .filter(|(_, msg)| lazy_verification || msg.proof.verify(&msg.dealing, pp, h))
        .take(pp.t as usize)
        .for_each(|(j, msg)| {
            pvss_dealings.insert(j, msg.dealing);
        });

    let pvss_result = JointPvssResult::new(
        pp,
        pvss_dealings
            .values()
            .take(pp.t as usize)
            .cloned()
            .collect(),
    );

    let my_share = pvss_result
        .shares_ciphertext
        .decrypt_mine(&pp.cl, my_id, my_cl_sk);

    let my_pub_share = G::generator() * &my_share;

    let dleq_proof = DleqNizk::prove(
        h,
        &pvss_result.curve_macs[&my_id],
        &G::generator(),
        &my_pub_share,
        &my_share,
    );

    let open_power_msg = OpenPowerMsg {
        point: my_pub_share.clone(),
        proof: dleq_proof,
    };

    // Round 2 interaction
    outgoing
        .send(Outgoing::broadcast(DkgMsg::PowOpen(open_power_msg.clone())))
        .await
        .map_err(Error::Round2Send)?;

    let open_power_messages = rounds.complete(round2).await.map_err(Error::Round2Recv)?;

    // Round 2 processing
    let mut pub_shares = BTreeMap::new();
    pub_shares.insert(my_id, open_power_msg.point);

    open_power_messages
        .into_iter_indexed()
        .map(|(inner_id, _, msg)| ((inner_id + 1) as Id, msg))
        .filter(|(id, msg)| {
            msg.proof
                .verify(h, &pvss_result.curve_macs[id], &G::generator(), &msg.point)
        })
        .for_each(|(j, msg)| {
            pub_shares.insert(j, msg.point);
        });

    let pk = pp.interpolate_on_curve(&pub_shares).unwrap();

    // Round 3 - Generate auxiliary information
    // 1. Compute (pe_{x,i}, st_{x,i}) ← NIM.Encode_B(crs, x_i)
    let mut nim = Nim::new(crs.clone());
    let (pe_x_i, st_x_i) = nim.encode_B(&Mpz::from(my_share.clone())).unwrap();
    
    // 2. Compute a ZKPoK that pe_{x,i} and X_i use the same x_i
    let thpp = setup(128, pp.n, pp.t);
    let nizk = NIZKAoK::new(&thpp);
    let mut rng = RandGen::new();
    let proof = nizk.prove_cl(&thpp, &mut rng, &pe_x_i, &pub_shares[&my_id], &st_x_i, &my_share);
    
    // 3. Broadcast pe_{x,i} and the ZKPoK
    let aux_info_msg = AuxInfoMsg {
        commitment: pe_x_i.clone(),
        proof,
    };
    
    outgoing
        .send(Outgoing::broadcast(DkgMsg::AuxInfo(aux_info_msg.clone())))
        .await
        .map_err(Error::Round3Send)?;
    
    let aux_info_messages = rounds.complete(round3).await.map_err(Error::Round3Recv)?;
    
    // 4. Set aux := ({pe_{x,i}}_{i∈[m]})
    let mut aux = BTreeMap::new();
    aux.insert(my_id, pe_x_i); // Add my commitment with my ID
    
    aux_info_messages
        .into_iter_indexed()
        .map(|(inner_id, _, msg)| ((inner_id + 1) as Id, msg))
        .filter(|(id, msg)| {
            // Verify the ZKPoK for each message
            nizk.verify_cl(&thpp, &msg.commitment, &pub_shares[id], &msg.proof)
        })
        .for_each(|(id, msg)| {
            aux.insert(id, msg.commitment);
        });

    Ok((pk, pub_shares, my_share, aux, st_x_i))
}


pub fn keygen_step_1(
    pp: &PublicParams,
) -> (G, Zq, BTreeMap<Id, Zq>, BTreeMap<Id, G>) {
    // 并行生成随机系数
    let coeffs: Vec<_> = (0..pp.t)
        .into_par_iter()
        .map(|_| Zq::random())
        .collect();
    
    let x = coeffs[0].clone();
    let polynomial = Polynomial { coeffs };
    
    let curve_polynomial = CurvePolynomial::from_exp(&polynomial, &pp.g);
    
    let (pub_shares, x_shares): (BTreeMap<_, _>, BTreeMap<_, _>) = (1..=pp.n)
        .into_par_iter()
        .map(|i| {
            let i_scalar = Zq::from(i as u64);
            let pub_share = curve_polynomial.eval(&i_scalar);
            let x_share = polynomial.eval(&i_scalar);
            ((i, pub_share), (i, x_share))
        })
        .unzip();

    let threshold_pk = ThresholdPubKey {
        pk: curve_polynomial.coeffs[0].clone(),
        pub_shares,
    };

    (threshold_pk.pk, x, x_shares, threshold_pk.pub_shares)
}

pub async fn keygen_step_2<M>(
    party: M,
    my_id: Id,
    pp: &PublicParams,
    pub_shares_map: &BTreeMap<Id, G>,
    my_share: &Zq,
) -> Result<(BTreeMap<Id, bicycl::CipherText>, Mpz), Error<M::ReceiveError, M::SendError>>
where
    M: Mpc<ProtocolMessage = DkgMsg>,
{
    let mut nim = Nim::new(pp.crs.clone());
    let (pe_x_i, st_x_i) = nim.encode_B(&Mpz::from(my_share.clone())).unwrap();

    let nizk = NIZKAoK::new(pp);
    let mut rng = RandGen::new();
    let proof = nizk.prove_cl(&pp, &mut rng, &pe_x_i, &pub_shares_map[&my_id], &st_x_i, &my_share);
    let aux_info_msg = AuxInfoMsg {
        commitment: pe_x_i.clone(),
        proof
    };

    let MpcParty { delivery, .. } = party.into_party();
    let (incoming, mut outgoing) = delivery.split();

    let i = (my_id - 1) as u16;
    let n = pp.n as u16;

    let mut rounds = RoundsRouter::<DkgMsg>::builder();
    let round1 = rounds.add_round(RoundInput::<AuxInfoMsg>::broadcast(i, n));
    let mut rounds = rounds.listen(incoming);
    
    outgoing
    .send(Outgoing::broadcast(DkgMsg::AuxInfo(aux_info_msg.clone())))
    .await
    .map_err(Error::Round3Send)?;

    let aux_info_messages = rounds.complete(round1).await.map_err(Error::Round3Recv)?;

    let mut aux = BTreeMap::new();

    aux_info_messages
    .into_iter_indexed()
    .map(|(inner_id, _, msg)| ((inner_id + 1) as Id, msg))
    .filter(|(id, msg)| {
        // Verify the ZKPoK for each message
        nizk.verify_cl(&pp, &msg.commitment, &pub_shares_map[id], &msg.proof)
    })
    .for_each(|(id, msg)| {
        aux.insert(id, msg.commitment);
    });

    aux.insert(my_id, pe_x_i); // Add my commitment with my ID

    Ok((aux, st_x_i))
}

pub async fn keygen_simulation(
    pp: &PublicParams,
) -> Result<(Zq, G, BTreeMap<Id, G>, BTreeMap<Id, bicycl::CipherText>, BTreeMap<Id, Zq>,  BTreeMap<Id, Mpz>), Error<round_based::rounds_router::simple_store::RoundInputError, round_based::rounds_router::simple_store::RoundInputError>>
{
    let (pk, x, share_map, pub_shares_map) = keygen_step_1(pp);
    
    let mut simulation = Simulation::<DkgMsg>::new();
    let mut party_output = vec![];

    for i in 1..=pp.n {
        let party = simulation.add_party();
        let result = keygen_step_2(party, i as Id, &pp, &pub_shares_map, &share_map[&i]);
        party_output.push(result);
    }

    let res = futures::future::try_join_all(party_output).await.unwrap();

    let aux_map = res[0].0.clone();
    let mut st_x_i_map = BTreeMap::new();

    for (i, (_aux, st_x_i)) in res.iter().enumerate() {
        // if _aux != &aux_map {
        //     return Err(Error::Step2Error("Auxiliary information mismatch".to_string()));
        // }
        st_x_i_map.insert((i + 1) as Id, st_x_i.clone());
    }

    Ok((x, pk, pub_shares_map, aux_map, share_map, st_x_i_map))
}


#[cfg(test)]
mod tests {
    use std::time::Instant;
    use super::*;
    use curv::arithmetic::Converter;

    const N: u8 = 3;
    const T: u8 = 2;
    #[tokio::test]
    async fn test_keygen() {
        let pp = setup(128, N, T);
        let (x, pk, pub_shares_map, aux_map, share_map, st_x_i_map) = keygen_simulation(&pp).await.unwrap();

        println!("\n=== keygen test result ===\n");
        println!("public key (pk):");
        println!("  {}", pk.to_bytes(false).iter().map(|b| format!("{:02x}", b)).collect::<String>());
        
        println!("main secret key (x):");
        println!("  {}", x.to_bigint().to_str_radix(16));
        
        println!("share map (share_map):");
        for (id, share) in share_map.iter() {
            println!("  party {}: {}", id, share.to_bigint().to_str_radix(16));
        }
        println!();
        
        println!("public key shares (pub_shares):");
        for (id, pub_share) in pub_shares_map.iter() {
            println!("  party {}: {}", id, pub_share.to_bytes(false).iter().map(|b| format!("{:02x}", b)).collect::<String>());
        }
        println!("\n========================\n");    

        println!("auxiliary information (aux):");
        for (id, ciphertext) in aux_map.iter() {
            println!("    pe_x{} (CipherText):", id);
            println!("      c1（QFI）: {:?}", ciphertext.c1());
            println!("      c2（QFI）: {:?}", ciphertext.c2());
        }
        println!("\n========================\n");
        
        println!("st_x_i_map:");
        for (id, st_x_i) in st_x_i_map.iter() {
            println!("  party {}: {}", id, st_x_i.to_bytes().iter().map(|b| format!("{:02x}", b)).collect::<String>());
        }
        println!("\n========================\n");
        
    }

    #[tokio::test]
    async fn test_keygen_step2_time() {
        let pp = setup(128, N, T);
        let (_pk, _x, share_map, pub_shares_map) = keygen_step_1(&pp);


        let mut simulation = Simulation::<DkgMsg>::new();
        let mut party_output = vec![];

        let start = Instant::now();

        for i in 1..=pp.n {
            let party = simulation.add_party();
            let result = keygen_step_2(party, i as Id, &pp, &pub_shares_map, &share_map[&i]);
            party_output.push(result);
        }

        let _res = futures::future::try_join_all(party_output).await.unwrap();

        let duration = start.elapsed();

        println!("threshold ({}, {}), broadcast phase time: {:?}", N, T, duration);
    }

    #[tokio::test]
    async fn test_dkg() {
        let _pp = setup(128, N, T);
        let crs = _pp.crs;

        let (pp, secret_keys) = simulate_pp(N, T);
        let h = G::base_point2();
        let mut simulation = Simulation::<DkgMsg>::new();
        let mut party_output = vec![];

        let start = Instant::now();

        for i in 1..=pp.n {
            let party = simulation.add_party();
            let result = dkg(party, i, &pp, h, &secret_keys[&i], false, &crs);
            party_output.push(result);
        }

        let res = futures::future::try_join_all(party_output).await.unwrap();

        let duration = start.elapsed();

        for (i, (pk, pub_shares, my_share, aux, st_x_i)) in res.iter().enumerate() {
            println!("party {} result:", i + 1);
            println!("  pk:");
            println!("    {}", pk.to_bytes(false).iter().map(|b| format!("{:02x}", b)).collect::<String>());
            println!("  pub_shares:");
            for (id, share) in pub_shares.iter() {
                println!("    {}", id);
                println!("    {}", share.to_bytes(false).iter().map(|b| format!("{:02x}", b)).collect::<String>());
            }
            println!("  my_share:");
            println!("    {}", my_share.to_bigint().to_str_radix(16));
            println!("  aux:");
            for (id, ciphertext) in aux.iter() {
                println!("    pe_x{} (CipherText):", id);
                println!("      c1（QFI）: {:?}", ciphertext.c1());
                println!("      c2（QFI）: {:?}", ciphertext.c2());
            }
            println!("  st_x{} (Mpz):", i + 1);
            println!("    {}", st_x_i.to_bytes().iter().map(|b| format!("{:02x}", b)).collect::<String>());
            println!();
        }


        println!("threshold ({}, {}), broadcast phase time: {:?}", N, T, duration);
    }
}
