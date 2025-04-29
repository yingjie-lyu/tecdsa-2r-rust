use bicycl::{CipherText, Mpz, QFI, RandGen};
use curv::{arithmetic::Converter, BigInt};
use futures::SinkExt;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

use crate::utils::*;
use crate::nim::Nim;
use round_based::
rounds_router::{
        simple_store::{RoundInput, RoundInputError},
        CompleteRoundError, RoundsRouter,
    };
use round_based::simulation::Simulation;
use round_based::{Delivery, Mpc, MpcParty, Outgoing, ProtocolMessage};
use thiserror::Error;
use crate::setup::*;
use crate::lagrange::*;
use crate::zk::*;

#[derive(Debug, Error)]
pub enum Error<RecvErr, SendErr> {
    #[error("sending, round 1")]
    Round1Send(#[source] SendErr),
    #[error("receiving, round 1")]
    Round1Recv(#[source] CompleteRoundError<RoundInputError, RecvErr>),
    #[error("step1 error")]
    Step1Error(String),
    #[error("step2 error")]
    Step2Error(String),
    #[error("insufficient valid messages")]
    InsufficientMessages,
}

// protocol message M_i^(1) := (K_i, Γ_i, pe_{k,i}, pe_{γ,i}, π_i)
#[derive(Clone, Debug, PartialEq, ProtocolMessage, Serialize, Deserialize)]
pub enum PresignMsg {
    Round1(PreSignM),
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PreSignM {
    pub K: G,          // K_i := g^k_i
    pub Gamma: G,      // Γ_i := g^γ_i
    pub pe_k: CipherText,     // pe_{k,i} 
    pub pe_gamma: QFI, // pe_{γ,i}
    pub ped_proof: PedProof,
    pub cl_proof: CLProof,
}

impl PreSignM {
    // convert PreSignM to bytes
    pub fn to_bytes(&self) -> Vec<u8> {
        // convert K and Gamma to bytes, using compressed format
        let k_bytes = self.K.to_bytes(true).to_vec();
        let gamma_bytes = self.Gamma.to_bytes(true).to_vec();
        
        // convert pe_k and pe_gamma to bytes
        let pe_k_bytes = self.pe_k.to_bytes();
        let pe_gamma_bytes = self.pe_gamma.to_bytes();
        
        // TODO: add ped_proof and cl_proof to bytes
        let ped_proof_bytes = self.ped_proof.to_bytes();
        let cl_proof_bytes = self.cl_proof.to_bytes();

        // concatenate all bytes together
        [k_bytes, gamma_bytes, pe_k_bytes, pe_gamma_bytes, ped_proof_bytes, cl_proof_bytes].concat()
    }
}

// secret state S_i := (k_i, γ_i, st_{k,i}, st_{γ,i})
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PreSignS {
    k: Mpz,          // k_i
    gamma: Mpz,      // γ_i
    st_k: Mpz,       // st_{k,i}
    st_gamma: Mpz,   // st_{γ,i}
}

// add function to describe the protocol
pub fn presign_step1(
    pp: &PublicParams,
) -> Result<(PreSignM, PreSignS), Box<dyn std::error::Error>> {
    // 1. sample k_i and γ_i from Zq
    let k = Zq::random();
    let gamma = Zq::random();

    // 2. calculate K_i := g^k_i and Γ_i := g^γ_i
    let K = G::generator() * &k;
    let Gamma = G::generator() * &gamma;

    // get Mpz representation of k and gamma
    let k_mpz = Mpz::from_bytes(&k.to_bytes());
    let gamma_mpz = Mpz::from_bytes(&gamma.to_bytes());

    // use NIM encoding - according to the protocol use crs
    let mut nim = Nim::new(pp.crs.clone());
    let (pe_k, st_k) = nim.encode_B(&k_mpz).unwrap();
    let (pe_gamma, st_gamma) = nim.encode_A(&gamma_mpz).unwrap();

    // 3. calculate zero-knowledge proof π_i (not implemented yet)
    // here should prove pe_k and K use the same k, pe_gamma and Gamma use the same gamma
    let nizk = NIZKAoK::new(&pp);
    let mut rng = RandGen::new();
    let proof_k = nizk.prove_cl(&pp, &mut rng, &pe_k, &K, &st_k, &k);
    let proof_gamma = nizk.prove_ped(&pp, &mut rng, &pe_gamma, &Gamma, &st_gamma, &gamma);

    // 4. output protocol message and secret state
    let message = PreSignM {
        K,
        Gamma,
        pe_k,
        pe_gamma,
        ped_proof: proof_gamma,
        cl_proof: proof_k,
    };

    let secret = PreSignS {
        k: k_mpz,
        gamma: gamma_mpz,
        st_k,
        st_gamma,
    };

    Ok((message, secret))
}


pub fn presign_step2(
    pp: &PublicParams,
    my_id: Id,
    my_share: &Zq,
    my_secret: &PreSignS,
    st_x: &Mpz,
    pe_x_j_map: &BTreeMap<Id, CipherText>,
    valid_messages: &BTreeMap<Id, PreSignM>,
) -> Result<(G, Zq, Zq), Box<dyn std::error::Error>> {
    // calculate Lagrange coefficients λ
    let parties: Vec<Id> = valid_messages.keys().cloned().collect();
    let lagrange_coefficients = compute_lagrange_coefficients(&parties)?;
    let my_lambda = lagrange_coefficients.get(&my_id).unwrap();
    
    // calculate K as the product of all K_j
    let mut K = G::generator() * &Zq::from(0u64);
    for (_, msg_j) in valid_messages.iter() {
        K = K + &msg_j.K;
    }
    
    let nim = Nim::new(pp.crs.clone());
    // process each j ∈ P\{i}
    let mut pw = Zq::from(0u64);
    let mut pu = Zq::from(0u64);
    for (j, msg_j) in valid_messages.iter() {
        if *j == my_id {
            continue;
        }

        // calculate α_{i,j}、β_{j,i}
        // α_{i,j} := NIM.Decode_B(crs, pe_{γ,j}, st_{k,i})
        let alpha_i_j = nim.decode_B(&msg_j.pe_gamma, &my_secret.st_k)
            .map_err(|e| Box::new(std::io::Error::new(
                std::io::ErrorKind::Other, 
                format!("Failed to decode alpha: {}", e)
            )))?;
        
        let alpha_i_j_mpz = pp.crs.cl.dlog_in_F(&alpha_i_j);
        let alpha_i_j_zq = Zq::from_bigint(&BigInt::from_bytes(&alpha_i_j_mpz.to_bytes()));
        pu = pu + alpha_i_j_zq;
        
        // β_{j,i} := NIM.Decode_A(crs, pe_{k,j}, st_{γ,i})
        let beta_j_i = nim.decode_A(&msg_j.pe_k, &my_secret.st_gamma, &my_secret.gamma)
            .map_err(|e| Box::new(std::io::Error::new(
                std::io::ErrorKind::Other, 
                format!("Failed to decode beta: {}", e)
            )))?;

        let beta_j_i_mpz = pp.crs.cl.dlog_in_F(&beta_j_i);  
        let beta_j_i_zq = Zq::from_bigint(&BigInt::from_bytes(&beta_j_i_mpz.to_bytes()));
        pu = pu + beta_j_i_zq;
        
        // μ_{i,j} := λ_{i,P} · NIM.Decode_B(crs, pe_{γ,j}, st_{x,i})
        let mu_decoded = nim.decode_B(&msg_j.pe_gamma, st_x)
            .map_err(|e| Box::new(std::io::Error::new(
                std::io::ErrorKind::Other, 
                format!("Failed to decode mu term: {}", e)
            )))?;
        
        let mu_decoded_mpz = pp.crs.cl.dlog_in_F(&mu_decoded);
        let mu_decoded_zq = Zq::from_bigint(&BigInt::from_bytes(&mu_decoded_mpz.to_bytes()));
        pw = pw + mu_decoded_zq * my_lambda;
        
        // ν_{j,i} := λ_{j,P} · NIM.Decode_A(crs, pe_{x,j}, st_{γ,i})
        let pe_x_j = pe_x_j_map.get(j).unwrap();
        let nu_decoded = nim.decode_A(pe_x_j, &my_secret.st_gamma,&my_secret.gamma)
            .map_err(|e| Box::new(std::io::Error::new(
                std::io::ErrorKind::Other, 
                format!("Failed to decode nu term: {}", e)
            )))?;
        
        let nu_decoded_mpz = pp.crs.cl.dlog_in_F(&nu_decoded);
        let nu_decoded_zq = Zq::from_bigint(&BigInt::from_bytes(&nu_decoded_mpz.to_bytes()));
        let lambda_j = lagrange_coefficients.get(j).unwrap();
        pw = pw + nu_decoded_zq * lambda_j;
    }

    let k_Zq = Zq::from_bigint(&BigInt::from_bytes(&my_secret.k.to_bytes()));
    let gamma_Zq = Zq::from_bigint(&BigInt::from_bytes(&my_secret.gamma.to_bytes()));

    pw = pw + my_lambda * my_share * &gamma_Zq;
    pu = pu + k_Zq * gamma_Zq;

    Ok((K, pw, pu))
}

pub async fn protocol<M>(
    party: M,
    party_num: u16,
    my_id: Id,
    pp: &PublicParams,
    my_share: &Zq,
    st_x: &Mpz,
    pe_x_j_map: &BTreeMap<Id, CipherText>,
) -> Result<(G, Zq, Zq, Zq, BTreeMap<Id, PreSignM>), Error<M::ReceiveError, M::SendError>>
where
    M: Mpc<ProtocolMessage = PresignMsg>,
{
    let (my_prsign_m, my_presign_s) = presign_step1(pp).map_err(|e| Error::Step1Error(e.to_string()))?;

    let MpcParty { delivery, .. } = party.into_party();
    let (incoming, mut outgoing) = delivery.split();

    let i = (my_id - 1) as u16;

    outgoing
        .send(Outgoing::broadcast(PresignMsg::Round1(my_prsign_m.clone())))
        .await
        .map_err(Error::Round1Send)?;

    let mut rounds = RoundsRouter::<PresignMsg>::builder();
    let round1 = rounds.add_round(RoundInput::<PreSignM>::broadcast(i, party_num));
    let mut rounds = rounds.listen(incoming);

    let presign_messages = rounds.complete(round1).await.map_err(Error::Round1Recv)?;

    let mut valid_presign_messages = BTreeMap::new();

    let nizk = NIZKAoK::new(&pp);
    presign_messages
        .into_iter_indexed()
        .map(|(inner_id, _, msg)| ((inner_id + 1) as Id, msg))
        .filter(|(_, msg)| {
            nizk.verify_cl(&pp, &msg.pe_k, &msg.K, &msg.cl_proof) && nizk.verify_ped(&pp, &msg.pe_gamma, &msg.Gamma, &msg.ped_proof)
        })
        .for_each(|(j, msg)| {
            valid_presign_messages.insert(j, msg);
        });
    
    valid_presign_messages.insert(my_id, my_prsign_m);   
    // check if there are enough valid messages (at least t)
    if valid_presign_messages.len() < pp.t as usize {
        return Err(Error::InsufficientMessages);
    }

    // signature pre-calculation step2, calculate K, α_{i,j}, β_{j,i}, μ_{i,j}, ν_{j,i}, λ_{i,P} * x_{i} * γ_{i}, k_{i} * γ_{i}
    let (K, pw, pu) = presign_step2(
        pp, my_id, my_share, &my_presign_s, st_x, pe_x_j_map, &valid_presign_messages)
        .map_err(|e| Error::Step2Error(e.to_string()))?;
    let gamma_i = Zq::from_bigint(&BigInt::from_bytes(&my_presign_s.gamma.to_bytes()));

    Ok((K, pw, pu, gamma_i, valid_presign_messages))
}

pub async fn presign_simulation(
    participant: &Vec<Id>,
    pp: &PublicParams,
    share_map: &BTreeMap<Id, Zq>,
    st_x_i_map: &BTreeMap<Id, Mpz>,
    aux_map: &BTreeMap<Id, CipherText>,
) -> Vec<(G, Zq, Zq, Zq, BTreeMap<Id, PreSignM>)>
{
    let mut simulation = Simulation::<PresignMsg>::new();

    let mut party_output = vec![];
    for (_, id) in participant.iter().enumerate() {
        let party = simulation.add_party(); 
        let result = protocol(
            party,
            participant.len() as u16,
            *id,
            &pp,
            &share_map[id],
            &st_x_i_map[id],
            &aux_map,
        );
        party_output.push(result);
    }

    futures::future::try_join_all(party_output).await.unwrap()
}

// get a random set of t participants
pub fn get_random_participant_set(t: u8, n: u8) -> Vec<Id> {
    use rand::seq::SliceRandom;
    
    let mut rng = rand::thread_rng();
    let mut participant_set = (1..=n).collect::<Vec<_>>();
    participant_set.shuffle(&mut rng);
    participant_set[0..(t as usize)].to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;
    use crate::setup::setup;
    use crate::keygen::*;

    const N: u8 = 3;
    const T: u8 = 2;

    #[test]
    fn test_get_random_participant_set() {
        let t = 3;
        let n = 5;
        
        // test multiple times to ensure the function has randomness
        for _ in 0..10 {
            let participants = get_random_participant_set(t, n);
            
            // verify the result length is correct
            assert_eq!(participants.len(), t as usize, "the number of participants should be t");
            
            // verify all IDs are in the range of 1 to n
            for &id in &participants {
                assert!(id >= 1 && id <= n, "participant ID should be in the range of 1 to n");
            }
            
            // verify no duplicate IDs
            let mut unique_ids = participants.clone();
            unique_ids.sort();
            unique_ids.dedup();
            assert_eq!(unique_ids.len(), participants.len(), "participant IDs should not be duplicated");
        }
        
        println!("random participant set test passed");
    }

    #[test]
    fn test_presign_step1() {
        let pp = setup(128, N, T);
        
        println!("测试presign_step1...");
        let now = Instant::now();
        
        let result = presign_step1(&pp);
        
        // verify the result
        assert!(result.is_ok(), "presign_step1 should successfully return the result");
        
        let (message, secret) = result.unwrap();
        
        // verify K and Gamma are generated correctly from k and gamma
        let k_bytes = secret.k.to_bytes();
        let gamma_bytes = secret.gamma.to_bytes();
        let k_zq = Zq::from_bytes(&k_bytes).expect("failed to convert k to Zq");
        let gamma_zq = Zq::from_bytes(&gamma_bytes).expect("failed to convert gamma to Zq");
        
        let expected_K = G::generator() * &k_zq;
        let expected_Gamma = G::generator() * &gamma_zq;
        
        assert_eq!(message.K, expected_K, "K verification failed");
        assert_eq!(message.Gamma, expected_Gamma, "Gamma verification failed");
        
        let elapsed = now.elapsed();

        println!("test passed! time: {:.2?}", elapsed);
    }

    #[tokio::test]
    async fn test_presign() {
        let pp = setup(128, N, T);

        let (_, _, _, aux_map, share_map, st_x_i_map) = keygen_simulation(&pp).await.unwrap();

        let participant = vec![1, 2];
        let now = Instant::now();
        let output = presign_simulation(&participant, &pp, &share_map, &st_x_i_map, &aux_map).await;
        let elapsed = now.elapsed();

        for (i, result) in output.into_iter().enumerate() {
            let (K, pw, pu, gamma_i, _) = result;
            println!("participant {} result: ", i + 1);
            println!("K: {:?}", K.to_bytes(false).iter().map(|b| format!("{:02x}", b)).collect::<String>());
            println!("pw: {:?}", pw.to_bigint().to_str_radix(16));
            println!("pu: {:?}", pu.to_bigint().to_str_radix(16));
            println!("γ_{}: {:?}", i + 1, gamma_i.to_bigint().to_str_radix(16));

            println!("γ_{}: {:?}", i + 1, gamma_i.to_bytes().iter().map(|b| format!("{:02x}", b)).collect::<String>());
        }

        println!("test passed! time: {:.2?}", elapsed);
    }

}