use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use rand::random;
use ecdsa_2r::setup::*;
use ecdsa_2r::keygen::*;
use ecdsa_2r::presign::*;
use ecdsa_2r::signature::*;

const N: u8 = 7;
const T: u8 = 5;
const MSG: &[u8] = b"Hello, world!";

fn bench_setup(c: &mut Criterion) {
    c.bench_function("Setup", |b| {
        b.iter(|| {
            setup(128, N, T)
        })
    });
}

fn bench_keygen(c: &mut Criterion) {
    let rt = tokio::runtime::Runtime::new().unwrap();
    let pp = setup(128, N, T);
    
    c.bench_function("KeyGen", |b| {
        b.iter(|| {
            rt.block_on(async {
                keygen_simulation(&pp).await.unwrap()
            })
        })
    });
}

fn bench_presign(c: &mut Criterion) {
    let rt = tokio::runtime::Runtime::new().unwrap();
    let pp = setup(128, N, T);
    let participant = (1..=N).collect::<Vec<_>>();
    
    let (_x, _pk, _pub_shares_map, aux_map, share_map, st_x_i_map) = rt.block_on(async {
        keygen_simulation(&pp).await.unwrap()
    });
    
    c.bench_function("PreSign", |b| {
        b.iter(|| {
            rt.block_on(async {
                presign_simulation(&participant, &pp, &share_map, &st_x_i_map, &aux_map).await
            })
        })
    });
}

fn bench_sign(c: &mut Criterion) {
    let rt = tokio::runtime::Runtime::new().unwrap();
    let pp = setup(128, N, T);
    let participant = (1..=N).collect::<Vec<_>>();
    
    let (_x, pk, _pub_shares_map, aux_map, share_map, st_x_i_map) = rt.block_on(async {
        keygen_simulation(&pp).await.unwrap()
    });
    
    c.bench_function("Signing", |b| {
        b.iter_with_setup(
            || {
                let presign_output = rt.block_on(async {
                    presign_simulation(&participant, &pp, &share_map, &st_x_i_map, &aux_map).await
                });
                presign_output
            },
            |presign| {
                let (w_i_list, u_i_list, r_list) = parallel_partial_sign(presign, &pp, &pk, &MSG);
                combine_sign(w_i_list, u_i_list, r_list)
            }
        )
    });
}

fn bench_verify(c: &mut Criterion) {
    let rt = tokio::runtime::Runtime::new().unwrap();
    let pp = setup(128, N, T);
    let participant = (1..=N).collect::<Vec<_>>();
    
    let (_x, pk, _pub_shares_map, aux_map, share_map, st_x_i_map) = rt.block_on(async {
        keygen_simulation(&pp).await.unwrap()
    });
    
    c.bench_function("Verify", |b| {
        b.iter_with_setup(|| {
            let presign_output = rt.block_on(async {
                presign_simulation(&participant, &pp, &share_map, &st_x_i_map, &aux_map).await
            });
            let (w_i_list, u_i_list, r_list) = parallel_partial_sign(presign_output, &pp, &pk, &MSG);
            let signature = combine_sign(w_i_list, u_i_list, r_list);
            signature
        }, |signature| {
            verify_sign(&pp, &pk, &MSG, signature)
        })
    });
}

fn bench_participant_scaling(c: &mut Criterion) {
    let rt = tokio::runtime::Runtime::new().unwrap();
    let pp = setup(128, N, T);
    
    let (_x, _pk, _pub_shares_map, aux_map, share_map, st_x_i_map) = rt.block_on(async {
        keygen_simulation(&pp).await.unwrap()
    });
    
    let participant_counts = (T..=N).collect::<Vec<_>>();
    let mut group = c.benchmark_group("PreSign_Participants");
    
    for &count in &participant_counts {
        let participants = (1..=count).collect::<Vec<_>>();
        group.bench_with_input(BenchmarkId::from_parameter(count), &count, |b, &_| {
            b.iter(|| {
                rt.block_on(async {
                    presign_simulation(&participants, &pp, &share_map, &st_x_i_map, &aux_map).await
                })
            })
        });
    }
    
    group.finish();
}

fn bench_partial_sign(c: &mut Criterion) {
    let rt = tokio::runtime::Runtime::new().unwrap();
    let pp = setup(128, N, T);
    let participant = (1..=N).collect::<Vec<_>>();
    
    let (_x, pk, _pub_shares_map, aux_map, share_map, st_x_i_map) = rt.block_on(async {
        keygen_simulation(&pp).await.unwrap()
    });
    let presign_output = rt.block_on(async {
        presign_simulation(&participant, &pp, &share_map, &st_x_i_map, &aux_map).await
    });
    
    c.bench_function("PartialSign", |b| {
        b.iter_with_setup(
            || {
                let i = random::<usize>() % presign_output.len();
                let (g_k, pw, pu, gamma_i, presign_m_list) = presign_output[i].clone();
                (g_k, pw, pu, gamma_i, presign_m_list)
            },
            |(g_k, pw, pu, gamma_i, presign_m_list)| {
                partial_sign(&pp, &pk, &MSG, &g_k, presign_m_list, pw, pu, gamma_i);
            }
        )
    });
}

fn bench_combine_sign_scaling(c: &mut Criterion) {
    let rt = tokio::runtime::Runtime::new().unwrap();
    let pp = setup(128, N, T);
    
    let (_x, pk, _pub_shares_map, aux_map, share_map, st_x_i_map) = rt.block_on(async {
        keygen_simulation(&pp).await.unwrap()
    });

    let participant_counts = (T..=N).collect::<Vec<_>>();
    let mut group = c.benchmark_group("CombineSign_Participants");

    for &count in &participant_counts {
        let participant = (1..=count).collect::<Vec<_>>();
        group.bench_with_input(BenchmarkId::from_parameter(count), &count, |b, &_| {
            b.iter_with_setup(
                || {
                    let presign_output = rt.block_on(async {
                        presign_simulation(&participant, &pp, &share_map, &st_x_i_map, &aux_map).await
                    });
                    let (w_i_list, u_i_list, r_list) = parallel_partial_sign(presign_output, &pp, &pk, &MSG);
                    (w_i_list, u_i_list, r_list)
                },
                |(w_i_list, u_i_list, r_list)| {
                    combine_sign(w_i_list, u_i_list, r_list)
                }
            )
        });
    }
}

fn bench_full_pipeline(c: &mut Criterion) {
    let rt = tokio::runtime::Runtime::new().unwrap();
    let participant = (1..=N).collect::<Vec<_>>();
    
    c.bench_function("Full_Pipeline", |b| {
        b.iter(|| {
            rt.block_on(async {
                // 1. Setup
                let pp = setup(128, N, T);
                
                // 2. Key Generation
                let (_x, pk, _pub_shares_map, aux_map, share_map, st_x_i_map) = keygen_simulation(&pp).await.unwrap();
                
                // 3. PreSign
                let output = presign_simulation(&participant, &pp, &share_map, &st_x_i_map, &aux_map).await;
                
                // 4. Sign
                let (w_i_list, u_i_list, r_list) = parallel_partial_sign(output, &pp, &pk, &MSG);
                let signature = combine_sign(w_i_list, u_i_list, r_list);
                
                // 5. Verify
                verify_sign(&pp, &pk, &MSG, signature)
            })
        })
    });
}

fn bench_nizk(c: &mut Criterion) {
    use ecdsa_2r::nim::*;
    use ecdsa_2r::zk::*;
    use ecdsa_2r::utils::Zq;
    use bicycl::{Mpz, RandGen};
    
    let pp = setup(128, N, T);
    let mut nim = Nim::new(pp.crs.clone());
    let mut rnd = RandGen::new();
    
    let x = Zq::random();
    let y = Zq::random();
    let x_mpz = Mpz::from(&x);
    let y_mpz = Mpz::from(&y);
    
    let (pe_x, st_x) = nim.encode_A(&x_mpz).unwrap();
    let (pe_y, st_y) = nim.encode_B(&y_mpz).unwrap();
    
    let nizk = NIZKAoK::new(&pp);
    let g_x = &pp.g * &x;
    let g_y = &pp.g * &y;
    
    c.bench_function("NIZK_CL_Prove", |b| {
        b.iter_with_setup(|| {
            let pe_y = pe_y.clone();
            let g_y = g_y.clone();
            let st_y = st_y.clone();
            let y = y.clone();
            (pe_y, g_y, st_y, y)
        }, |(pe_y, g_y, st_y, y)| {
            nizk.prove_cl(&pp, &mut rnd, &pe_y, &g_y, &st_y, &y)
        })
    });
    
    let proof_cl = nizk.prove_cl(&pp, &mut rnd, &pe_y, &g_y, &st_y, &y);
    
    c.bench_function("NIZK_CL_Verify", |b| {
        b.iter_with_setup(|| {
            let pe_y = pe_y.clone();
            let g_y = g_y.clone();
            let proof_cl = proof_cl.clone();
            (pe_y, g_y, proof_cl)
        }, |(pe_y, g_y, proof_cl)| {
            nizk.verify_cl(&pp, &pe_y, &g_y, &proof_cl)
        })
    });
    
    c.bench_function("NIZK_Pedersen_Prove", |b| {
        b.iter_with_setup(|| {
            let pe_x = pe_x.clone();
            let g_x = g_x.clone();
            let st_x = st_x.clone();
            let x = x.clone();
            (pe_x, g_x, st_x, x)
        }, |(pe_x, g_x, st_x, x)| {
            nizk.prove_ped(&pp, &mut rnd, &pe_x, &g_x, &st_x, &x)
        })
    });
    
    let proof_ped = nizk.prove_ped(&pp, &mut rnd, &pe_x, &g_x, &st_x, &x);
    
    c.bench_function("NIZK_Pedersen_Verify", |b| {
        b.iter_with_setup(|| {
            let pe_x = pe_x.clone();
            let g_x = g_x.clone();
            let proof_ped = proof_ped.clone();
            (pe_x, g_x, proof_ped)
        }, |(pe_x, g_x, proof_ped)| {
            nizk.verify_ped(&pp, &pe_x, &g_x, &proof_ped)
        })
    });
}

criterion_group!(
    benches,
    bench_setup,
    bench_keygen,
    bench_presign,
    bench_sign,
    bench_verify,
    bench_participant_scaling,
    bench_partial_sign,
    bench_combine_sign_scaling,
    bench_full_pipeline,
    bench_nizk
);
criterion_main!(benches); 