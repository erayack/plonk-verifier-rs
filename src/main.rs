use curve25519_dalek::ristretto::CompressedRistretto;
use curve25519_dalek::traits::Identity;
use curve25519_dalek::scalar::Scalar;

fn verify_proof(setup: &str, group_order: u64, vk: &VerificationKey, proof: &Proof, public: Vec<Scalar>, optimized: bool) -> Result<(), &'static str> {
    let (A_pt, B_pt, C_pt, Z_pt, T1_pt, T2_pt, T3_pt, W_z_pt, W_zw_pt, A_ev, B_ev, C_ev, S1_ev, S2_ev, Z_shifted_ev) = proof;

    let Ql_pt = &vk.Ql;
    let Qr_pt = &vk.Qr;
    let Qm_pt = &vk.Qm;
    let Qo_pt = &vk.Qo;
    let Qc_pt = &vk.Qc;
    let S1_pt = &vk.S1;
    let S2_pt = &vk.S2;
    let S3_pt = &vk.S3;
    let X2 = &vk.X_2;

    // Compute challenges (should be same as those computed by prover)
    let mut buf = Vec::new();
    buf.extend_from_slice(A_pt.compress().as_bytes());
    buf.extend_from_slice(B_pt.compress().as_bytes());
    buf.extend_from_slice(C_pt.compress().as_bytes());

    let beta = binhash_to_f_inner(keccak256(&buf));
    let gamma = binhash_to_f_inner(keccak256(keccak256(&buf)));

    let alpha = binhash_to_f_inner(keccak256(Z_pt.compress().as_bytes()));

    let mut buf2 = Vec::new();
    buf2.extend_from_slice(T1_pt.compress().as_bytes());
    buf2.extend_from_slice(T2_pt.compress().as_bytes());
    buf2.extend_from_slice(T3_pt.compress().as_bytes());
    let zed = binhash_to_f_inner(keccak256(&buf2));

    let mut buf3 = Vec::new();
    buf3.extend_from_slice(&A_ev.as_bytes());
    buf3.extend_from_slice(&B_ev.as_bytes());
    buf3.extend_from_slice(&C_ev.as_bytes());
    buf3.extend_from_slice(&S1_ev.as_bytes());
    buf3.extend_from_slice(&S2_ev.as_bytes());
    buf3.extend_from_slice(&Z_shifted_ev.as_bytes());

    fn binhash_to_f_inner(hash: &str) -> f64 {
        // Convert the binary hash value to some other value
        // and return the result
    }
    
    fn keccak256(buf: &[u8]) -> String {
        // Implement the Keccak-256 hash function and return the
        // resulting hash value as a string
    }
    
    fn main() {
        let buf = [0u8, 1, 2, 3];
        let buf2 = [4u8, 5, 6, 7];
        let buf3 = [8u8, 9, 10, 11];
    
        let v = binhash_to_f_inner(&keccak256(&buf3));
    
        // Does not need to be standardized, only needs to be unpredictable
        let u = binhash_to_f_inner(&keccak256(&[&buf, &buf2, &buf3].concat()));
    
        let ZH_ev = (zed.pow(group_order) - 1.0) as f64;
    
        let root_of_unity = get_root_of_unity(group_order);
    
        let L1_ev = (ZH_ev / (group_order * (zed - 1.0))) as f64;
    
        let PI_ev = barycentric_eval_at_point(
            [
                f_inner(-x) for x in public
            ] + [
                f_inner(0) for _ in 0..group_order - public.len()
            ],
            zed
        );
    }

    if !optimized {
        // Basic, easier-to-understand version of what's going on
    
        // Recover the commitment to the linearization polynomial R,
        // exactly the same as what was created by the prover
        let R_pt = ec_lincomb(vec![
            (Qm_pt, A_ev * B_ev),
            (Ql_pt, A_ev),
            (Qr_pt, B_ev), 
            (Qo_pt, C_ev), 
            (&b.G1, PI_ev),
            (Qc_pt, 1.0),
            (Z_pt, (
                (A_ev + beta * zed + gamma) *
                (B_ev + beta * 2.0 * zed + gamma) *
                (C_ev + beta * 3.0 * zed + gamma) *
                alpha
            )),
            (S3_pt, (
                -(A_ev + beta * S1_ev + gamma) * 
                (B_ev + beta * S2_ev + gamma) *
                beta *
                alpha * Z_shifted_ev
            )),
            (&b.G1, (
                -(A_ev + beta * S1_ev + gamma) * 
                (B_ev + beta * S2_ev + gamma) *
                (C_ev + gamma) *
                alpha * Z_shifted_ev
            )),
            (Z_pt, L1_ev * alpha.powi(2)),
            (&b.G1, -L1_ev * alpha.powi(2)),
            (T1_pt, -ZH_ev),
            (T2_pt, -ZH_ev * zed.pow(group_order as u32)),
            (T3_pt, -ZH_ev * zed.pow(2 * group_order as u32)),
        ]);
    
        println!("verifier R_pt: {:?}", R_pt);
    }
