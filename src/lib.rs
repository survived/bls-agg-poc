use ff::Field;
use pairing_plus::bls12_381::{Bls12, Fq12, Fr, G1Affine as G1, G2Affine as G2, G1 as G1Proj};
use pairing_plus::hash_to_curve::HashToCurve;
use pairing_plus::hash_to_field::{hash_to_field, ExpandMsgXmd};
use pairing_plus::{CurveAffine, CurveProjective, Engine};

use rand::RngCore;

type PublicKey = G2;
type AggregatedPublicKey = G2;
type SecretKey = Fr;
type LocalSignature = G1;
type Signature = G1;

/// Domain separation tag
const TAG: &[u8] = b"BLS aggregation";

pub fn keygen<R: RngCore>(rng: &mut R) -> (PublicKey, SecretKey) {
    let sk = SecretKey::random(rng);
    let pk = G2::one().mul(sk).into_affine();
    (pk, sk)
}

pub fn key_aggregation(pk: &[PublicKey]) -> AggregatedPublicKey {
    let mut labels = LabelsGenerator::new(pk);

    let mut apk = G2::zero().into_projective();
    for pk_i in pk {
        let label_i = labels.derive_label(pk_i);
        apk.add_assign(&pk_i.mul(label_i))
    }
    apk.into_affine()
}

pub fn sign(sk_i: &SecretKey, pk: &[PublicKey], msg: &[u8]) -> LocalSignature {
    let pk_i = G2::one().mul(*sk_i).into_affine();
    let a_i = LabelsGenerator::new(pk).derive_label(&pk_i);
    let msg_on_curve =
        <G1Proj as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(msg, TAG).into_affine();

    let mut s_i = msg_on_curve.into_projective();
    s_i.mul_assign(a_i);
    s_i.mul_assign(*sk_i);
    s_i.into_affine()
}

pub fn combine(sigs: &[LocalSignature]) -> Signature {
    let mut sig = G1Proj::zero();
    for sig_i in sigs {
        sig.add_assign_mixed(sig_i);
    }
    sig.into_affine()
}

pub fn verify(
    apk: &AggregatedPublicKey,
    sig: &Signature,
    msg: &[u8],
) -> Result<(), InvalidSignature> {
    let msg_on_curve =
        <G1Proj as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(msg, TAG).into_affine();

    let a = Bls12::pairing(*sig, G2::one());
    let b = Bls12::pairing(msg_on_curve, *apk);

    if a == b {
        Ok(())
    } else {
        Err(InvalidSignature)
    }
}

pub fn verify_optimised(
    apk: &AggregatedPublicKey,
    sig: &Signature,
    msg: &[u8],
) -> Result<(), InvalidSignature> {
    // Instead of performing two pairings that involves 2 Miller loops and 2 final exponentiations,
    // we can do a well-known trick that saves 1 final exponentiation
    //
    // Unoptimised comparison: `e'(sig, g_2)^x == e'(H_0(m), apk)^x`
    // Optimised comparison  : `(e'(sig, -g_2) e'(H_0(m), apk))^x == 1`

    let msg_on_curve =
        <G1Proj as HashToCurve<ExpandMsgXmd<sha2::Sha256>>>::hash_to_curve(msg, TAG).into_affine();

    let mut neg_g_2 = G2::one();
    neg_g_2.negate();

    let result = Bls12::pairing_product(*sig, neg_g_2, msg_on_curve, *apk);
    if result == Fq12::one() {
        Ok(())
    } else {
        Err(InvalidSignature)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct InvalidSignature;

/// Derives a "label" for a party
///
/// Label is a hashed concatenation of all the public keys in form: `label = H_1(pk_i || pk_1 || ... || pk_n)`
struct LabelsGenerator {
    label_template: Vec<u8>,
}

impl LabelsGenerator {
    pub fn new(parties: &[PublicKey]) -> Self {
        let mut label_template = vec![0u8; 96 * (parties.len() + 1)];
        for (label_i, pk_i) in label_template.chunks_mut(96).skip(1).zip(parties) {
            label_i.copy_from_slice(pk_i.into_compressed().as_ref());
        }
        Self { label_template }
    }

    /// Derives label for party with public key `pk_i`
    pub fn derive_label(&mut self, pk_i: &PublicKey) -> Fr {
        self.label_template[0..96].copy_from_slice(pk_i.into_compressed().as_ref());
        hash_to_field::<Fr, ExpandMsgXmd<sha2::Sha256>>(&self.label_template, TAG, 1)[0]
    }
}

#[cfg(test)]
mod tests {
    use rand::rngs::OsRng;

    use super::*;

    #[test]
    fn sign_and_verify() {
        let (pk1, sk1) = keygen(&mut OsRng);
        let (pk2, sk2) = keygen(&mut OsRng);
        let (pk3, sk3) = keygen(&mut OsRng);
        let (pk4, sk4) = keygen(&mut OsRng);

        let pk = &[pk1, pk2, pk3, pk4];
        let apk = key_aggregation(pk);

        let msg = b"message to sign";

        let sig1 = sign(&sk1, pk, msg);
        let sig2 = sign(&sk2, pk, msg);
        let sig3 = sign(&sk3, pk, msg);
        let sig4 = sign(&sk4, pk, msg);

        let sig = combine(&[sig1, sig2, sig3, sig4]);

        verify(&apk, &sig, msg).expect("verification failed");
        verify_optimised(&apk, &sig, msg).expect("optimised verification failed");
    }
}
