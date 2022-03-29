#[test]
fn bls_signature(){
    use bellman::bls12_381::{Fr, Bls12};
    use bellman::{PrimeField, PrimeFieldRepr, Engine, CurveAffine, CurveProjective};
    use rescue_poseidon::rescue_hash;
    use rand::Rand;
    let message = "North Korea's secret nuclear bomb launch code is 69.";
    let mut repr = <Fr as PrimeField>::Repr::default();

    let bytes = message.as_bytes();
    repr.read_le(&bytes[..]).unwrap();
    let messagefr = Fr::from_repr(repr).unwrap();

    let mut rng = rand::thread_rng();
    let generator_one = <Bls12 as Engine>::G1Affine::one();
    let generator_two = <Bls12 as Engine>::G2Affine::one();
    let private_key = Fr::rand(&mut rng);
    let public_key = generator_one.mul(private_key);

    let result = rescue_hash::<Bls12, 1>(&[messagefr]);
    let hash = generator_two.mul(result[1]);

    // Private key * HASH
    let alfa = hash.into_affine().mul(private_key);

    let pair_one = hash.into_affine().pairing_with(&public_key.into_affine());
    let pair_two = generator_one.pairing_with(&alfa.into_affine());
    
    assert_eq!(pair_one, pair_two);
}
