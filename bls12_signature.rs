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

#[test]
fn aggregated_same_message(){
    use bellman::bls12_381::{Fr, Bls12};
    use bellman::{PrimeField, PrimeFieldRepr, Engine, CurveAffine, CurveProjective};
    use rescue_poseidon::rescue_hash;
    use rand::Rand;
    let message = "South Park is the best cartoon.";
    let mut repr = <Fr as PrimeField>::Repr::default();

    let bytes = message.as_bytes();
    repr.read_le(&bytes[..]).unwrap();
    let messagefr = Fr::from_repr(repr).unwrap();

    let mut rng = rand::thread_rng();
    let generator_one = <Bls12 as Engine>::G1Affine::one();
    let generator_two = <Bls12 as Engine>::G2Affine::one();
    let private_key_one = Fr::rand(&mut rng);
    let private_key_two = Fr::rand(&mut rng);
    let private_key_three = Fr::rand(&mut rng);

    let public_key_one = generator_one.mul(private_key_one);
    let public_key_two = generator_one.mul(private_key_two);
    let public_key_three = generator_one.mul(private_key_three);

    let mut aggregated_pub_keys = public_key_one;
    aggregated_pub_keys.add_assign(&public_key_two);
    aggregated_pub_keys.add_assign(&public_key_three);

    let result = rescue_hash::<Bls12, 1>(&[messagefr]);
    let hash = generator_two.mul(result[1]);
    let hash_affine = hash.into_affine();

    let signature_one = hash_affine.mul(private_key_one);
    let signature_two = hash_affine.mul(private_key_two);
    let signature_three = hash_affine.mul(private_key_three);

    let mut aggregated_signature = signature_one;
    aggregated_signature.add_assign(&signature_two);
    aggregated_signature.add_assign(&signature_three);

    let pair_one = hash_affine.pairing_with(&aggregated_pub_keys.into_affine());
    let pair_two = generator_one.pairing_with(&aggregated_signature.into_affine());
    
    assert_eq!(pair_one, pair_two);
}

#[test]
fn multi_message_bls_signature(){
    use bellman::bls12_381::{Fr, Bls12, G1, G2, G2Affine};
    use bellman::{PrimeField, PrimeFieldRepr, Engine, CurveAffine, CurveProjective, Field};
    use rescue_poseidon::rescue_hash;
    use rand::{Rand, Rng};
    let mut repr = <Fr as PrimeField>::Repr::default();
    let mut rng = rand::thread_rng();
    let random_number = (&rng.gen::<u32>())%101;
    let mut counter = 0;
    let mut vector_counter = 0;
    let generator_one = <Bls12 as Engine>::G1Affine::one();
    let generator_two = <Bls12 as Engine>::G2Affine::one();
    let mut private_key_vec: Vec<Fr> = Vec::new();
    let mut public_key_vec:  Vec<G1> = Vec::new();
    let mut message_vec:     Vec<G2Affine> = Vec::new();
    let mut signature_vec:   Vec<G2> = Vec::new();

    println!("{:#?} Private Keys will be used", random_number);

    //PRIVATE KEYS CREATED AND STORED
    while counter < random_number {
        let private_key = Fr::rand(&mut rng);
        private_key_vec.push(private_key);
        counter += 1;
    }
    counter = 0;
    
    //PUBLIC KEYS CREATED AND STORED
    while counter < random_number{
        let public_key = generator_one.mul(private_key_vec[vector_counter]);
        public_key_vec.push(public_key);
        vector_counter +=1;
        counter +=1;
    }
    counter = 0;
    vector_counter = 0;

    //RANDOM MESSAGES CREATED, HASHED AND STORED
    while counter < random_number{
        let random_str = (&rng.gen::<u32>()).to_string();
        let bytes = random_str.as_bytes();
        repr.read_le(&bytes[..]);
        let messagefr = Fr::from_repr(repr).unwrap();
        let result = rescue_hash::<Bls12, 1>(&[messagefr]);
        let hash = generator_two.mul(result[1]);
        message_vec.push(hash.into_affine());
        counter += 1;
    }
    counter = 0;
 

    //ALL PUBLIC KEYS CORRESPONDING TO THEIR MESSAGES COMBINED AND CREATED SIGNATURES SEPERETALY
    while counter < random_number{
        let signature = message_vec[vector_counter].mul(private_key_vec[vector_counter]);
        signature_vec.push(signature);
        counter     += 1;
        vector_counter += 1;
    }
    counter = 0;
    vector_counter = 0;

    //AGGREGATING ALL OF THE SIGNATURES
    let mut aggregated_signature = signature_vec[vector_counter];
    while counter < random_number - 1{
        aggregated_signature.add_assign(&signature_vec[vector_counter + 1]);
        counter     += 1;
        vector_counter += 1;
    }
    counter = 0;
    vector_counter = 0;

    //PAIRING AND MULTIPLYING CORRESPONDING PK AND HASHES
    let mut pair_one = message_vec[vector_counter].pairing_with(&public_key_vec[vector_counter].into_affine());
    while counter < random_number - 1{
        pair_one.mul_assign(&message_vec[vector_counter + 1].pairing_with(&public_key_vec[vector_counter + 1].into_affine()));
        counter     += 1;
        vector_counter += 1;
    }

    let pair_two = generator_one.pairing_with(&aggregated_signature.into_affine());

    assert_eq!(pair_one, pair_two);

}


