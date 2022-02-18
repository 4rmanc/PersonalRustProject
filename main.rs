#[test]
fn schnorring(){
   use bellman::bn256::{Fr, Bn256};
    use bellman::{PrimeField, PrimeFieldRepr, Engine, CurveAffine, CurveProjective, Field};
    use rescue_poseidon::rescue_hash;
    use rand::Rand;

    let message = "North Korea's secret nuclear bomb launch code is 69.";
    let mut rng = rand::thread_rng();
    let mut ascii = vec![0; 32];
    let mut buf = vec![0; 32];
    let generator = <Bn256 as Engine>::G1Affine::one();
    let mut repr = <Fr as PrimeField>::Repr::default();
    let mut i = 0;
    let mut counter = 1;
    const L: usize = 1;

    //Storing ASCII values of the message in a way
    for element in message.chars() {
        let a = element as u8;
        if counter < 32 {
            ascii[i] = a;
            i = i+1;
            counter = counter+1;
        } 
        else {
            if i == 32  {i = 0;}
            ascii[i] = (ascii[i] + a)%32;
            i = i+1;
            }
        }


    repr.read_le(&ascii[..]).unwrap();
    let messagefr = Fr::from_repr(repr).unwrap();
    println!("Message to Fr :{:#?}", messagefr);

    //r value created here
    let nonce = Fr::rand(&mut rng);
    println!("Nonce Value :{}", nonce);

    //r*G = R value created here and X coordinate value of R is taken here
    let generated_nonce = generator.mul(nonce);
    println!("Generated Nonce Value :{}", generated_nonce);
    let generated_nonce_affine = generated_nonce.into_affine();
    let (x, _y) = generated_nonce_affine.as_xy();
    let z = x.into_repr();
    z.write_le(&mut buf[..]).unwrap();
    repr.read_le(&buf[..]).unwrap();
    let mut generated_nonce_x = Fr::from_repr(repr).unwrap();
    println!("Generated Nonce's X coordinate :{:#?}", generated_nonce_x);
    
    //Conconating message and R value to a hash function here
    generated_nonce_x.mul_assign(&messagefr);
    let input = [generated_nonce_x; L];
    let result = rescue_hash::<Bn256, L>(&input);
    let mut conconated = result[0];
    conconated.mul_assign(&result[1]);
    println!("Conconated hash value is :{:#?}", conconated);


    //Creating private and public key here
    let private_key = Fr::rand(&mut rng);
    let mut private_key_carrier = private_key;
    println!("Private key :{:#?}", private_key);
    let public_key = generator.mul(private_key);
    let mut public_key_carrier = public_key;      
    println!("Public key :{}", public_key);


    // Schnorr signature p*c+r
    private_key_carrier.mul_assign(&conconated);
    private_key_carrier.add_assign(&nonce);

    //Validating section (Checking s*G =? P*c+R)

    //s*G
    let signature = generator.mul(private_key_carrier);

    //P*c+R
    public_key_carrier.mul_assign(conconated.into_repr());
    public_key_carrier.add_assign(&generated_nonce);

    println!("Signature is :{}", signature);
    println!("Confirmation of signature is :{}", public_key_carrier);
    assert_eq!(public_key_carrier, signature);

}
