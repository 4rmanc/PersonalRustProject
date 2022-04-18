struct Polynomial{
    coefficients: Vec<u64>
}

impl Polynomial{
    //Shows polynomial in the string form
    fn building(&self, x:&str) -> String{
        let mut counter = 0;
        let mut poly= String::new();
        for elements in &self.coefficients{
            let y = format!("{}{}{}", x,"^",counter);
            if counter == 0{
            poly.push_str(&format!("{}{}{}", elements,"*", y));} 
            else{
            if elements >= &0{
            poly.push_str(&format!("{}{}{}{}","+", elements,"*", y));}
            else{
            poly.push_str(&format!("{}{}{}",elements,"*", y));}}
            counter += 1;
        }
        poly
    }
    //Evaluates polynomial in u64 bounds.
    fn eval_poly(&self, x:u64) -> u64{
        let mut poly = Vec::new();
        let mut counter = 0;
        let mut answer:u64 = 0;
        for elements in &self.coefficients{
            let temp = elements * (x.pow(counter));
            poly.push(temp);
            counter += 1;
        }
        for elements in poly{
        answer = answer + elements;
        }
        answer
    }
    //Converts polynomial's coefficients to the corresponding field elements
    fn to_field(&self) -> Vec<bellman::bls12_381::Fr>{
        use bellman::bls12_381::Fr;
        use bellman::{PrimeField, PrimeFieldRepr};
        let mut buffer:Vec<Fr> = Vec::new();
        let mut repr = <Fr as PrimeField>::Repr::default();
        let mut buf: Vec<u8> = vec![0; 32];

        for elements in &self.coefficients{
        buf[0] = *elements as u8;
        repr.read_le(&buf[..]).unwrap();
        let coef_fr = Fr::from_repr(repr).unwrap();
        buffer.push(coef_fr);
        }
        buffer
    }
}
//Evaluates polynomial in the field
fn eval_poly_fr(poly: Vec<bellman::bls12_381::Fr>, x: bellman::bls12_381::Fr) -> bellman::bls12_381::Fr{
    use bellman::bls12_381::{Fr, FrRepr};
    use bellman::Field;
    let temp_x = x;
    let mut counter = 0;
    let mut counter2: u64 = 0;
    let mut answer: Fr = Fr::zero();
    for _i in 0..poly.len(){
        let pow_u64: u64 = counter2;
        let pow_repr = FrRepr::from(pow_u64);
        let power = pow_repr.as_ref();
        let mut temp = poly[counter];
        let temp_y = temp_x.pow(power);
        temp.mul_assign(&temp_y);
        answer.add_assign(&temp);
        counter += 1;
        counter2 += 1;
    }
    answer
}

//FFT operation
fn fft(poly: &Vec<bellman::bls12_381::Fr>, mut rou: bellman::bls12_381::Fr) -> Vec<bellman::bls12_381::Fr>{
    use bellman::bls12_381::{Fr, FrRepr};
    use bellman::Field;
    let n = poly.len()/2;
    let temp_rou: Fr = rou;
    if poly.len() == 1{
        poly.to_vec()
    }
    else{
    let mut counter: u64 = 0;
    let mut evens: Vec<Fr> = Vec::new();
    let mut odds: Vec<Fr> = Vec::new();
    for elements in poly.to_vec(){
        if counter % 2 == 0{
            evens.push(elements);
        }
        else{
            odds.push(elements);
        }
        counter += 1;
    }
    let mut p: Vec<Fr> = vec![Fr::zero(); n*2];
    rou.square();
    let e = fft(&evens, rou);
    let mut o = fft(&odds, rou);
    for i in 0..n {
        let pow_u64: u64 = i as u64;
        let pow_repr = FrRepr::from(pow_u64);
        let power = pow_repr.as_ref();
        let temp_rou2 = temp_rou.pow(power);
        o[i].mul_assign(&temp_rou2);
        p[i].add_assign(&e[i]);
        p[i].add_assign(&o[i]);
        p[i + n].add_assign(&e[i]);
        p[i + n].sub_assign(&o[i]);
    }
    p
    }
}

//Barycentric evaluation
fn barycentric(x:bellman::bls12_381::Fr, eval_poly: &Vec<bellman::bls12_381::Fr>, rou: bellman::bls12_381::Fr) -> bellman::bls12_381::Fr {
    use bellman::bls12_381::{Fr, FrRepr};
    use bellman::{Field, PrimeField, PrimeFieldRepr};
    let n = eval_poly.len();
    let one = Fr::one();
    let mut eval = Fr::zero();
    let mut repr = <Fr as PrimeField>::Repr::default();
    let mut vector:Vec<Fr> = vec![Fr::zero(); n];
    let pow_u64: u64 = n as u64;
    let pow_repr = FrRepr::from(pow_u64);
    let power = pow_repr.as_ref();
    let mut x_carrier = x.pow(power);
    x_carrier.sub_assign(&one);

    for i in 0..n{
        let pow_u64: u64 = i as u64;
        let pow_repr = FrRepr::from(pow_u64);
        let power = pow_repr.as_ref();
        let rou_carrier = rou.pow(power);
        let mut x_carrier = x;
        vector[i] = eval_poly[i];
        vector[i].mul_assign(&rou_carrier);
        x_carrier.sub_assign(&rou_carrier);
        x_carrier = x_carrier.inverse().unwrap();
        vector[i].mul_assign(&x_carrier);
        eval.add_assign(&vector[i]);
    }
    let mut buf: Vec<u8> = vec![0; 32];
    buf[0] = n as u8;
    repr.read_le(&buf[..]).unwrap();
    let mut n_fr = Fr::from_repr(repr).unwrap();
    n_fr = n_fr.inverse().unwrap();
    x_carrier.mul_assign(&n_fr);
    eval.mul_assign(&x_carrier);
    eval
}

//Multiplicator for the field polynomials
fn polynomial_multiplicator(poly_1: &Vec<bellman::bls12_381::Fr>, poly_2: &Vec<bellman::bls12_381::Fr>) ->  Vec<bellman::bls12_381::Fr>{
    use bellman::bls12_381::Fr;
    use bellman::Field;
    let poly_1_len = poly_1.len();
    let poly_2_len = poly_2.len();
    let mut result = vec![Fr::zero(); poly_1_len+poly_2_len-1];
    for i in 0..poly_1_len{
        for j in 0..poly_2_len{
            let mut poly_1_coef = poly_1.to_vec();
            poly_1_coef[i].mul_assign(&poly_2.to_vec()[j]);
            result[i + j].add_assign(&poly_1_coef[i]);
        }
    }
    result
}

//KZG by using barycentric evaluation
fn barycentric_commitment(trusted_secret_g1_vector: &Vec<bellman::bls12_381::G1>, eval_poly: &Vec<bellman::bls12_381::Fr>, rou: bellman::bls12_381::Fr, z: bellman::bls12_381::Fr, y: bellman::bls12_381::Fr) -> (bellman::bls12_381::G1, bellman::bls12_381::G1) {
    use bellman::bls12_381::{Fr, FrRepr, G1};
    use bellman::{Field, PrimeField, PrimeFieldRepr, CurveAffine, CurveProjective};
    let n = eval_poly.len();
    let mut repr = <Fr as PrimeField>::Repr::default();
    let mut vector_px:Vec<Fr> = vec![Fr::zero(); n];
    let mut vector_qx:Vec<Fr> = vec![Fr::zero(); n];
    let mut buf: Vec<u8> = vec![0; 32];
    let mut commit_carrier_px:G1 = G1::zero();
    let mut commit_carrier_qx:G1 = G1::zero();
    buf[0] = n as u8;
    repr.read_le(&buf[..]).unwrap();
    let mut n_fr = Fr::from_repr(repr).unwrap();
    n_fr = n_fr.inverse().unwrap();

    //P(X) and Q(X)
    for i in 0..n{
        let pow_u64: u64 = i as u64;
        let pow_repr = FrRepr::from(pow_u64);
        let power = pow_repr.as_ref();
        let mut rou_carrier = rou.pow(power);
        let mut poly_1 = vec![Fr::one()];
        for j in 0..n{
            let pow_u64: u64 = j as u64;
            let pow_repr = FrRepr::from(pow_u64);
            let power = pow_repr.as_ref();
            let mut rou_carrier_2 = rou.pow(power);
            rou_carrier_2.negate();
            let poly_2 = vec![rou_carrier_2, Fr::one()];
            poly_1 = polynomial_multiplicator(&poly_1, &poly_2);
        }
        rou_carrier.negate();
        let mut poly_px = vec![rou_carrier, Fr::one()];
        rou_carrier.negate();
        poly_px = compute_quotient(poly_1, &poly_px);
        let mut poly_qx = poly_px.clone();
        for a in 0..n{
            let mut temp_px = rou_carrier;
            let mut temp_qx = rou_carrier;
            let mut rou_qx = rou_carrier;
            temp_px.mul_assign(&n_fr);
            temp_qx.mul_assign(&n_fr);
            let carrier_px = eval_poly[i];
            let mut carrier_qx = eval_poly[i];
            carrier_qx.sub_assign(&y);
            rou_qx.sub_assign(&z);
            carrier_qx.mul_assign(&rou_qx.inverse().unwrap());
            temp_qx.mul_assign(&carrier_qx);
            temp_px.mul_assign(&carrier_px);
            poly_qx[a].mul_assign(&temp_qx);
            poly_px[a].mul_assign(&temp_px);
        }
        for i in 0..n{
            vector_px[i].add_assign(&poly_px[i]);
            vector_qx[i].add_assign(&poly_qx[i]);
        }  
    }

    let mut holder_px:G1 = G1::zero();
    let mut holder_qx:G1 = G1::zero();
    for i in 0..n{
        let temp_px = trusted_secret_g1_vector[i].into_affine().mul(vector_px[i]);
        let temp_qx = trusted_secret_g1_vector[i].into_affine().mul(vector_qx[i]);
        holder_px.add_assign(&temp_px); 
        holder_qx.add_assign(&temp_qx); 

    }
    commit_carrier_px.add_assign(&holder_px);
    commit_carrier_qx.add_assign(&holder_qx);

    (commit_carrier_px, commit_carrier_qx)
}

//Long Division for the polynomials
fn compute_quotient(mut dividend: Vec<bellman::bls12_381::Fr>, divisor: &Vec<bellman::bls12_381::Fr>) -> Vec<bellman::bls12_381::Fr> {
    use bellman::Field;
    use bellman::bls12_381::Fr;
    let mut coefficients:Vec<Fr> = Vec::new();
    let mut dividend_pos = dividend.len() - 1;
    let divisor_pos = divisor.len() - 1;
    let mut difference = dividend_pos as isize - divisor_pos as isize;
    while difference >= 0{
        let mut term_quotient = dividend[dividend_pos];
        term_quotient.mul_assign(&divisor[divisor_pos].inverse().unwrap());
        coefficients.push(term_quotient);
        for i in (0..=divisor_pos).rev() {
            let mut x = divisor[i];
            x.mul_assign(&term_quotient);
            dividend[difference as usize + i].sub_assign(&x);
        }
        dividend_pos -= 1;
        difference -= 1;
    }
    coefficients.reverse();
    coefficients

}

//KZG by using coefficients
fn kzg(trusted_secret_g1_vector: &Vec<bellman::bls12_381::G1>, poly_coef: &Vec<bellman::bls12_381::Fr>, mut z: bellman::bls12_381::Fr, y: bellman::bls12_381::Fr) -> (bellman::bls12_381::G1, bellman::bls12_381::G1) {
    use bellman::bls12_381::{Fr, G1};
    use bellman::{Field, CurveAffine, CurveProjective};
    let mut holder: G1 = G1::zero();
    //P(x)
    for i in 0..poly_coef.len(){
    let temp = trusted_secret_g1_vector[i].into_affine().mul(poly_coef[i]);
    holder.add_assign(&temp);
    }
    //Q(x)
    let mut holder_quotient: G1 = G1::zero();
    z.negate();
    let divisor: Vec<Fr> = vec![z, Fr::one()];
    let mut poly_coef_y = poly_coef.to_vec();
    poly_coef_y[0].sub_assign(&y);
    let quotient = compute_quotient(poly_coef_y, &divisor);
    for i in 0..quotient.len(){
        let temp = trusted_secret_g1_vector[i].into_affine().mul(quotient[i]);
        holder_quotient.add_assign(&temp);
    }
    (holder, holder_quotient)
}
#[test]
fn tester(){
    use bellman::bls12_381::{G1, G2, Fr, FrRepr, Bls12};
    use bellman::{PrimeField, Field, Engine, CurveAffine, CurveProjective};
    use rand::Rand;
    let mut rng = rand::thread_rng();
    let generator_one = <Bls12 as Engine>::G1Affine::one();
    let generator_two = <Bls12 as Engine>::G2Affine::one();
    let challenge = Fr::rand(&mut rng);  
    let trusted_setup_secret = Fr::rand(&mut rng);
    let coef:Vec<u64> = vec![2, 6, 231, 221, 26, 98, 10, 11];
    let mut rou = <Fr as PrimeField>::root_of_unity();
    let mut len = coef.len().clone();
    let mut counter = 0;

    while len > 0 {
        len = len/2;
        counter +=1;
    }
    //Calculates correct rou (roots of unity) value to use for this polynomial
    for _i in 0..(33 - counter){
    rou.square();
    }

    let polynom = Polynomial{
        coefficients: coef
    };

    //TRUSTED SETUP
    let mut trusted_setup_g1_vector: Vec<G1> = Vec::new();
    let mut trusted_setup_g2_vector: Vec<G2> = Vec::new();
    for i in 0..polynom.coefficients.len(){
        let pow_u64: u64 = i as u64;
        let pow_repr = FrRepr::from(pow_u64);
        let power = pow_repr.as_ref();
        let trusted_secret_power = trusted_setup_secret.pow(power);
        let trusted_secret_g1 = generator_one.mul(trusted_secret_power.into_repr());
        let trusted_secret_g2 = generator_two.mul(trusted_secret_power.into_repr());
        trusted_setup_g1_vector.push(trusted_secret_g1);
        trusted_setup_g2_vector.push(trusted_secret_g2);
    }

    let fft:Vec<Fr> = fft(&polynom.to_field(), rou);
    let barycentric = barycentric(challenge, &fft, rou);
    let mut barycentric_commitment = barycentric_commitment(&trusted_setup_g1_vector, &fft, rou, challenge, barycentric);
    let mut kzg = kzg(&trusted_setup_g1_vector, &polynom.to_field(), challenge, barycentric);

    //KZG pairing
    let mut trusted_carrier_kzg = trusted_setup_g2_vector[1];
    let challenge_g2 = generator_two.mul(challenge);
    trusted_carrier_kzg.sub_assign(&challenge_g2);
    let pair_one_kzg = kzg.1.into_affine().pairing_with(&trusted_carrier_kzg.into_affine());
    kzg.0.sub_assign(&generator_one.mul(barycentric));
    let pair_two_kzg = kzg.0.into_affine().pairing_with(&generator_two); 

    //Barycentric pairing
    let mut trusted_carrier_bary = trusted_setup_g2_vector[1];
    let challenge_g2 = generator_two.mul(challenge);
    trusted_carrier_bary.sub_assign(&challenge_g2);
    let pair_one_bary = barycentric_commitment.1.into_affine().pairing_with(&trusted_carrier_bary.into_affine());
    barycentric_commitment.0.sub_assign(&generator_one.mul(barycentric));
    let pair_two_bary = barycentric_commitment.0.into_affine().pairing_with(&generator_two); 

    //Test
    assert_eq!("2*x^0+6*x^1+231*x^2+221*x^3+26*x^4+98*x^5+10*x^6+11*x^7", polynom.building("x"));
    assert_eq!(17073191237357, polynom.eval_poly(55));
    assert_eq!(barycentric, eval_poly_fr(polynom.to_field(), challenge));
    assert_eq!(pair_one_kzg, pair_two_kzg);
    assert_eq!(pair_one_bary, pair_two_bary); 
}
