struct Polynomial{
    coefficients: Vec<i64>
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
    //Evaluates polynomial in u8 bounds.
    fn eval_poly(&self, x:i64) -> i64{
        let mut poly = Vec::new();
        let mut counter = 0;
        let mut answer:i64 = 0;
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
    //Converts u8 polynomial's coefficients to the corresponding field elements
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
fn eval_poly_fr(poly: &Vec<bellman::bls12_381::Fr>, x: &bellman::bls12_381::Fr) -> bellman::bls12_381::Fr{
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

fn polynomial_adder(poly_1: &Vec<bellman::bls12_381::Fr>, poly_2: &Vec<bellman::bls12_381::Fr>) ->  Vec<bellman::bls12_381::Fr>{
    use bellman::bls12_381::Fr;
    use bellman::Field;
    let poly_1_len = poly_1.len();
    let poly_2_len = poly_2.len();
    if poly_2_len < poly_1_len{
    let mut result = poly_1.clone();
    for i in 0..poly_2_len{
        result[i].add_assign(&poly_2[i]);
    }
    result
    }
    else{
        let mut result = poly_2.clone();
        for i in 0..poly_1_len{
        result[i].add_assign(&poly_1[i]);
        }
    result
    }
}

fn polynomial_substroctor(poly_1: &Vec<bellman::bls12_381::Fr>, poly_2: &Vec<bellman::bls12_381::Fr>) ->  Vec<bellman::bls12_381::Fr>{
    use bellman::bls12_381::Fr;
    use bellman::Field;
    let poly_1_len = poly_1.len();
    let poly_2_len = poly_2.len();
    
    let mut result = poly_1.clone();
    for i in 0..poly_2_len{
        result[i].sub_assign(&poly_2[i]);
    }
    result
    
}


//Computer for the quotient polynomial
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

fn lagrange_polynomial(eval_poly: &Vec<bellman::bls12_381::Fr>, rou: bellman::bls12_381::Fr) -> Vec<bellman::bls12_381::Fr>{
    use bellman::bls12_381::{Fr, FrRepr};
    use bellman::{Field, PrimeField, PrimeFieldRepr};
    let n = eval_poly.len();
    let mut repr = <Fr as PrimeField>::Repr::default();
    let mut vector_px:Vec<Fr> = vec![Fr::zero(); n];

    //Calculating our N by calculating the length of the polynomial and taking inverse immediately to use it in the future
    let mut buf: Vec<u8> = vec![0; 32];
    buf[0] = n as u8;
    repr.read_le(&buf[..]).unwrap();
    let mut n_fr = Fr::from_repr(repr).unwrap();
    n_fr = n_fr.inverse().unwrap();

    //Calculating our Lagrange Polynomial for given eval points
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
            //making our xi negative
            rou_carrier_2.negate();
            //poly_2 is basically (x-xi) and xi is our roots of unity
            let poly_2 = vec![rou_carrier_2, Fr::one()];
            //Here calculating M(x) which is poly_1
            poly_1 = polynomial_multiplicator(&poly_1, &poly_2);
        }

        
        rou_carrier.negate();
        //poly_px is x-wi
        let mut poly_px = vec![rou_carrier, Fr::one()];
        rou_carrier.negate();
        //calculating M(X)/x-wi
        poly_px = compute_quotient(poly_1, &poly_px);
        for a in 0..n{
            let mut temp_px = rou_carrier;
            temp_px.mul_assign(&n_fr);
            let carrier_px = eval_poly[i];
            temp_px.mul_assign(&carrier_px);
            poly_px[a].mul_assign(&temp_px);
        }
        for i in 0..n{
            vector_px[i].add_assign(&poly_px[i]);
        }  
    }
    vector_px
}

fn commit_maker(poly: &Vec<bellman::bls12_381::Fr>, trusted_setup_g1_vector: &Vec<bellman::bls12_381::G1>, trusted_setup_g2_vector: &Vec<bellman::bls12_381::G2>) -> (bellman::bls12_381::G1, bellman::bls12_381::G2){
    use bellman::bls12_381::{Fr, G1, G2};
    use bellman::{CurveAffine, CurveProjective};
    let mut commit_result_g1: G1 = G1::zero();
    let mut commit_result_g2: G2 = G2::zero();
    let temp: Vec<Fr> = poly.clone();


    for i in 0..poly.len(){
        let temp_g1 = trusted_setup_g1_vector[i].into_affine().mul(temp[i]);
        commit_result_g1.add_assign(&temp_g1);
        let temp_g2 = trusted_setup_g2_vector[i].into_affine().mul(temp[i]);
        commit_result_g2.add_assign(&temp_g2);
    }

    (commit_result_g1,
     commit_result_g2)

}


fn main (){
    use bellman::plonk::better_cs::generator::make_non_residues;
    use bellman::bls12_381::{G1, G2, Fr, FrRepr, Bls12};
    use bellman::{PrimeField, Field, Engine, CurveAffine, CurveProjective};
    use rand::Rand;
    let mut rng = rand::thread_rng();

    //G1 and G2
    let generator_one = <Bls12 as Engine>::G1Affine::one();
    let generator_two = <Bls12 as Engine>::G2Affine::one();

    //Our Polynomial
    let coef:Vec<i64> = vec![2, 76, 11, 47];
    let polynom = Polynomial{
        coefficients: vec![2, 76, 11, 47]
    };

    //Calculates correct rou (roots of unity) value to use for this polynomial
    let mut len = coef.len();
    let mut counter = 0;
    while len > 0 {
        len = len/2;
        counter +=1;
    }

    let mut rou = <Fr as PrimeField>::root_of_unity();
    for _i in 0..(33 - counter){
    rou.square();
    }

    //Creating k1 and k2. We will use them to shift our ROU
    let x: usize = 2;
    let k_values: Vec<Fr> = make_non_residues(x);

    //TRUSTED SETUP
    let trusted_setup_secret = Fr::rand(&mut rng);
    let mut trusted_setup_g1_vector: Vec<G1> = Vec::new();
    let mut trusted_setup_g2_vector: Vec<G2> = Vec::new();
    for i in 0..polynom.coefficients.len() + 5{
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

    let poly_fr  = polynom.to_field();
    let poly_lag = lagrange_polynomial(&fft, rou);

    //Public Circuit Values (added two more zeroes to end of the each q because I want to hit 16 which is 2^4)
    let mut negone = Fr::one();
    negone.negate();

    let q_left =     Polynomial{ coefficients: vec![0, 0, 0, 1] };
    let q_right =    Polynomial{ coefficients: vec![0, 0, 0, 1] };
    let q_output: Vec<Fr> = vec![negone, negone, negone, negone];
    let q_constant = Polynomial{ coefficients: vec![0, 0, 0, 0] };
    let q_multi =    Polynomial{ coefficients: vec![1, 1, 1, 0] };

    //Public Circuit Values to Lagrange
    let q_left_lag: Vec<Fr> =     lagrange_polynomial(&q_left.to_field(), rou);
    let q_right_lag: Vec<Fr> =    lagrange_polynomial(&q_right.to_field(), rou);
    let q_output_lag: Vec<Fr> =   lagrange_polynomial(&q_output, rou);
    let q_constant_lag: Vec<Fr> = lagrange_polynomial(&q_constant.to_field(), rou);
    let q_multi_lag: Vec<Fr> =    lagrange_polynomial(&q_multi.to_field(), rou);

    
    //Witness values
    let a_coef = vec![3, 4, 5, 9];
    let b_coef = vec![3, 4, 5, 16];
    let c_coef = vec![9, 16, 25, 25];
 
    let a = Polynomial{ coefficients: a_coef.clone() };
    let b = Polynomial{ coefficients: b_coef.clone() };
    let c = Polynomial{ coefficients: c_coef.clone() };

    //Our Z(x) Polynomial.
    let z_h: Vec<Fr> = vec![negone, Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()];



    //Blinding values
    let blinding_scalars: Vec<Fr> = vec![Fr::rand(&mut rng), Fr::rand(&mut rng), Fr::rand(&mut rng), Fr::rand(&mut rng),
    Fr::rand(&mut rng), Fr::rand(&mut rng), Fr::rand(&mut rng), Fr::rand(&mut rng), Fr::rand(&mut rng)];

    //Witness Values to Lagrange
    let temp_a_blind: Vec<Fr> = polynomial_multiplicator(&vec![blinding_scalars[0], blinding_scalars[1]], &z_h);
    let temp_b_blind: Vec<Fr> = polynomial_multiplicator(&vec![blinding_scalars[2], blinding_scalars[3]], &z_h);
    let temp_c_blind: Vec<Fr> = polynomial_multiplicator(&vec![blinding_scalars[4], blinding_scalars[5]], &z_h);
    let a_carrier: Vec<Fr> = lagrange_polynomial(&a.to_field(), rou);
    let b_carrier: Vec<Fr> = lagrange_polynomial(&b.to_field(), rou);
    let c_carrier: Vec<Fr> = lagrange_polynomial(&c.to_field(), rou);

    let a_lag: Vec<Fr> = polynomial_adder(&temp_a_blind, &a_carrier);
    let b_lag: Vec<Fr> = polynomial_adder(&temp_b_blind, &b_carrier);
    let c_lag: Vec<Fr> = polynomial_adder(&temp_c_blind, &c_carrier);


   //Calculating Copy Constraints
   let mut rou_vec: Vec<Fr> = vec![Fr::zero(); q_output_lag.len()];
   let mut k1_vec:  Vec<Fr> = vec![Fr::zero(); q_output_lag.len()];
   let mut k2_vec:  Vec<Fr> = vec![Fr::zero(); q_output_lag.len()];

   for i in 0..q_output_lag.len(){
       let pow_u64: u64 = i as u64;
       let pow_repr     = FrRepr::from(pow_u64);
       let power        = pow_repr.as_ref();
       rou_vec[i] = rou.pow(power);
       k1_vec[i]  = rou.pow(power);
       k2_vec[i]  = rou.pow(power);
       k1_vec[i].mul_assign(&k_values[0]);
       k2_vec[i].mul_assign(&k_values[1]);
   }
   let mut rou_vec_temp: Vec<Fr> = rou_vec.clone();
   let mut k1_vec_temp: Vec<Fr>  = k1_vec.clone();
   let mut k2_vec_temp: Vec<Fr>  = k2_vec.clone();

   //We will use this one in round 3
   let rou_vec_temp2: Vec<Fr> = rou_vec.clone();
   let k1_vec_temp2: Vec<Fr>  = k1_vec.clone();
   let k2_vec_temp2: Vec<Fr>  = k2_vec.clone();

   //I wont do a hard coding here. I will do a little loop to construct our sigma by using witness values. In a normal scenario this should be precalculated without using witness values
   for i in 0..q_left_lag.len(){
        let mut temp: Fr = Fr::zero();
        let flag = 0;
       for j in 0..q_left_lag.len(){
           //These & operations checks the algorithm to change every same value once so it wont be same value again
            if (a_coef[i] == a_coef[j]) & (rou_vec[i] != rou_vec_temp[j]){
                temp = rou_vec[i];
                rou_vec[i] = rou_vec[j];
                rou_vec[j] = temp;
            }
            if (b_coef[i] == b_coef[j]) & (k1_vec[i] != k1_vec_temp[j]){
                temp = k1_vec[i];
                k1_vec[i] = k1_vec[j];
                k1_vec[j] = temp;
            }
            if (c_coef[i] == c_coef[j]) & (k2_vec[i] != k2_vec_temp[j]){
                temp = k2_vec[i];
                k2_vec[i] = k2_vec[j];
                k2_vec[j] = temp;
            }
            
            if (a_coef[i] == c_coef[j]) & (rou_vec[i] != k2_vec_temp[j]){
                temp = rou_vec[i];
                rou_vec[i] = k2_vec[j];
                k2_vec[j] = temp;
            }
            if (a_coef[i] == b_coef[j]) & (rou_vec[i] != k1_vec_temp[j]){
                temp = rou_vec[i];
                rou_vec[i] = k1_vec[j];
                k1_vec[j] = temp;
            }
            if (c_coef[i] == b_coef[j]) & (k2_vec[i] != k1_vec_temp[j]){
                temp = k2_vec[i];
                k2_vec[i] = k1_vec[j];
                k1_vec[j] = temp;  
            }
            
       }
   }
   let rou_copy_temp: Vec<Fr> = rou_vec.clone();
   let k1_copy_temp: Vec<Fr>  = k1_vec.clone();
   let k2_copy_temp: Vec<Fr>  = k2_vec.clone();

   //Calculating our sigma polynomials. Evaluation points are the vectors we calculated above and domains are the ROUs
   let sigma1 = lagrange_polynomial(&rou_vec, rou);
   let sigma2 = lagrange_polynomial(&k1_vec,  rou);
   let sigma3 = lagrange_polynomial(&k2_vec,  rou);

   //ROUND 1 Calculating commitments for the a, b and c polynomials. We will have both g1 and g2 commitments for the values.
   let commit_a = commit_maker(&a_lag, &trusted_setup_g1_vector, &trusted_setup_g2_vector);
   let commit_b = commit_maker(&b_lag, &trusted_setup_g1_vector, &trusted_setup_g2_vector);
   let commit_c = commit_maker(&c_lag, &trusted_setup_g1_vector, &trusted_setup_g2_vector);

   //ROUND 2 Calculating the grand product Z(X)
   //Finding L1(x) which is 1 for first element in rou and 0 for others
    let mut lag_coef: Vec<Fr> = vec![Fr::zero(); q_left_lag.len()];
    lag_coef[0] = Fr::one();
    let lag1 = lagrange_polynomial(&lag_coef, rou);
    lag_coef[0] = Fr::zero();

    let gamma = Fr::rand(&mut rng); 
    let beta  = Fr::rand(&mut rng);  

    //Calculating accumulator
    let mut carrier: Vec<Fr> = vec![Fr::zero(); q_left_lag.len()];
    for i in 0..q_left_lag.len(){
    let mut temp_a = Fr::zero();
    let mut temp_b = Fr::zero();
    let mut temp_c = Fr::zero();
    let mut accumulator_carrier = Fr::one();
        
    lag_coef[i] = Fr::one();
    let mut lagranges: Vec<Fr> = lagrange_polynomial(&lag_coef, rou);
    lag_coef[i] = Fr::zero();

    for j in 0..i{
        rou_vec_temp = rou_vec_temp2.clone();
        k1_vec_temp = k1_vec_temp2.clone();
        k2_vec_temp = k2_vec_temp2.clone();
        //(wj+beta*rou+gamma)
        //Numerator
        temp_a = eval_poly_fr(&a_lag, &rou_vec_temp2[j]);
        rou_vec_temp[j].mul_assign(&beta);
        temp_a.add_assign(&rou_vec_temp[j]);
        temp_a.add_assign(&gamma);

        temp_b = eval_poly_fr(&b_lag, &rou_vec_temp2[j]);
        k1_vec_temp[j].mul_assign(&beta);
        temp_b.add_assign(&k1_vec_temp[j]);
        temp_b.add_assign(&gamma);
    
        temp_c = eval_poly_fr(&c_lag, &rou_vec_temp2[j]);
        k2_vec_temp[j].mul_assign(&beta);
        temp_c.add_assign(&k2_vec_temp[j]);
        temp_c.add_assign(&gamma);

        temp_a.mul_assign(&temp_b);
        temp_a.mul_assign(&temp_c);

        accumulator_carrier.mul_assign(&temp_a);

        //Denominator
        rou_vec = rou_copy_temp.clone();
        k1_vec = k1_copy_temp.clone();
        k2_vec = k2_copy_temp.clone();

        temp_a = eval_poly_fr(&a_lag, &rou_vec_temp2[j]);
        rou_vec[j].mul_assign(&beta);
        temp_a.add_assign(&rou_vec[j]);
        temp_a.add_assign(&gamma);

        temp_b = eval_poly_fr(&b_lag, &rou_vec_temp2[j]);
        k1_vec[j].mul_assign(&beta);
        temp_b.add_assign(&k1_vec[j]);
        temp_b.add_assign(&gamma);

        temp_c = eval_poly_fr(&c_lag, &rou_vec_temp2[j]);
        k2_vec[j].mul_assign(&beta);
        temp_c.add_assign(&k2_vec[j]);
        temp_c.add_assign(&gamma);

        temp_a.mul_assign(&temp_b);
        temp_a.mul_assign(&temp_c);
        temp_a.inverse().unwrap();

        accumulator_carrier.mul_assign(&temp_a);

        }
        for i in 0..lagranges.len(){
            lagranges[i].mul_assign(&accumulator_carrier);
        }
        //Here adding every new polynomial to carrier to keep them as a grand product
        for i in 0..q_left_lag.len(){
            carrier[i].add_assign(&lagranges[i]);
        }
    }
    let mut grand_blinder = polynomial_multiplicator(&vec![blinding_scalars[6], blinding_scalars[7], blinding_scalars[8]], &z_h);
    for i in 0..carrier.len(){
        grand_blinder[i].add_assign(&carrier[i]);
    }

    //Grand Product Z(X)
    let grand_product_poly = grand_blinder;
    let grand_product = commit_maker(&grand_product_poly, &trusted_setup_g1_vector, &trusted_setup_g2_vector);

    //ROUND 3
    let alfa = Fr::rand(&mut rng);
    let teta = Fr::rand(&mut rng);
    let mut quotient_polynomial: Vec<Fr> = Vec::new();

    let mut a_b_qm = polynomial_multiplicator(&a_lag, &b_lag);
    a_b_qm = polynomial_multiplicator(&a_b_qm, &q_multi_lag);
    let a_ql = polynomial_multiplicator(&a_lag, &q_left_lag);
    let b_qr = polynomial_multiplicator(&b_lag, &q_right_lag);
    let c_qo = polynomial_multiplicator(&c_lag, &q_output_lag);

    let mut quot_first_row = polynomial_adder(&a_b_qm, &a_ql);
    quot_first_row = polynomial_adder(&quot_first_row, &b_qr);
    quot_first_row = polynomial_adder(&quot_first_row, &c_qo);
    quot_first_row = polynomial_adder(&quot_first_row, &q_constant_lag); 

    let mut temp_a_row2 = a_lag.clone();
    temp_a_row2[1].add_assign(&beta);
    temp_a_row2[0].add_assign(&gamma);

    let mut temp_b_row2 = b_lag.clone();
    let mut temp_beta_k1 = beta.clone();
    temp_beta_k1.mul_assign(&k_values[0]);
    temp_b_row2[1].add_assign(&temp_beta_k1);
    temp_b_row2[0].add_assign(&gamma);

    let mut temp_c_row2 = c_lag.clone();
    let mut temp_beta_k2 = beta.clone();
    temp_beta_k2.mul_assign(&k_values[1]);
    temp_c_row2[1].add_assign(&temp_beta_k2);
    temp_c_row2[0].add_assign(&gamma);

    let mut quot_second_row = polynomial_multiplicator(&grand_product_poly, &temp_a_row2);
    quot_second_row = polynomial_multiplicator(&quot_second_row, &temp_b_row2);
    quot_second_row = polynomial_multiplicator(&quot_second_row, &temp_c_row2);
    quot_second_row = polynomial_multiplicator(&quot_second_row, &vec![alfa]);

    let mut temp_a_row3 = polynomial_adder(&a_lag, &polynomial_multiplicator(&sigma1, &vec![beta]));
    temp_a_row3[0].add_assign(&gamma);

    let mut temp_b_row3 = polynomial_adder(&b_lag, &polynomial_multiplicator(&sigma2, &vec![beta]));
    temp_b_row3[0].add_assign(&gamma);

    let mut temp_c_row3 = polynomial_adder(&c_lag, &polynomial_multiplicator(&sigma3, &vec![beta]));
    temp_c_row3[0].add_assign(&gamma);

    let mut grand_product_poly_w = grand_product_poly.clone();
    for i in 0..grand_product_poly_w.len(){
        if i < rou_vec_temp2.len(){
        grand_product_poly_w[i].mul_assign(&rou_vec_temp2[i]);
        }
        else{
        grand_product_poly_w[i].mul_assign(&rou_vec_temp2[i - rou_vec_temp2.len()]);
        }
    }

    let mut quot_third_row = polynomial_multiplicator(&temp_a_row3, &temp_b_row3);
    quot_third_row = polynomial_multiplicator(&quot_third_row, &temp_c_row3);
    quot_third_row = polynomial_multiplicator(&quot_third_row, &grand_product_poly_w);
    quot_third_row = polynomial_multiplicator(&quot_third_row, &vec![alfa]);

    let mut quot_fourth_row = grand_product_poly.clone();
    quot_fourth_row = polynomial_multiplicator(&quot_fourth_row, &lag1);
    quot_fourth_row = polynomial_substroctor(&quot_fourth_row, &lag1);
    quot_fourth_row = polynomial_multiplicator(&quot_fourth_row, &vec![alfa]);
    quot_fourth_row = polynomial_multiplicator(&quot_fourth_row, &vec![alfa]);

    quotient_polynomial = polynomial_adder(&quot_first_row, &quot_second_row);
    quotient_polynomial = polynomial_substroctor(&quotient_polynomial, &quot_third_row);
    quotient_polynomial = polynomial_adder(&quotient_polynomial, &quot_fourth_row);
    quotient_polynomial = compute_quotient(quotient_polynomial, &z_h);

    let mut quotient_low  = vec![Fr::zero(); quotient_polynomial.len() / 3];
    let mut quotient_mid  = vec![Fr::zero(); quotient_polynomial.len() / 3];
    let mut quotient_high = vec![Fr::zero(); quotient_polynomial.len() / 3];
    for i in 0..(quotient_polynomial.len() / 3){
        quotient_low[i]  = quotient_polynomial[i];
        quotient_mid[i]  = quotient_polynomial[i + (quotient_polynomial.len() / 3)];
        quotient_high[i] = quotient_polynomial[i + ((quotient_polynomial.len() / 3) * 2)];
    }
    let quotient_low_commit  = commit_maker(&quotient_low, &trusted_setup_g1_vector, &trusted_setup_g2_vector);
    let quotient_mid_commit  = commit_maker(&quotient_mid, &trusted_setup_g1_vector, &trusted_setup_g2_vector);
    let quotient_high_commit = commit_maker(&quotient_high, &trusted_setup_g1_vector, &trusted_setup_g2_vector);

    //ROUND 4
    let a_open = eval_poly_fr(&a_lag, &teta);
    let b_open = eval_poly_fr(&b_lag, &teta);
    let c_open = eval_poly_fr(&c_lag, &teta);
    let sigma1_open = eval_poly_fr(&sigma1, &teta);
    let sigma2_open = eval_poly_fr(&sigma2, &teta);
    let quotient_open = eval_poly_fr(&quotient_polynomial, &teta);
    let mut rou_teta = teta.clone();
    rou_teta.mul_assign(&rou);
    let grand_product_open = eval_poly_fr(&grand_product_poly, &rou_teta);

    //ROUND 5
    let delta = Fr::rand(&mut rng);

    let mut ab_qm_linear = polynomial_multiplicator(&q_multi_lag, &vec![a_open]);
    ab_qm_linear = polynomial_multiplicator(&ab_qm_linear, &vec![b_open]);

    let a_ql_linear = polynomial_multiplicator(&q_left_lag, &vec![a_open]);

    let b_qr_linear = polynomial_multiplicator(&q_right_lag, &vec![b_open]);

    let c_qo_linear = polynomial_multiplicator(&q_output_lag, &vec![c_open]);

    let mut linear_first_row = polynomial_adder(&ab_qm_linear, &a_ql_linear);
    linear_first_row = polynomial_adder(&linear_first_row, &b_qr_linear);
    linear_first_row = polynomial_adder(&linear_first_row, &c_qo_linear);
    linear_first_row = polynomial_adder(&linear_first_row, &q_constant_lag);

    let mut temp_a_teta = teta.clone();
    temp_a_teta.mul_assign(&beta);
    temp_a_teta.add_assign(&a_open);
    temp_a_teta.add_assign(&gamma);

    let mut temp_b_teta = teta.clone();
    temp_b_teta.mul_assign(&beta);
    temp_b_teta.mul_assign(&k_values[0]);
    temp_b_teta.add_assign(&b_open);
    temp_b_teta.add_assign(&gamma);

    let mut temp_c_teta = teta.clone();
    temp_c_teta.mul_assign(&beta);
    temp_c_teta.mul_assign(&k_values[1]);
    temp_c_teta.add_assign(&c_open);
    temp_c_teta.add_assign(&gamma);

    //DIDNT MULTIPLIED WITH ALFA YET!!
    let mut linear_second_row = polynomial_multiplicator(&grand_product_poly, &vec![temp_a_teta]);
    linear_second_row = polynomial_multiplicator(&linear_second_row, &vec![temp_b_teta]);
    linear_second_row = polynomial_multiplicator(&linear_second_row, &vec![temp_c_teta]);

    let mut temp_a_sigma = sigma1_open.clone();
    temp_a_sigma.mul_assign(&beta);
    temp_a_sigma.add_assign(&a_open);
    temp_a_sigma.add_assign(&gamma);

    let mut temp_b_sigma = sigma2_open.clone();
    temp_b_sigma.mul_assign(&beta);
    temp_b_sigma.add_assign(&b_open);
    temp_b_sigma.add_assign(&gamma);

    let temp_c_sigma = polynomial_multiplicator(&sigma3, &vec![beta]);

    let mut linear_third_row = polynomial_multiplicator(&temp_c_sigma, &vec![temp_a_sigma]);
    linear_third_row = polynomial_multiplicator(&linear_third_row, &vec![temp_b_sigma]);
    linear_third_row = polynomial_multiplicator(&linear_third_row, &vec![grand_product_open]);

    let linear_second_third_merge = polynomial_multiplicator(&polynomial_substroctor(&linear_second_row, &linear_third_row), &vec![alfa]);

    let mut linear_grand_temp = grand_product_poly.clone();
    linear_grand_temp = polynomial_multiplicator(&linear_grand_temp, &vec![eval_poly_fr(&lag1, &teta)]);
    linear_grand_temp = polynomial_multiplicator(&linear_grand_temp, &vec![alfa]);
    let linear_fourth_row = polynomial_multiplicator(&linear_grand_temp, &vec![alfa]);

    let pow_u64: u64 = (quotient_high.len()*2) as u64;
    let pow_repr     = FrRepr::from(pow_u64);
    let power        = pow_repr.as_ref();
    let teta_over_2n = teta.pow(power);

    let pow_u64: u64 = quotient_mid.len() as u64;
    let pow_repr     = FrRepr::from(pow_u64);
    let power        = pow_repr.as_ref();
    let teta_over_n = teta.pow(power);
 
    let mut linearisation_polynomial = polynomial_adder(&linear_first_row, &linear_second_third_merge);
    linearisation_polynomial = polynomial_adder(&linearisation_polynomial, &linear_fourth_row);
  
    let r_open = eval_poly_fr(&linearisation_polynomial, &teta);

    //Compute opening proof polynomial
    let mut temp_opening_proof_linearisation = linearisation_polynomial.clone();
    let mut temp_opening_proof_a = a_lag.clone();
    let mut temp_opening_proof_b = b_lag.clone();
    let mut temp_opening_proof_c = c_lag.clone();
    let mut temp_opening_proof_sigma1 = sigma1.clone();
    let mut temp_opening_proof_sigma2 = sigma2.clone();

    temp_opening_proof_linearisation[0].sub_assign(&r_open);
    temp_opening_proof_a[0].sub_assign(&a_open);
    temp_opening_proof_b[0].sub_assign(&b_open);
    temp_opening_proof_c[0].sub_assign(&c_open);
    temp_opening_proof_sigma1[0].sub_assign(&sigma1_open);
    temp_opening_proof_sigma2[0].sub_assign(&sigma2_open);
    
    //delta^n
    let mut delta_vec = Vec::new();
    for i in 1..7{
        let pow_u64: u64 = i as u64;
        let pow_repr     = FrRepr::from(pow_u64);
        let power        = pow_repr.as_ref();
        delta_vec.push(delta.pow(power));
    }

    temp_opening_proof_linearisation = polynomial_multiplicator(&temp_opening_proof_linearisation, &vec![delta_vec[0]]);
    temp_opening_proof_a = polynomial_multiplicator(&temp_opening_proof_a, &vec![delta_vec[1]]);
    temp_opening_proof_b = polynomial_multiplicator(&temp_opening_proof_b, &vec![delta_vec[2]]);
    temp_opening_proof_c = polynomial_multiplicator(&temp_opening_proof_c, &vec![delta_vec[3]]);
    temp_opening_proof_sigma1 = polynomial_multiplicator(&temp_opening_proof_sigma1, &vec![delta_vec[4]]);
    temp_opening_proof_sigma2 = polynomial_multiplicator(&temp_opening_proof_sigma2, &vec![delta_vec[5]]);

    let opening_quot_high = polynomial_multiplicator(&quotient_high, &vec![teta_over_2n]);
    let opening_quot_mid  = polynomial_multiplicator(&quotient_mid, &vec![teta_over_n]);
    let mut temp_opening_quot = polynomial_adder(&quotient_low, &opening_quot_mid);
    temp_opening_quot = polynomial_adder(&temp_opening_quot, &opening_quot_high);
    temp_opening_quot[0].sub_assign(&quotient_open);

    let mut opening_proof_polynomial = polynomial_adder(&temp_opening_quot, &temp_opening_proof_linearisation);
    opening_proof_polynomial = polynomial_adder(&opening_proof_polynomial, &temp_opening_proof_a);
    opening_proof_polynomial = polynomial_adder(&opening_proof_polynomial, &temp_opening_proof_b);
    opening_proof_polynomial = polynomial_adder(&opening_proof_polynomial, &temp_opening_proof_c);
    opening_proof_polynomial = polynomial_adder(&opening_proof_polynomial, &temp_opening_proof_sigma1);
    opening_proof_polynomial = polynomial_adder(&opening_proof_polynomial, &temp_opening_proof_sigma2);

    let mut negated_teta = teta.clone();
    negated_teta.negate();

    opening_proof_polynomial = compute_quotient(opening_proof_polynomial, &vec![negated_teta, Fr::one()]);

    let mut opening_proof_polynomial_rou = grand_product_poly.clone();
    opening_proof_polynomial_rou[0].sub_assign(&grand_product_open);

    let mut negated_rou_teta = teta.clone();
    negated_rou_teta.mul_assign(&rou);
    negated_rou_teta.negate();

    opening_proof_polynomial_rou = compute_quotient(opening_proof_polynomial_rou, &vec![negated_rou_teta, Fr::one()]);


    let opening_proof_polynomial_commit = commit_maker(&opening_proof_polynomial, &trusted_setup_g1_vector, &trusted_setup_g2_vector);
    let opening_proof_polynomial_rou_commit = commit_maker(&opening_proof_polynomial_rou, &trusted_setup_g1_vector, &trusted_setup_g2_vector);


    //VERIFIER SECTION
    //Calculating R0
    let upsilon = Fr::rand(&mut rng); 
    let mut r0 = r_open.clone();

    let mut r0_lag = alfa.clone();
    r0_lag.mul_assign(&alfa);
    r0_lag.mul_assign(&eval_poly_fr(&lag1, &teta));
    
    let mut temp_r0_a = beta.clone();
    temp_r0_a.mul_assign(&sigma1_open);
    temp_r0_a.add_assign(&a_open);
    temp_r0_a.add_assign(&gamma);

    let mut temp_r0_b = beta.clone();
    temp_r0_b.mul_assign(&sigma2_open);
    temp_r0_b.add_assign(&b_open);
    temp_r0_b.add_assign(&gamma);

    let mut temp_r0_c = gamma.clone();
    temp_r0_c.add_assign(&c_open);

    temp_r0_a.mul_assign(&temp_r0_b);
    temp_r0_a.mul_assign(&temp_r0_c);
    temp_r0_a.mul_assign(&grand_product_open);
    temp_r0_a.mul_assign(&alfa);

    r0.sub_assign(&r0_lag);
    r0.sub_assign(&temp_r0_a);
    let zh_divid = eval_poly_fr(&z_h, &teta);

    let mut q_multi_commit    = commit_maker(&q_multi_lag, &trusted_setup_g1_vector, &trusted_setup_g2_vector);
    let mut q_left_commit     = commit_maker(&q_left_lag, &trusted_setup_g1_vector, &trusted_setup_g2_vector);
    let mut q_right_commit    = commit_maker(&q_right_lag, &trusted_setup_g1_vector, &trusted_setup_g2_vector);
    let mut q_output_commit   = commit_maker(&q_output_lag, &trusted_setup_g1_vector, &trusted_setup_g2_vector);
    let mut q_constant_commit     = commit_maker(&q_constant_lag, &trusted_setup_g1_vector, &trusted_setup_g2_vector);

    q_multi_commit.0.mul_assign(a_open.into_repr());
    q_multi_commit.0.mul_assign(b_open.into_repr());
    q_multi_commit.0.mul_assign(delta.into_repr());
    
    q_left_commit.0.mul_assign(a_open.into_repr());
    q_left_commit.0.mul_assign(delta.into_repr());
    
    q_right_commit.0.mul_assign(b_open.into_repr());
    q_right_commit.0.mul_assign(delta.into_repr());

    q_output_commit.0.mul_assign(c_open.into_repr());
    q_output_commit.0.mul_assign(delta.into_repr());

    q_constant_commit.0.mul_assign(delta.into_repr());

    let mut batched_first_row = q_multi_commit.0.clone();
    batched_first_row.add_assign(&q_left_commit.0);
    batched_first_row.add_assign(&q_right_commit.0);
    batched_first_row.add_assign(&q_output_commit.0);
    batched_first_row.add_assign(&q_constant_commit.0);
    
    let mut batched_a_temp = teta.clone();
    batched_a_temp.mul_assign(&beta);
    batched_a_temp.add_assign(&a_open);
    batched_a_temp.add_assign(&gamma);

    let mut batched_b_temp = teta.clone();
    batched_b_temp.mul_assign(&beta);
    batched_b_temp.mul_assign(&k_values[0]);
    batched_b_temp.add_assign(&b_open);
    batched_b_temp.add_assign(&gamma);

    let mut batched_c_temp = teta.clone();
    batched_c_temp.mul_assign(&beta);
    batched_c_temp.mul_assign(&k_values[1]);
    batched_c_temp.add_assign(&c_open);
    batched_c_temp.add_assign(&gamma);

    let mut batched_lag1 = eval_poly_fr(&lag1, &teta);
    batched_lag1.mul_assign(&alfa);
    batched_lag1.mul_assign(&alfa);
    batched_lag1.mul_assign(&delta);
    batched_lag1.add_assign(&upsilon);

    batched_a_temp.mul_assign(&batched_b_temp);
    batched_a_temp.mul_assign(&batched_c_temp);
    batched_a_temp.mul_assign(&alfa);
    batched_a_temp.mul_assign(&delta);

    batched_lag1.add_assign(&batched_a_temp);

    let mut batched_second_row = grand_product.0.clone();
    batched_second_row.mul_assign(batched_lag1);

    let mut batched_sigma1 = sigma1_open.clone();
    batched_sigma1.mul_assign(&beta);
    batched_sigma1.add_assign(&a_open);
    batched_sigma1.add_assign(&gamma);

    let mut batched_sigma2 = sigma2_open.clone();
    batched_sigma2.mul_assign(&beta);
    batched_sigma2.add_assign(&b_open);
    batched_sigma2.add_assign(&gamma);

    batched_sigma1.mul_assign(&batched_sigma2);
    batched_sigma1.mul_assign(&alfa);
    batched_sigma1.mul_assign(&delta);
    batched_sigma1.mul_assign(&beta);
    batched_sigma1.mul_assign(&grand_product_open);

    let sigma3_commit = commit_maker(&sigma3, &trusted_setup_g1_vector, &trusted_setup_g2_vector);
    let mut batched_third_row = sigma3_commit.0.clone();
    batched_third_row.mul_assign(batched_sigma1);

    let mut batched_polynomial_commitment = batched_first_row;
    batched_polynomial_commitment.add_assign(&batched_second_row);
    batched_polynomial_commitment.sub_assign(&batched_third_row);

    let sigma1_commit = commit_maker(&sigma1, &trusted_setup_g1_vector, &trusted_setup_g2_vector);
    let sigma2_commit = commit_maker(&sigma2, &trusted_setup_g1_vector, &trusted_setup_g2_vector);

    let mut full_batch_quotient_low = quotient_low_commit.0.clone();
    let mut full_batch_quotient_mid = quotient_mid_commit.0.clone();
    let mut full_batch_quotient_high = quotient_high_commit.0.clone();
    full_batch_quotient_mid.mul_assign(teta_over_n);
    full_batch_quotient_high.mul_assign(teta_over_2n);
    full_batch_quotient_low.add_assign(&full_batch_quotient_mid);
    full_batch_quotient_low.add_assign(&full_batch_quotient_high);

    let mut full_batch_a = commit_a.0.clone();
    full_batch_a.mul_assign(delta_vec[1]);

    let mut full_batch_b = commit_b.0.clone();
    full_batch_b.mul_assign(delta_vec[2]);

    let mut full_batch_c = commit_c.0.clone();
    full_batch_c.mul_assign(delta_vec[3]);

    let mut full_batch_sigma1 = sigma1_commit.0.clone();
    full_batch_sigma1.mul_assign(delta_vec[4]);

    let mut full_batch_sigma2 = sigma2_commit.0.clone();
    full_batch_sigma2.mul_assign(delta_vec[5]);

    let mut full_batched_polynomial_commitment = batched_polynomial_commitment.clone();
    full_batched_polynomial_commitment.add_assign(&full_batch_a);
    full_batched_polynomial_commitment.add_assign(&full_batch_b);
    full_batched_polynomial_commitment.add_assign(&full_batch_c);
    full_batched_polynomial_commitment.add_assign(&full_batch_sigma1);
    full_batched_polynomial_commitment.add_assign(&full_batch_sigma2);
    full_batched_polynomial_commitment.add_assign(&full_batch_quotient_low);

    let mut group_encoded_r = r_open.clone();
    group_encoded_r.mul_assign(&delta_vec[0]);

    let mut group_encoded_a = a_open.clone();
    group_encoded_a.mul_assign(&delta_vec[1]);

    let mut group_encoded_b = b_open.clone();
    group_encoded_b.mul_assign(&delta_vec[2]);

    let mut group_encoded_c = c_open.clone();
    group_encoded_c.mul_assign(&delta_vec[3]);

    let mut group_encoded_sigma1 = sigma1_open.clone();
    group_encoded_sigma1.mul_assign(&delta_vec[4]);

    let mut group_encoded_sigma2 = sigma2_open.clone();
    group_encoded_sigma2.mul_assign(&delta_vec[5]);

    let mut group_encoded_grand_open = grand_product_open.clone();
    group_encoded_grand_open.mul_assign(&upsilon);

    let mut temp_group_encoded = Fr::zero();
    temp_group_encoded.add_assign(&quotient_open);
    temp_group_encoded.add_assign(&group_encoded_r);
    temp_group_encoded.add_assign(&group_encoded_a);
    temp_group_encoded.add_assign(&group_encoded_b);
    temp_group_encoded.add_assign(&group_encoded_c);
    temp_group_encoded.add_assign(&group_encoded_sigma1);
    temp_group_encoded.add_assign(&group_encoded_sigma2);
    temp_group_encoded.add_assign(&group_encoded_grand_open);
    let group_encoded_polynomial_commitment = generator_one.mul(temp_group_encoded);

    let pair1_opening_proof_polynomial = opening_proof_polynomial_commit.0.clone();
    let mut pair1_opening_proof_polynomial_rou = opening_proof_polynomial_rou_commit.0.clone();
    pair1_opening_proof_polynomial_rou.mul_assign(upsilon);
    pair1_opening_proof_polynomial_rou.add_assign(&pair1_opening_proof_polynomial);
    let pairing_dummy_x = commit_maker(&vec![Fr::zero(), Fr::one()], &trusted_setup_g1_vector, &trusted_setup_g2_vector);

    let pair_one = pair1_opening_proof_polynomial_rou.into_affine().pairing_with(&pairing_dummy_x.1.into_affine());
    
    let mut pair2_opening_proof_polynomial = opening_proof_polynomial_commit.0.clone();
    let mut pair2_opening_proof_polynomial_rou = opening_proof_polynomial_rou_commit.0.clone();

    pair2_opening_proof_polynomial.mul_assign(teta);
    let mut teta_ups_rou = teta.clone();
    teta_ups_rou.mul_assign(&upsilon);
    teta_ups_rou.mul_assign(&rou);
    pair2_opening_proof_polynomial_rou.mul_assign(teta_ups_rou);

    pair2_opening_proof_polynomial_rou.add_assign(&pair2_opening_proof_polynomial);
    pair2_opening_proof_polynomial_rou.add_assign(&full_batched_polynomial_commitment);
    pair2_opening_proof_polynomial_rou.sub_assign(&group_encoded_polynomial_commitment);

    let pair_two = pair2_opening_proof_polynomial_rou.into_affine().pairing_with(&generator_two);

    assert_eq!(pair_one, pair_two);

}
