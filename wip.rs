 {
        let mut carrier: Vec<Fr> = vec![Fr::one(); q_left_lag.len()];
        let mut finders = Fr::one();
        let mut shit_gamma = Fr::one();
        let mut shit_beta = Fr::one();
        shit_gamma.add_assign(&Fr::one());
        for i in 0..11{
        shit_gamma.add_assign(&Fr::one());
        shit_beta.add_assign(&Fr::one());
        }
        
        for i in 0..q_left_lag.len(){
        let mut temp_a = Fr::zero();
        let mut temp_b = Fr::zero();
        let mut temp_c = Fr::zero();
        let mut accumulator_carrier = Fr::one();
        let aw = Polynomial{ coefficients: vec![1, 4, 16, 13] };
        let baw = Polynomial{ coefficients: vec![2, 8, 15, 9] };
        let caw = Polynomial{ coefficients: vec![3, 12, 14, 5] };

        let aaw = Polynomial{ coefficients: vec![2, 8, 15, 3] };
        let baaw = Polynomial{ coefficients: vec![1, 4, 16, 12] };
        let caaw = Polynomial{ coefficients: vec![13, 9, 5, 14] };

        let shittest_rou = aw.to_field();
        let mut shit_rou1 = aaw.to_field();
        let mut shit_k11 = baaw.to_field();
        let mut shit_k21 = caaw.to_field();

        let mut shit_rou = aw.to_field();
        let mut shit_k1 = baw.to_field();
        let mut shit_k2 = caw.to_field();

        for j in 0..i{
            shit_rou = aw.to_field();
            shit_k1 = baw.to_field();
            shit_k2 = caw.to_field();
            //(wj+beta*rou+gamma)
            //Numerator
            temp_a = eval_poly_fr(&a_lag, &rou_vec_temp2[j]);
            shit_rou[j].mul_assign(&shit_beta);
            temp_a.add_assign(&shit_rou[j]);
            temp_a.add_assign(&shit_gamma);
            println!("{:#?}qqq", temp_a);
            
            temp_b = eval_poly_fr(&b_lag, &rou_vec_temp2[j]);
            shit_k1[j].mul_assign(&shit_beta);
            temp_b.add_assign(&shit_k1[j]);
            temp_b.add_assign(&shit_gamma);
            println!("{:#?}eee", temp_b);


            temp_c = eval_poly_fr(&c_lag, &rou_vec_temp2[j]);
            shit_k2[j].mul_assign(&shit_beta);
            temp_c.add_assign(&shit_k2[j]);
            temp_c.add_assign(&shit_gamma);
            println!("{:#?}rrr", temp_c);

    
            accumulator_carrier.mul_assign(&temp_a);
            accumulator_carrier.mul_assign(&temp_b);
            accumulator_carrier.mul_assign(&temp_c);
    
            //Denominator
            shit_rou1 = aaw.to_field();
            shit_k11 = baaw.to_field();
            shit_k21 = caaw.to_field();
    
            temp_a = eval_poly_fr(&a_lag, &rou_vec_temp2[j]);
            shit_rou1[j].mul_assign(&shit_beta);
            temp_a.add_assign(&shit_rou1[j]);
            temp_a.add_assign(&shit_gamma);
            println!("{:#?}ccc", temp_a);

    
            temp_b = eval_poly_fr(&b_lag, &rou_vec_temp2[j]);
            shit_k11[j].mul_assign(&shit_beta);
            temp_b.add_assign(&shit_k11[j]);
            temp_b.add_assign(&shit_gamma);
            println!("{:#?}bbb", temp_b);
    
            temp_c = eval_poly_fr(&c_lag, &rou_vec_temp2[j]);
            shit_k21[j].mul_assign(&shit_beta);
            temp_c.add_assign(&shit_k21[j]);
            temp_c.add_assign(&shit_gamma);
            println!("{:#?}ddd", temp_c);

    
            accumulator_carrier.mul_assign(&temp_a.inverse().unwrap());
            accumulator_carrier.mul_assign(&temp_b.inverse().unwrap());
            accumulator_carrier.mul_assign(&temp_c.inverse().unwrap());
            println!("{:#?}aa", accumulator_carrier);

            }
            //Here adding every new polynomial to carrier to keep them as a grand product
            finders.mul_assign(&accumulator_carrier);
            carrier[i] = finders;
        }
        println!("{:#?}", carrier);

    }
