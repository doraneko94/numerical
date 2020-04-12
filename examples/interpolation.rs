use numerical::interpolation::*;

fn main() {
    println!("*** Lagrange interpolation ***");
    let i: LagrangeItp<f64> = LagrangeItp::new(&vec![   0.5,    1.0,    1.5,    2.0], 
                                               &vec![0.3734, 0.5104, 0.4712, 0.3345]);
    println!("x = [   0.5,    1.0,    1.5,    2.0]");
    println!("y = [0.3734, 0.5104, 0.4712, 0.3345]");
    println!("f(0.5) = {:.5}", i.calc(0.5));
    println!("f(1.0) = {:.5}", i.calc(1.0));
    println!("f(1.5) = {:.5}", i.calc(1.5));
    println!("f(2.0) = {:.5}", i.calc(2.0));
    println!("f(0.8) = {:.5}", i.calc(0.8));
    println!("");

    let mut i: LagrangeItp<f64> = LagrangeItp::new(&vec![   0.5,    1.5], 
                                                   &vec![0.3734, 0.4712]);
    i.push_vec(&vec![   1.0,    2.0], 
               &vec![0.5104, 0.3345]);
    println!("x = [   0.5,    1.5,    1.0,    2.0]");
    println!("y = [0.3734, 0.4712, 0.5104, 0.3345]");
    println!("f(0.5) = {:.5}", i.calc(0.5));
    println!("f(1.0) = {:.5}", i.calc(1.0));
    println!("f(1.5) = {:.5}", i.calc(1.5));
    println!("f(2.0) = {:.5}", i.calc(2.0));
    println!("f(0.8) = {:.5}", i.calc(0.8));
    println!("");

    println!("*** Newton's Divided Difference interpolation ***");
    let i: NewtonDivItp<f64> = NewtonDivItp::new(&vec![   0.2,    0.5,    1.0,    1.5,    2.0,    3.0], 
                                                 &vec![0.0793, 0.1915, 0.3413, 0.4332, 0.4772, 0.4987]);
    println!("x = [   0.2,    0.5,    1.0,    1.5,    2.0,    3.0]");
    println!("y = [0.0793, 0.1915, 0.3413, 0.4332, 0.4772, 0.4987]");
    println!("f(0.7) = {:.5}", i.calc(0.7));
    println!("");

    let mut i: NewtonDivItp<f64> = NewtonDivItp::new(&vec![   0.2,    1.0,    2.0], 
                                                     &vec![0.0793, 0.3413, 0.4772]);
    i.push_vec(&vec![   0.5,    1.5,    3.0], 
               &vec![0.1915, 0.4332, 0.4987]);
    println!("x = [   0.2,    1.0,    2.0,    0.5,    1.5,    3.0]");
    println!("y = [0.0793, 0.3413, 0.4772, 0.1915, 0.4332, 0.4987]");
    println!("f(0.7) = {:.5}", i.calc(0.7));
}