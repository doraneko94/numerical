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
}