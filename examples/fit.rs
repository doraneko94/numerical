use numerical::fit::*;

fn main() {
    println!("*** 3D spline curve fitting ***");
    let sp: Spline3d<f64> = Spline3d::new(&vec![0.0, 1.0, 1.5, 2.0, 3.0],
                                          &vec![2.0, 4.0, 3.0, 1.0, 2.0],
                                          5.0, 3.0);
    println!("x = [0.0, 1.0, 1.5, 2.0, 3.0]");
    println!("y = [2.0, 4.0, 3.0, 1.0, 2.0]");
    println!("C0 = 5, Cn = 3");
    println!("** result **");
    println!("f(1.75) = {:.5}", sp.calc(1.75).unwrap());
    println!("");

    println!("*** Least-square method fitting ***");
    let lsm: LSM<f64> = LSM::new(&vec![ 0.2,  0.5,  1.0,  2.0,  4.0,  8.0, 10.0],
                                 &vec![12.1,  4.9,  2.9,  2.1,  2.1,  3.4,  4.3],
                                 vec![|x: f64| x, |x: f64| 1.0 / x]);
    println!("x = [ 0.2,  0.5,  1.0,  2.0,  4.0,  8.0, 10.0]");
    println!("y = [12.1,  4.9,  2.9,  2.1,  2.1,  3.4,  4.3]");
    println!("** equation **");
    println!("y = ax + (b/x)");
    println!("** result **");
    println!("a = {:.5}", lsm.coeff[0]);
    println!("b = {:.5}", lsm.coeff[1]);
    println!("f(6) = {:.5}", lsm.calc(6.0));
}