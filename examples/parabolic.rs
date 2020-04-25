use numerical::pde::parabolic::*;
use numerical::explicit::*;
use numerical::implicit::*;
use ndarray::*;
use eom::*;

fn main() {
    let eom = Parabolic::new((0.0, 1.0), 1.0, 0.01, 6);
    let mut teo = Diff::from(eom.clone());
    let ts = adaptor::time_series(Array1::from(
        (0..eom.x_part+1).map(|x| {
            let x = x as f64 * eom.dx + eom.range.0;
            if x <= 0.5 { x }
            else { 1.0 - x }
        }).collect::<Vec<f64>>()
    ), &mut teo);
    let end_time = 10;
    
    println!("*** ∂u/∂t = ∂^2/∂x^2      (0 < x < 1)     ***");
    println!("*** u(x, 0) = {{ x         ( 0 <= x <=1/2) ***");
    println!("***           {{ 1 - x     (1/2<= x <= 1 ) ***");
    println!("*** u(0, t) = u(1, t) = 0 (t >= 0)        ***");
    println!("** Difference scheme **");
    for (_, v) in ts.take(end_time).enumerate() {
        println!("x0 = {:>7.4}, x1 = {:>7.4}, x2 = {:>7.4}, x3 = {:>7.4}, x4 = {:>7.4}, x5 = {:>7.4}, x6 = {:>7.4}",
        v[0], v[1], v[2], v[3], v[4], v[5], v[6]);
    }

    let eom = Parabolic::new((0.0, 4.0), 2.0, 0.1, 8);
    let mut teo = DiffIm::from(eom.clone());
    let ts = adaptor::time_series(Array1::from(
        (0..eom.x_part+1).map(|x| {
            if x == 0 || x == 8 { 0.0 }
            else { 1.0 }
        }).collect::<Vec<f64>>()
    ), &mut teo);
    let end_time = 10;
    
    println!("*** ∂u/∂t = 2∂^2/∂x^2     (0 < x < 4)     ***");
    println!("*** u(x, 0) = {{ 0         (x = 0, 4) ***");
    println!("***           {{ 1         (else) ***");
    println!("** Difference scheme (Implicit)**");
    for (_, v) in ts.take(end_time).enumerate() {
        println!("x0 = {:>7.4}, x1 = {:>7.4}, x2 = {:>7.4}, x3 = {:>7.4}, x4 = {:>7.4}, x5 = {:>7.4}, x6 = {:>7.4}, x7 = {:>7.4}, x8 = {:>7.4}",
        v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
    }
}