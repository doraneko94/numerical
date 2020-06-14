use numerical::fem::poisson1d::*;
use ndarray::*;

fn main() {
    println!("*** d^2u/dx^2 = 1 (-1 <= x <= 1) ***");
    println!("** u(0) = 1, du/dx(1) = -1 **");
    let p = Poisson1dFEM::<f64>::default();
    let x_arr: Array1<f64> = Array::from((0..4)
                                        .map(|i| -1.0 + 2.0 / 3.0 * i as f64)
                                        .collect::<Vec<f64>>());
    let u_arr = p.calc().unwrap();
    for (x, u) in x_arr.iter().zip(u_arr.iter()) {
        println!("x = {:6.3}, u = {:6.3}", x, u);
    }
}