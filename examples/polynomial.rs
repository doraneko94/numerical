use numerical::polynomial::*;

fn main() {
    let chev: Chebyshev<f64> = Chebyshev::new(|x: f64| x.exp(), 5);
    println!("*** Interpolation by Chebyshev polynomial ***");
    println!("y = exp(x) [-1, 1]; ~T5(x)");
    println!("{:?}", chev.c);
    println!("exp(0.0) = {:.6}", chev.calc(0.0).unwrap());
    println!("exp(0.5) = {:.6}", chev.calc(0.5).unwrap());
    println!("exp(1.0) = {:.6}", chev.calc(1.0).unwrap());
}