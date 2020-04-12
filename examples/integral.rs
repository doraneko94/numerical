use numerical::integral::*;

fn main() {
    let i: Integral<f64> = Integral::new(|x: f64| (1.0 - x) * (-x).exp());
    println!("*** f(x) = (1 - x) * e^(-x) ***");
    println!("** trapezoid [-1, 1], n = 10 **");
    println!("n = 10: I = {:.8}", i.trapezoid(-1.0, 1.0, 10));
    println!("n = 20: I = {:.8}", i.trapezoid(-1.0, 1.0, 20));
    println!("n = 30: I = {:.8}\n", i.trapezoid(-1.0, 1.0, 30));

    println!("** Simpson [-1, 1], n = 10 **");
    println!("n = 5: I = {:.8}\n", i.simpson(-1.0, 1.0, 10));

    println!("** Gauss-Legendre [-1, 1], n = 10 **");
    println!("n = 10: I = {:.8}", i.gauss_legendre(-2.0, 1.0, 10));
}