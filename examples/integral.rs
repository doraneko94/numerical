use numerical::integral::*;

fn main() {
    let i: Integral<f64> = Integral::new(|x: f64| (1.0 - x*x).sqrt());
    println!("*** f(x) = (1 - x^2).sqrt() ***");
    println!("*** I = \\int_0^1 f(x) dx ***\n");

    println!("** trapezoid (n = 20) **");
    println!("n = 20: I = {:.6}\n", i.trapezoid(0.0, 1.0, 20));

    println!("** Simpson (n = 20) **");
    println!("n = 20: I = {:.6}\n", i.simpson(0.0, 1.0, 20));

    println!("** Gauss-Legendre (n = 10) **");
    println!("n = 10: I = {:.6}\n", i.gauss_legendre(0.0, 1.0, 10));

    println!("** Chebyshev-Horinouchi (n = 10) **");
    println!("n = 10: I = {:.6}\n", i.chebyshev(0.0, 1.0, 10));

    println!("** Double Exponential Formula (n = 20, h = 0.1) **");
    println!("n = 20: I = {:.6}", i.def(0.0, 1.0, 20, 0.1));
}