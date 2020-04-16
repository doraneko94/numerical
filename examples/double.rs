use numerical::double::*;

fn main() {
    let i = DIntegral::<f64>::new(|x: f64, y:f64| x + y * y);
    println!("*** f(x) = x - y^2 ***");
    println!("*** I = \\int\\int_D f(x) dxdy ***");
    println!("*** D: 0 <= x <= 2, 1 <= y <= 1 + 0.5x ***\n");

    println!("** trapezoid (n = 20) **");
    println!("n = 20: I = {:.6}\n", i.trapezoid((0.0, 2.0), (|_x: f64| 1.0, |x: f64| 1.0 + 0.5 * x), 20));

    println!("** Simpson (n = 20) **");
    println!("n = 20: I = {:.6}\n", i.simpson((0.0, 2.0), (|_x: f64| 1.0, |x: f64| 1.0 + 0.5 * x), 20));

    println!("** Gauss-Legendre (n = 10) **");
    println!("n = 10: I = {:.6}\n", i.gauss_legendre((0.0, 2.0), (|_x: f64| 1.0, |x: f64| 1.0 + 0.5 * x), 10));

    println!("** Chebyshev-Horinouchi (n = 10) **");
    println!("n = 10: I = {:.6}\n", i.chebyshev((0.0, 2.0), (|_x: f64| 1.0, |x: f64| 1.0 + 0.5 * x), 10));

    println!("** Double Exponential Formula (n = 20, h = 0.1) **");
    println!("n = 20: I = {:.6}\n", i.def((0.0, 2.0), (|_x: f64| 1.0, |x: f64| 1.0 + 0.5 * x), 20, 0.1));
}