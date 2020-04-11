use numerical::integral::*;

fn main() {
    let i: Integral<f64> = Integral::new(|x: f64| (1.0 - x) * (-x).exp());
    println!("n = 10: I = {}", i.trapezoid(-1.0, 1.0, 10));
    println!("n = 20: I = {}", i.trapezoid(-1.0, 1.0, 20));
    println!("n = 30: I = {}", i.trapezoid(-1.0, 1.0, 30));

    println!("n = 5: I = {}", i.simpson(-1.0, 1.0, 10));
}