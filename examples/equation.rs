use numerical::equation::*;

fn main() {
    let eq = Equation::new(|x: f64| x*x*x - 3.0*x*x + 9.0*x - 8.0, 1e-7);
    println!("*** f(x) = x^3 - 3x^2 + 9x - 8 ***");
    println!("f(1) = {:.7}", eq.calc(1.0));
    println!("** bisection **");
    println!("x1 = {:.7}", eq.bisection(1.0, 2.0).unwrap());
    println!("** newton **");
    println!("x1 = {:.7}\n", eq.newton(1.0).unwrap());

    let eq = Equation::new(|x: f64| x*x + 1_000_000_000_000_000.0*x + 100_000_000_000_000.0, 1e-7);
    println!("*** f(x) = x^2 - 10^15x + 10^14 ***");
    println!("f(1) = {:.7}", eq.calc(1.0));
    println!("** bisection **");
    println!("x1 = {:.7}", eq.bisection(-1.0, 1.0).unwrap());
    println!("** newton **");
    println!("x1 = {:.7}\n", eq.newton(1.0).unwrap());

    let eq = Equation::new(|x: f64| x*x*x*x - 3.0*x + 1.0, 1e-7);
    println!("*** f(x) = x^4 - 3x + 1 ***");
    println!("** bisection **");
    println!("x = {:.7}, {:.7}", eq.bisection(0.0, 1.0).unwrap(), eq.bisection(1.0, 2.0).unwrap());
    println!("** newton **");
    println!("x = {:.7}, {:.7}\n", eq.newton(0.0).unwrap(), eq.newton(1.5).unwrap());

    let eq = Equation::new(|x: f64| x.cos() - x, 1e-7);
    println!("*** f(x) = cosx - x ***");
    println!("** bisection **");
    println!("x = {:.7}", eq.bisection(0.0, 1.0).unwrap());
    println!("** newton **");
    println!("x = {:.7}\n", eq.newton(0.0).unwrap());

    let eq = Equation::new(|x: f64| x.exp() - 1.0 / x, 1e-7);
    println!("*** f(x) = exp(x) - 1/x ***");
    println!("** bisection **");
    println!("x = {:.7}", eq.bisection(0.0, 1.0).unwrap());
    println!("** newton **");
    println!("x = {:.7}\n", eq.newton(1.0).unwrap());
}