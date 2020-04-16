use eom::traits::*;
use eom::*;
use ndarray::*;

#[derive(Clone, Copy, Debug)]
struct ODETest {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

impl Default for ODETest {
    fn default() -> Self {
        ODETest {
            a: 1.0,
            b: -12.0,
            c: 3.0,
        }
    }
}

impl ModelSpec for ODETest {
    type Scalar = f64;
    type Dim = Ix1;
    fn model_size(&self) -> usize {
        2
    }
}

impl Explicit for ODETest {
    fn rhs<'a, S>(&mut self, v: &'a mut ArrayBase<S, Ix1>) -> &'a mut ArrayBase<S, Ix1>
    where
        S: DataMut<Elem = f64>,
    {
        let x = v[0];
        let y = v[1];
        v[0] = 1.0;
        v[1] = self.a * y + self.b * x + self.c;
        v
    }
}

#[derive(Clone, Copy, Debug)]
struct SimODETest {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

impl Default for SimODETest {
    fn default() -> Self {
        SimODETest {
            a: 1.0,
            b: -1.0,
            c: 1.0,
        }
    }
}

impl ModelSpec for SimODETest {
    type Scalar = f64;
    type Dim = Ix1;
    fn model_size(&self) -> usize {
        3
    }
}

impl Explicit for SimODETest {
    fn rhs<'a, S>(&mut self, v: &'a mut ArrayBase<S, Ix1>) -> &'a mut ArrayBase<S, Ix1>
    where
        S: DataMut<Elem = f64>,
    {
        let x = v[0];
        let y = v[1];
        let z = v[2];
        v[0] = 1.0;
        v[1] = self.c * z;
        v[2] = self.b * x * (self.a * y - self.c * z);
        v
    }
}

fn main() {
    let h = 0.1;
    let eom = ODETest::default();
    let mut teo = explicit::Euler::new(eom, h);
    let ts = adaptor::time_series(arr1(&[0.0, 1.0]), &mut teo);
    let end_time = 10;

    let ans = |x: f64| 12.0 * x - 8.0 * x.exp() + 9.0;
    
    println!("***      dy/dx = y(x) - 12x + 3      ***");
    println!("*** answer: y(x) = 12x - 8exp(x) + 9 ***\n");
    println!("** Eular scheme **");
    println!("x = 0.0, y =  1.0000, truth =  1.0000");
    for (_, v) in ts.take(end_time).enumerate() {
        println!("x = {:.1}, y = {:>7.4}, truth = {:>7.4}", v[0], v[1], ans(v[0]));
    }

    let mut teo = explicit::RK4::new(eom, h);
    let ts = adaptor::time_series(arr1(&[0.0, 1.0]), &mut teo);

    println!("\n** Runge-Kutta 4 scheme **");
    for (_, v) in ts.take(end_time).enumerate() {
        println!("x = {:.1}, y = {:>7.4}, truth = {:>7.4}", v[0], v[1], ans(v[0]));
    }

    let eom = SimODETest::default();
    let mut teo = explicit::RK4::new(eom, h);
    let ts = adaptor::time_series(arr1(&[0.0, 1.0, 1.0]), &mut teo);

    println!("\n*** d^2y/dx^2 - x dy/dx + xy = 0 ***");
    println!("***         z(x) = dy/dx         ***");
    println!("***      {{ dy/dx = z             ***");
    println!("***      {{ dz/dx = x(z - y)      ***\n");

    println!("** Runge-Kutta 4 scheme **");
    println!("x = 0.0, y =  1.0000, z =  1.0000");
    for (_, v) in ts.take(end_time).enumerate() {
        println!("x = {:.1}, y = {:>7.4}, z = {:>7.4}", v[0], v[1], v[2]);
    }
}