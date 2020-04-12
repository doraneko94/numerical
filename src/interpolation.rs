use num_traits::float::Float;
use num_traits::{Zero, One};

pub trait Itp<F: Float> {
    fn new(x: &Vec<F>, y: &Vec<F>) -> Self;
    fn push(&mut self, x: F, y: F);
    fn push_vec(&mut self, x: &Vec<F>, y: &Vec<F>);
    fn calc(&self, x: F) -> F;
}

pub struct LagrangeItp<F: Float> {
    pub x: Vec<F>,
    pub y: Vec<F>,
}

impl<F: Float> Itp<F> for LagrangeItp<F> {
    fn new(x: &Vec<F>, y: &Vec<F>) -> Self {
        if x.len() != y.len() {
            panic!("x and y have different number of elements!");
        }
        Self { x: x.clone(), y: y.clone() }
    }

    fn push(&mut self, x: F, y: F) {
        self.x.push(x);
        self.y.push(y);
    }

    fn push_vec(&mut self, x: &Vec<F>, y: &Vec<F>) {
        if x.len() != y.len() {
            panic!("x and y have different number of elements!");
        }
        for (&xi, &yi) in x.iter().zip(y.iter()) {
            self.x.push(xi);
            self.y.push(yi);
        }
    }

    fn calc(&self, z: F) -> F {
        let mut ans: F = Zero::zero();
        for (i, (&xi, &yi)) in self.x.iter().zip(self.y.iter()).enumerate() {
            let mut l: F = One::one();
            for (j, &xj) in self.x.iter().enumerate() {
                if i != j {
                    l = l * (z - xj) / (xi - xj);
                }
            }
            ans = ans + l * yi;
        }
        ans
    }
}

pub struct NewtonDivItp<F: Float> {
    pub x: Vec<F>,
    pub y: Vec<F>,
    pub diff: Vec<Vec<F>>,
}

impl<F: Float> Itp<F> for NewtonDivItp<F> {
    fn new(x: &Vec<F>, y: &Vec<F>) -> Self {
        let n = x.len();
        if n != y.len() {
            panic!("x and y have different number of elements!");
        }
        let x = x.clone();
        let y = y.clone();
        let mut diff: Vec<Vec<F>> = vec![(0..n-1).map(|i| (y[i] - y[i+1]) / (x[i] - x[i+1])).collect()];
        for i in 1..n-1 {
            diff.push((0..n-i-1).map(|j|  (diff[i-1][j] - diff[i-1][j+1]) / (x[j] - x[j+i+1])).collect());
        }
        Self { x, y, diff }
    }

    fn push(&mut self, x: F, y: F) {
        let n = self.x.len();
        self.x.push(x);
        self.y.push(y);
        self.diff[0].push((self.y[n-1] - self.y[n]) / (self.x[n-1] - self.x[n]));
        for i in 1..n-1 {
            let d_diff = self.diff[i-1][n-1-i] - self.diff[i-1][n-i];
            self.diff[i].push(d_diff / (self.x[n-i-1] - self.x[n]));
        }
        let d_diff = self.diff[n-2][0] - self.diff[n-2][1];
        self.diff.push(vec![d_diff / (self.x[0] - self.x[n])]);
    }

    fn push_vec(&mut self, x: &Vec<F>, y: &Vec<F>) {
        if x.len() != y.len() {
            panic!("x and y have different number of elements!");
        }
        for (&xi, &yi) in x.iter().zip(y.iter()) {
            self.push(xi, yi);
        }
    }

    fn calc(&self, z: F) -> F {
        let mut ans: F = self.y[0];
        let mut m: F = One::one();
        for i in 0..self.diff.len() {
            m = m * (z - self.x[i]);
            ans = ans + m * self.diff[i][0];
        }
        ans
    }
}

pub fn pnx<F: Float>(n: usize, x: F) -> F {
    let one: F = One::one();
    if n == 0 {
        one
    } else if n == 1 {
        x
    } else {
        let mut p_old = one;
        let mut p_now = x;
        let mut p_new = one;
        for i in 1..n {
            p_new = (F::from(2 * i + 1).unwrap() * x * p_now - F::from(i).unwrap() * p_old) / F::from(i + 1).unwrap();
            p_old = p_now;
            p_now = p_new;
        }
        p_new
    }
}