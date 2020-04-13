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
            diff.push((0..n-i-1).map(|j| (diff[i-1][j] - diff[i-1][j+1]) / (x[j] - x[j+i+1])).collect());
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

pub struct NewtonForItp<F: Float> {
    pub x: Vec<F>,
    pub y: Vec<F>,
    pub diff: Vec<Vec<F>>,
}

impl<F: Float> Itp<F> for NewtonForItp<F> {
    fn new(x: &Vec<F>, y: &Vec<F>) -> Self {
        let n = x.len();
        if n != y.len() {
            panic!("x and y have different number of elements!");
        }
        if n < 2 {
            panic!("length of x should be more than 2!");
        }
        let dh = x[1] - x[0];
        let epsilon = F::from(1e-5f64).unwrap();
        for i in 1..n-1 {
            if (x[i+1] - x[i] - dh).abs() > epsilon {
                panic!("x should be equally spaced!");
            }
        }
        let x = x.clone();
        let y = y.clone();
        let mut diff: Vec<Vec<F>> = vec![(0..n-1).map(|i| y[i+1] - y[i]).collect()];
        for i in 1..n-1 {
            diff.push((0..n-i-1).map(|j| diff[i-1][j+1] - diff[i-1][j]).collect());
        }
        Self { x, y, diff }
    }

    fn push(&mut self, x: F, y: F) {
        let dh = self.x[1] - self.x[0];
        let epsilon = F::from(1e-5f64).unwrap();
        let n = self.x.len();
        if (dh - (x - self.x[n-1])).abs() > epsilon {
            panic!("x should be equally spaced!");
        } 
        self.x.push(x);
        self.y.push(y);
        self.diff[0].push(self.y[n] - self.y[n-1]);
        for i in 1..n-1 {
            let d_diff = self.diff[i-1][n-i] - self.diff[i-1][n-i-1];
            self.diff[i].push(d_diff);
        }
        let d_diff = self.diff[n-2][1] - self.diff[n-2][0];
        self.diff.push(vec![d_diff]);
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
        let one: F = One::one();
        let mut ans: F = self.y[0];
        let mut c: F = one;
        let dh = self.x[1] - self.x[0];
        let k = (z - self.x[0]) / dh;
        for i in 0..self.diff.len() {
            let i_f = F::from(i).unwrap();
            c = c * (k - i_f) / (i_f + one);
            ans = ans + c * self.diff[i][0];
        }
        ans
    }
}