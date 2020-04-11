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