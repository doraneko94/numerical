use num_traits::float::Float;
use num_traits::{Zero, One};

pub struct Integral<F: Float> {
    pub rhs: fn(F) -> F,
}

impl<F:Float> Integral<F> {
    pub fn new(rhs: fn(F)->F) -> Self {
        Self { rhs }
    }

    pub fn trapezoid(&self, a: F, b: F, n: usize) -> F {
        if a == b {
            return Zero::zero();
        }
        let (a, b, reverse) = if a > b {
            (b, a, true)
        } else {
            (a, b, false)
        };
        let dh = (b - a) / F::from(n).unwrap();
        let zero: F = Zero::zero();
        let one: F = One::one();
        let two: F = one + one;
        let s1 = ((self.rhs)(a) + (self.rhs)(b)) / two;
        let s2 = (1..n).map(|i| (self.rhs)(a + F::from(i).unwrap() * dh)).fold(zero, |m, j| m + j);
        if reverse {
            zero - (s1 + s2) * dh
        } else {
            (s1 + s2) * dh
        }
    }

    pub fn simpson(&self, a: F, b: F, n: usize) -> F {
        if a == b {
            return Zero::zero();
        }
        let (a, b, reverse) = if a > b {
            (b, a, true)
        } else {
            (a, b, false)
        };
        let n2 = n + n;
        let dh = (b - a) / F::from(n2).unwrap();
        let zero: F = Zero::zero();
        let one: F = One::one();
        let two: F = one + one;
        let three: F = two + one;
        let four: F = two + two;
        let f: Vec<F> = (0..n2+1).map(|i| (self.rhs)(a + F::from(i).unwrap() * dh)).collect();
        let s1 = f[0] + f[n2];
        let s2 = (1..n).map(|i| f[2 * i] * two).fold(zero, |m, j| m + j);
        let s4 = (0..n).map(|i| f[2 * i + 1] * four).fold(zero, |m, j| m + j);
        if reverse {
            zero - (s1 + s2 + s4) * dh / three
        } else {
            (s1 + s2 + s4) * dh / three
        }
    }
}