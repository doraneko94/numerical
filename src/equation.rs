use num_traits::float::Float;
use num_traits::{Zero, One};

pub struct Equation<F: Float> {
    pub rhs: fn(F) -> F,
    pub dh: F,
}

impl<F: Float> Equation<F> {

    pub fn new(rhs: fn(F) -> F, dh: F) -> Self {
        Self { rhs, dh, }
    }

    pub fn calc(&self, x: F) -> F {
        (self.rhs)(x)
    }

    pub fn bisection(&self, mut inf: F, mut sup: F) -> Option<F> {
        if inf > sup {
            return None;
        }
        let f_inf = self.calc(inf);
        let f_sup = self.calc(sup);
        if f_inf.abs() < self.dh {
            return Some(inf);
        }
        if f_sup.abs() < self.dh {
            return Some(sup);
        }
        if f_inf * f_sup > Zero::zero() {
            return None;
        }
        let rev = if f_inf < Zero::zero() {
            false
        } else {
            true
        };
        let one: F = One::one();
        let two = one + one;
        while sup - inf >= self.dh {
            let mid = (sup + inf) / two;
            let f_mid = self.calc(mid);
            if f_mid < Zero::zero() {
                if rev { sup = mid; }
                else { inf = mid; }
            } else {
                if rev { inf = mid; }
                else { sup = mid; }
            }
        }
        Some(inf)
    }

    pub fn newton(&self, mut x: F) -> Option<F> {
        let mut dx = Float::infinity();
        while self.calc(x).abs() >= self.dh {
            let f_x = self.calc(x);
            let grad = (self.calc(x + self.dh) - f_x) / self.dh;
            let next_x = x - f_x / grad;
            if (next_x - x).abs() > dx {
                return None;
            }
            dx = (next_x - x).abs();
            x = next_x;
        }
        Some(x)
    }
}