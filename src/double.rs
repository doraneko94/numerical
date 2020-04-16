use num_traits::float::{Float, FloatConst};
use num_traits::{Zero, One};
use crate::polynomial::{zeta_vec};
use crate::integral::*;

pub struct DIntegral<F: Float + FloatConst> {
    pub rhs: fn(F, F) -> F,
}

impl<F:Float + FloatConst> DIntegral<F> {
    pub fn new(rhs: fn(F, F)->F) -> Self {
        Self { rhs }
    }
    
    pub fn trapezoid(&self, ab1: (F, F), ab2: (fn(F)->F, fn(F)->F), n: usize) -> F {
        let (a1, b1) = ab1;
        if a1 == b1 {
            return Zero::zero();
        }
        let (a1, b1, reverse) = if a1 > b1 {
            (b1, a1, true)
        } else {
            (a1, b1, false)
        };
        let dh1 = (b1 - a1) / F::from(n).unwrap();
        let zero: F = Zero::zero();
        let one: F = One::one();
        let two: F = one + one;
        let ia = {
            let (a2, b2) = ((ab2.0)(a1), (ab2.1)(a1));
            _trapezoid(self.rhs, a1, a2, b2, n)
        };
        let ib = {
            let (a2, b2) = ((ab2.0)(b1), (ab2.1)(b1));
            _trapezoid(self.rhs, b1, a2, b2, n)
        };
        let s1 = (ia + ib) / two;
        let s2 = (1..n).map(|i| {
            let ai = a1 + dh1 * F::from(i).unwrap();
            let (a2, b2) = ((ab2.0)(ai), (ab2.1)(ai));
            _trapezoid(self.rhs, ai, a2, b2, n)
        }).fold(zero, |m, j| m + j);
        if reverse {
            zero - (s1 + s2) * dh1
        } else {
            (s1 + s2) * dh1
        }
    }

    pub fn simpson(&self, ab1: (F, F), ab2: (fn(F)->F, fn(F)->F), n: usize) -> F {
        let (a1, b1) = ab1;
        if a1 == b1 {
            return Zero::zero();
        }
        let (a1, b1, reverse) = if a1 > b1 {
            (b1, a1, true)
        } else {
            (a1, b1, false)
        };
        let n2 = n + n;
        let dh1 = (b1 - a1) / F::from(n2).unwrap();
        let zero: F = Zero::zero();
        let one: F = One::one();
        let two: F = one + one;
        let three: F = two + one;
        let four: F = two + two;
        let f1: Vec<F> = (0..n2+1).map(|i| {
            let ai = a1 + F::from(i).unwrap() * dh1;
            let (a2, b2) = ((ab2.0)(ai), (ab2.1)(ai));
            if a2 == b2 {
                return zero;
            }
            let (a2, b2, reverse) = if a2 > b2 {
                (b2, a2, true)
            } else {
                (a2, b2, false)
            };
            let dh2 = (b2 - a2) / F::from(n2).unwrap();
            let f2: Vec<F> = (0..n2+1).map(|i| (self.rhs)(ai, a2 + F::from(i).unwrap() * dh2)).collect();
            let s1 = f2[0] + f2[n2];
            let s2 = (1..n).map(|i| f2[2 * i] * two).fold(zero, |m, j| m + j);
            let s4 = (0..n).map(|i| f2[2 * i + 1] * four).fold(zero, |m, j| m + j);
            if reverse {
                zero - (s1 + s2 + s4) * dh2 / three
            } else {
                (s1 + s2 + s4) * dh2 / three
            }
        }).collect();
        let s1 = f1[0] + f1[n2];
        let s2 = (1..n).map(|i| f1[2 * i] * two).fold(zero, |m, j| m + j);
        let s4 = (0..n).map(|i| f1[2 * i + 1] * four).fold(zero, |m, j| m + j);
        if reverse {
            zero - (s1 + s2 + s4) * dh1 / three
        } else {
            (s1 + s2 + s4) * dh1 / three
        }
    }

    pub fn gauss_legendre(&self, ab1: (F, F), ab2: (fn(F)->F, fn(F)->F), n: usize) -> F {
        let alpha: Vec<F> = solve_pnx(n);
        let w: Vec<F> = calc_weights(&alpha);
        let zero: F = Zero::zero();
        let one: F = One::one();
        let two: F = one + one;
        let (a1, b1) = ab1;
        let tau = |t: F, a: F, b: F| ((b - a) * t  + (b + a)) / two;
        let s1 = w.iter().zip(alpha.iter()).map(|(&w1i, &a1i)| {
            let x1i = tau(a1i, a1, b1);
            let (a2, b2) = ((ab2.0)(x1i), (ab2.1)(x1i));
            let s2 = w.iter().zip(alpha.iter()).map(|(&w2i, &a2i)| {
                let x2i = tau(a2i, a2, b2);
                w2i * (self.rhs)(x1i, x2i)
            }).fold(zero, |m, e| m + e);
            w1i * s2 * (b2 - a2)
        }).fold(zero, |m, e| m + e);
        (b1 - a1) * s1 / (two * two)
    }

    pub fn chebyshev(&self, ab1: (F, F), ab2: (fn(F)->F, fn(F)->F), n: usize) -> F {
        let (a1, b1) = ab1;
        let n2 = n / 2;
        let n1f: F = F::from(2*(n+1)).unwrap();
        let theta: Vec<F> = (0..n+1).map(|k: usize| F::from(2*k+1).unwrap() / n1f * F::PI()).collect();
        let zeta = zeta_vec(n);
        let zero: F = Zero::zero();
        let one: F = One::one();
        let two: F = one + one;
        let tau = |t: F, a: F, b: F| ((b - a) * t + (b + a)) / two;
        let denom: Vec<F> = (1..n2+1).map(|j| F::from(2*j*2*j-1).unwrap()).collect();
        let n1_hnk = |k: usize| {
            let s: F = (1..n2+1).map(|j| (F::from(2 * j).unwrap() * theta[k]).cos() / denom[j-1]).fold(zero, |m, i| m + i);
            one - two * s
        };
        let s1 = (0..n+1).map(|k| {
            let z1 = tau(zeta[k], a1, b1);
            let (a2, b2) = ((ab2.0)(z1), (ab2.1)(z1));
            let s2 = (0..n+1).map(|l| {
                let z2 = tau(zeta[l], a2, b2);
                n1_hnk(l) * (self.rhs)(z1, z2)
            }).fold(zero, |m, i| m + i);
            n1_hnk(k) * s2 * (b2 - a2)}).fold(zero, |m, i| m + i);
        let n1 = F::from(n + 1).unwrap();
        (b1 - a1) / (n1 * n1) * s1
    }

    pub fn def(&self, ab1: (F, F), ab2: (fn(F)->F, fn(F)->F), n: usize, h: F) -> F {
        let (a1, b1) = ab1;
        let zero: F = Zero::zero();
        let one: F = One::one();
        let two: F = one + one;
        let nf: F = F::from(n).unwrap();
        let tau = |t: F, a: F, b: F| ((b - a) * t + (b + a)) / two;
        let pi: F = F::PI();
        let phi = |t: F| (pi / two * t.sinh()).tanh();
        let phip = |t: F| pi * t.cosh() / (one + (pi * t.sinh()).cosh());
        let tn: F = h * nf;
        let t0: F = zero - tn;
        let gtau = |t: F, a: F, b: F| {
            let z1 = tau(phi(t), a, b);
            let (a2, b2) = ((ab2.0)(z1), (ab2.1)(z1));
            let mut s2 = ((self.rhs)(z1, tau(phi(t0), a2, b2)) * phip(t0) + (self.rhs)(z1, tau(phi(tn), a2, b2)) * phip(tn)) / two;
            for i in 1..2*n {
                let t: F = (F::from(i).unwrap() - nf) * h;
                s2 = s2 + (self.rhs)(z1, tau(phi(t), a2, b2)) * phip(t);
            }
            s2 * h * (b2 - a2) / two * phip(t)
        };
        let mut s1 = (gtau(t0, a1 ,b1) + gtau(tn, a1, b1)) / two;
        for i in 1..2*n {
            let t: F = (F::from(i).unwrap() - nf) * h;
            s1 = s1 + gtau(t, a1, b1);
        }
        s1 * h * (b1 - a1) / two
    }
}

fn _trapezoid<F:Float + FloatConst>(rhs: fn(F, F) -> F, x: F, a: F, b: F, n: usize) -> F {
    let zero: F = Zero::zero();
    let one: F = One::one();
    let two: F = one + one;
    if a == b {
        zero
    } else {
        let (a, b, reverse) = if a > b {
            (b, a, true)
        } else {
            (a, b, false)
        };
        let dh = (b - a) / F::from(n).unwrap();
        let s1 = (rhs(x, a) + rhs(x, b)) / two;
        let s2 = (1..n).map(|i| rhs(x, a + F::from(i).unwrap() * dh)).fold(zero, |m, j| m + j);
        if reverse {
            zero - (s1 + s2) * dh
        } else {
            (s1 + s2) * dh
        }
    }
}