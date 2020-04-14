use num_traits::float::{Float, FloatConst};
use num_traits::{Zero, One};

pub struct Chebyshev<F: Float + FloatConst> {
    pub rhs: fn(F) -> F,
    pub c: Vec<F>,
}

impl<F: Float + FloatConst> Chebyshev<F> {
    pub fn new(rhs: fn(F)->F, n: usize) -> Self {
        let zero: F = Zero::zero();
        let zeta = zeta_vec(n);
        let f: Vec<F> = zeta.iter().map(|&z| rhs(z)).collect();
        let mut c = Vec::with_capacity(n+1);
        c.push(F::from(1.0 / (n + 1) as f64).unwrap() * f.iter().fold(zero, |m, &i| m + i));
        let tnx: Vec<Vec<F>> = zeta.iter().map(|&z| tnx_vec(n, z)).collect();
        for i in 1..n+1 {
            c.push(F::from(2.0 / (n + 1) as f64).unwrap() * (0..n+1).map(|j| f[j] * tnx[j][i]).fold(zero, |m, k| m + k));
        }
        Self { rhs, c }
    }

    pub fn calc(&self, x: F) -> Result<F, &str> {
        let one: F = One::one();
        if x.abs() > one {
            return Err("x should be in [-1, 1]");
        }
        let zero: F = Zero::zero();
        let n = self.c.len() - 1;
        let tnx: Vec<F> = tnx_vec(n, x);
        Ok(self.c.iter().zip(tnx.iter()).map(|(&ci, &tix)| ci * tix).fold(zero, |m, i| m + i))
    }
}

pub fn tnx_vec<F: Float + FloatConst>(n: usize, x: F) -> Vec<F> {
    let one: F = One::one();
    let two: F = one + one;
    if n == 0 {
        vec![one]
    } else if n == 1 {
        vec![one, x]
    } else {
        let mut ret = Vec::with_capacity(n+1);
        ret.push(one);
        ret.push(x);
        for i in 1..n+1 {
            ret.push(two * x * ret[i] - ret[i-1]);
        }
        ret
    }
}

pub fn zeta_vec<F: Float + FloatConst>(n: usize) -> Vec<F> {
    let mut ret = Vec::with_capacity(n+1);
    let n1f = F::from(2 * (n+1)).unwrap();
    for i in 0..n+1 {
        ret.push((F::from(2 * i + 1).unwrap() / n1f * F::PI()).cos());
    }
    ret
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