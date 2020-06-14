use eom::traits::*;
use ndarray::*;
use ndarray_linalg::*;
//use crate::implicit::*;
use crate::difference::*;
use ndarray_linalg::solve::LUFactorized;

#[derive(Clone)]
pub struct Parabolic {
    pub range: (f64, f64),
    pub c: f64,
    pub dt: f64,
    pub x_part: usize,
    pub dx: f64,
    pub lu: LUFactorized<ndarray::OwnedRepr<f64>>,
}

fn calc_lu(range: (f64, f64), c: f64, dt: f64, x_part: usize) -> LUFactorized<ndarray::OwnedRepr<f64>> {
    let n = x_part - 1;
    let (a, b) = range;
    let dx = (b - a) / x_part as f64;
    let alpha = dt * c / (dx * dx);
    println!("{}", alpha);
    let v = (0..n*n).map(|i| {
        let j = i / n;
        let k = i % n;
        if k + 1 == j || k == j + 1 { -alpha }
        else if k == j { 1.0 + 2.0 * alpha }
        else { 0.0 }
    }).collect::<Vec<f64>>();
    let a = Array::from(v).into_shape((n, n)).unwrap();
    a.factorize_into().unwrap()
}

impl Default for Parabolic {
    fn default() -> Self {
        Self {
            range: (0.0, 1.0),
            c: 1.0,
            dt: 0.01,
            x_part: 6,
            dx: 1.0 / 6.0,
            lu: calc_lu((0.0, 1.0), 1.0, 0.01, 6)
        }
    }
}

impl Parabolic {
    pub fn new(range: (f64, f64), c: f64, dt: f64, x_part: usize) -> Self {
        let (a, b) = range;
        if a >= b {
            panic!("invalid range!");
        }
        let l = b - a;
        let dx = l / x_part as f64;
        let lu = calc_lu(range, c, dt, x_part);
       Self { range, c, dt, x_part, dx, lu }
    }
}

impl ModelSpec for Parabolic {
    type Scalar = f64;
    type Dim = Ix1;
    fn model_size(&self) -> usize {
        self.x_part + 1
    }
}

impl Difference for Parabolic {
    fn recr<'a, S>(&mut self, v: &'a mut ArrayBase<S, Ix1>) -> &'a mut ArrayBase<S, Ix1>
    where
        S: DataMut<Elem = Self::Scalar>
    {
        let v_pre = v.slice(s![1..self.x_part]).to_owned();
        let v_post = self.lu.solve_into(v_pre).unwrap();

        for i in 1..self.x_part {
            v[i] = v_post[i-1];
        }

        v
    }
}

#[derive(Clone)]
pub struct ParabolicEx {
    pub range: (f64, f64),
    pub c: f64,
    pub dt: f64,
    pub x_part: usize,
    pub dx: f64,
}

impl Default for ParabolicEx {
    fn default() -> Self {
        Self {
            range: (0.0, 1.0),
            c: 1.0,
            dt: 0.01,
            x_part: 6,
            dx: 1.0 / 6.0,
        }
    }
}

impl ParabolicEx {
    pub fn new(range: (f64, f64), c: f64, dt: f64, x_part: usize) -> Self {
        let (a, b) = range;
        if a >= b {
            panic!("invalid range!");
        }
        let l = b - a;
        let dx = l / x_part as f64;
        Self { range, c, dt, x_part, dx }
    }
}

impl ModelSpec for ParabolicEx {
    type Scalar = f64;
    type Dim = Ix1;
    fn model_size(&self) -> usize {
        self.x_part + 1
    }
}

impl Difference for ParabolicEx {
    fn recr<'a, S>(&mut self, v: &'a mut ArrayBase<S, Ix1>) -> &'a mut ArrayBase<S, Ix1>
    where
        S: DataMut<Elem = f64>,
    {
        let dx2 = (self.x_part * self.x_part) as f64;
        let xn = self.x_part;
        let v_pre = v.to_owned();
        v[0] = v_pre[0];
        for i in 1..self.x_part {
            v[i] = self.c * dx2 * (v_pre[i+1] - 2.0 * v_pre[i] + v_pre[i-1]) * self.dt + v_pre[i];
        }
        v[xn] = v_pre[xn];
        v
    }
}