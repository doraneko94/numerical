use num_traits::float::Float;
use num_traits::{Zero, One};
use ndarray::*;
use ndarray_linalg::*;
use ndarray_linalg::lapack::Lapack;

pub struct Spline3d<F: Float + Lapack> {
    pub x: Vec<F>,
    pub y: Vec<F>,
    pub abc: Vec<(F, F, F)>,
}

impl<F: Float + Lapack> Spline3d<F> {
    pub fn new(x: &Vec<F>, y: &Vec<F>, c0: F, cn: F) -> Self {
        let zero: F = Zero::zero();
        let one: F = One::one();
        let two: F = one + one;
        let three: F = two + one;
        let n = x.len();
        if n != y.len() {
            panic!("x and y have different number of elements!");
        }
        let mut xy: Vec<(F, F)> = x.iter().zip(y.iter()).map(|(&xi, &yi)| (xi, yi)).collect();
        xy.sort_by(|a, b| (a.0).partial_cmp(&b.0).unwrap());
        let h: Vec<F> = (0..n-1).map(|i| xy[i+1].0 - xy[i].0).collect();
        let u: Vec<F> = (0..n-1).map(|i| (xy[i+1].1 - xy[i].1) / h[i]).collect();
        let mut a: Array2<F> = Array2::zeros((n-2, n-2));
        let mut b: Vec<F> = vec![zero; n-2];
        for i in 0..n-2 {
            if i == 0 {
                b[i] = b[i] - h[i+1] * c0;
                a[[i, i+1]] = h[i];
            } else if i == n-3 {
                b[i] = b[i] - h[i] * cn;
                a[[i, i-1]] = h[i+1];
            } else {
                a[[i, i-1]] = h[i+1];
                a[[i, i+1]] = h[i];
            }
            a[[i, i]] = two * (h[i] + h[i+1]);
            b[i] = b[i] + three * (h[i+1] * u[i] + h[i] * u[i+1]);
        }

        let b = Array1::from(b);
        let mut c = a.solve(&b).unwrap().to_vec();
        c.insert(0, c0);
        c.push(cn);
        let b: Vec<F> = (0..n-1).map(|i| (three * u[i] - two * c[i] - c[i+1]) / h[i]).collect();
        let a: Vec<F> = (0..n-1).map(|i| (u[i] - h[i] * b[i] - c[i]) / (h[i] * h[i])).collect();
        let abc: Vec<(F, F, F)> = (0..n-1).map(|i| (a[i], b[i], c[i])).collect();

        Self { x: x.clone(), y: y.clone(), abc }
    }

    pub fn calc(&self, z: F) -> Result<F, &str> {
        let n = self.x.len();
        match self.x.binary_search_by(|xi| xi.partial_cmp(&z).unwrap()) {
            Ok(index) => Ok(self.y[index]),
            Err(index) => {
                if index == 0 || index == n {
                    Err("out of range!")
                } else {
                    let (a, b, c) = self.abc[index-1];
                    let mut m: F = One::one();
                    let mut ans = self.y[index-1];
                    m = m * (z - self.x[index-1]);
                    ans = ans + c * m;
                    m = m * (z - self.x[index-1]);
                    ans = ans + b * m;
                    m = m * (z - self.x[index-1]);
                    ans = ans + a * m;
                    Ok(ans)
                }
            }
        }
    }
}

pub struct LSM<F: Float + Lapack> {
    pub x: Vec<F>,
    pub y: Vec<F>,
    pub items: Vec<fn(F)->F>,
    pub coeff: Vec<F>,
}

impl<F: Float + Lapack> LSM<F> {
    pub fn new(x: &Vec<F>, y: &Vec<F>, items: Vec<fn(F)->F>) -> Self {
        let n = x.len();
        if n != y.len() {
            panic!("x and y have different number of elements!");
        }
        let nf = items.len();
        if nf == 0 {
            panic!("items should contain at least 1 function!");
        }
        let mut a = Array2::zeros((n, nf));
        let mut b = Array1::zeros(n);
        for i in 0..n {
            for j in 0..nf {
                a[[i, j]] = items[j](x[i]);
            }
            b[i] = y[i];
        }
        let taa = a.t().dot(&a);
        let tab = b.dot(&a);
        let coeff = taa.solve(&tab).unwrap().to_vec();

        Self { x: x.clone(), y: y.clone(), items, coeff, }
    }

    pub fn calc(&self, z: F) -> F {
        let zero: F = Zero::zero();
        self.items.iter().zip(self.coeff.iter()).map(|(&i, &c)| c * i(z)).fold(zero, |m, j| m + j)
    }
}