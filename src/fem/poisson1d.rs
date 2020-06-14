use ndarray::*;
use ndarray_linalg::*;
use ndarray_linalg::error::Result;
use num_traits::{Float, FloatConst};
use num_traits::{Zero, One};
use cauchy::Scalar;

/// c d^2u/(dx)^2 = f(x)
pub struct Poisson1dFEM<A: Float + FloatConst + Scalar + Lapack> {
    pub c: A,
    pub f: fn(A)->A,
    pub range: (A, A),
    pub n: usize,
    pub dirichlet: (A, A),
    pub neumann: (A, A),
}

impl<A> Default for Poisson1dFEM<A>
where
    A: Float + FloatConst + Scalar + Lapack
{
    fn default() -> Self {
        let one = A::one();
        let nan = A::nan();
        Self { 
            c: one, 
            f: |_x: A| One::one(),
            range: (-one, one),
            n: 3,
            dirichlet: (one, nan),
            neumann: (nan, -one), }
    }
}

impl<A> Poisson1dFEM<A>
where
    A: Float + FloatConst + Scalar + Lapack
{
    pub fn calc(&self) -> Result<Array1<A>> {
        let mut u: Array1<A> = Array::zeros(self.n + 1);

        let n = A::from(self.n).unwrap();
        let (x0, x1) = self.range;
        let dx = (x1 - x0) / n;
        let dx_inv = n / (x1 - x0);
        let dx_2 = dx / (A::one() + A::one());

        let x: Array1<A> = Array::from((0..self.n+1)
                                        .map(|i| x0 + dx * A::from(i).unwrap())
                                        .collect::<Vec<A>>());
        let fx: Array1<A> = Array::from((0..self.n)
                                        .map(|i| -(self.f)(x[i] + dx_2) / self.c)
                                        .collect::<Vec<A>>());
                                        
        let mut du: Array1<A> = Array::from(vec![-dx_inv; self.n]);
        let mut d: Array1<A> = Array::from(vec![dx_inv + dx_inv; self.n+1]);
        d[0] = dx_inv;
        d[self.n] = dx_inv;
        let mut dl: Array1<A> = Array::from(vec![-dx_inv; self.n]);

        
        let mut f: Array1<A> = Array::from((0..self.n+1)
                                            .map(|i| {
                                                if i == 0 { dx_2 * fx[0] }
                                                else if i == self.n { dx_2 * fx[self.n-1] }
                                                else { dx_2 * (fx[i-1] + fx[i]) }
                                            })
                                            .collect::<Vec<A>>());
        
        let (a1, a2) = self.dirichlet;
        if !a1.is_nan() {
            du = du.slice(s![1..]).to_owned();
            d = d.slice(s![1..]).to_owned();
            f = f.slice(s![1..]).to_owned();
            f[0] = f[0] - dl[0] * a1;
            dl = dl.slice(s![1..]).to_owned();
            u[0] = a1;
        }
        if !a2.is_nan() {
            let n_col = d.len();
            dl = dl.slice(s![..n_col-2]).to_owned();
            d = d.slice(s![..n_col-1]).to_owned();
            f = f.slice(s![..n_col-1]).to_owned();
            f[n_col-2] = f[n_col-2] - du[n_col-2] * a2;
            du = du.slice(s![..n_col-2]).to_owned();
            u[self.n] = a2;
        }
        let (b1, b2) = self.neumann;
        if !b1.is_nan() {
            if !a1.is_nan() {
                panic!();
            }
            f[0] = f[0] + b1;
        }
        if !b2.is_nan() {
            if !a2.is_nan() {
                panic!();
            }
            let n_col = f.len();
            f[n_col-1] = f[n_col-1] + b2;
        }

        let u_ster = if a1.is_nan() && a2.is_nan() {
            u.slice_mut(s![..])
        } else if a1.is_nan() {
            u.slice_mut(s![..u.len()-1])
        } else if a2.is_nan() {
            u.slice_mut(s![1..])
        } else {
            u.slice_mut(s![1..u.len()-1])
        };

        let n_col = d.len() as i32;
        let l = MatrixLayout::C((n_col, n_col));
        let tridiag = TriDiagonal { l: l, n1: <A as Scalar>::Real::zero(), dl: dl, d: d, du: du, };
        let u_solve = tridiag.solve_tridiagonal_into(f)?;

        Zip::from(u_ster)
            .and(&u_solve)
            .apply(|a, &b| *a = b);
        
        Ok(u)
    }
}