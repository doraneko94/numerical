use num_traits::float::{Float, FloatConst};
use num_traits::{Zero, One};
use crate::polynomial::{pnx, zeta_vec};
use std::collections::VecDeque;

pub struct Integral<F: Float + FloatConst> {
    pub rhs: fn(F) -> F,
}

impl<F:Float + FloatConst> Integral<F> {
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

    pub fn gauss_legendre(&self, a: F, b: F, n: usize) -> F {
        let alpha: Vec<F> = solve_pnx(n);
        let w: Vec<F> = calc_weights(&alpha);
        let zero: F = Zero::zero();
        let one: F = One::one();
        let two: F = one + one;
        let bma2 = (b - a) / two;
        let bpa2 = (b + a) / two;
        let tau = |t: F| bma2 * t + bpa2;
        let s = w.iter().zip(alpha.iter()).map(|(&wi, &ai)| {
            let xi = tau(ai);
            wi * (self.rhs)(xi)
        }).fold(zero, |m, e| m + e);
        bma2 * s
    }

    pub fn chebyshev(&self, a: F, b: F, n: usize) -> F {
        let n2 = n / 2;
        let n1f: F = F::from(2*(n+1)).unwrap();
        let theta: Vec<F> = (0..n+1).map(|k: usize| F::from(2*k+1).unwrap() / n1f * F::PI()).collect();
        let zeta = zeta_vec(n);
        let zero: F = Zero::zero();
        let one: F = One::one();
        let two: F = one + one;
        let bma2 = (b - a) / two;
        let bpa2 = (b + a) / two;
        let tau = |t: F| bma2 * t + bpa2;
        let f: Vec<F> = zeta.iter().map(|&z| (self.rhs)(tau(z))).collect();
        let denom: Vec<F> = (1..n2+1).map(|j| F::from(2*j*2*j-1).unwrap()).collect();
        let n1_hnk = |k: usize| {
            let s: F = (1..n2+1).map(|j| (F::from(2 * j).unwrap() * theta[k]).cos() / denom[j-1]).fold(zero, |m, i| m + i);
            one - two * s
        };
        let s = (0..n+1).map(|k| n1_hnk(k) * f[k]).fold(zero, |m, i| m + i);
        F::from(b - a).unwrap() / F::from(n + 1).unwrap() * s
    }
}

pub fn solve_pnx<F: Float>(n: usize) -> Vec<F> {
    let zero: F = Zero::zero();
    let one: F = One::one();
    let two: F = one + one;
    let dh: F = F::from(1e-7).unwrap();
    let mut ans = Vec::with_capacity(n);
    if n % 2 == 1 {
        ans.push(zero);
    }
    let f = |x: F| pnx(n, x);
    let mut que: VecDeque<(F, F)> = VecDeque::new();
    que.push_back((zero+dh, one-dh));
    while ans.len() < n && que.len() > 0 {
        let (r_inf, r_sup) = que.pop_front().unwrap();
        let mut inf = r_inf;
        let mut sup = r_sup;
        if inf > sup {
            continue;
        }
        let f_inf = f(inf);
        let f_sup = f(sup);
        if f_inf * f_sup > zero {
            let mid = (sup + inf) / two;
            que.push_back((inf, mid));
            que.push_back((mid, sup));
            continue;
        }
        let rev = if f_inf < zero {
            false
        } else {
            true
        };
        while sup - inf >= dh {
            let mid = (sup + inf) / two;
            let f_mid = f(mid);
            if f_mid < zero {
                if rev { sup = mid; }
                else { inf = mid; }
            } else {
                if rev { inf = mid; }
                else { sup = mid; }
            }
        }
        ans.push(inf);
        ans.push(-inf);
        que.push_back((r_inf, inf-dh));
        que.push_back((sup+dh, r_sup));
    }
    ans.sort_by(|a, b| a.partial_cmp(&b).unwrap());
    ans
}

fn rec_sum<F: Float>(i: usize, j: usize, nar:usize, dep: usize, now: usize, ar: &Vec<F>, p_sum: &Vec<F>) -> F {
    if dep == now + 1 {
        if dep == 1 {
            p_sum[0]
        } else {
            ar[j] * p_sum[j+2-dep]
        }
    } else {
        let zero: F = Zero::zero();
        let start = if now == 0 { j } else { j+1 };
        let s = (start..nar).map(|k| {
            if nar + now >= dep + k {
                rec_sum(i, k, nar, dep, now+1, ar, p_sum)
            } else {
                zero
            }
        }).fold(zero, |m, l| m + l);
        if now == 0 { s } else { s * ar[j] }
    }
}

pub fn calc_weights<F: Float>(a: &Vec<F>) -> Vec<F> {
    let zero: F = Zero::zero();
    let one: F = One::one();
    let two: F = one + one;
    let n = a.len();
    let nw = (n + 1) / 2;
    let mut w = vec![zero; n];
    for i in 0..nw {
        let mut ar = a.clone();
        let ai = ar.remove(i);
        let nar = n - 1;
        let denom = ar.iter().map(|&aj| ai - aj).fold(one, |m, ak| m * ak);
        let numer = (0..nw).map(|k| {
            let dep = nar - 2 * k;
            if dep == 0 {
                one
            } else {
                let mut p_sum = vec![zero; nar+1-dep];
                let mut s = zero;
                for l in 0..nar+1-dep {
                    s = s + ar[nar - 1 - l];
                    p_sum[nar + 1 - dep - 1 - l] = s;
                }
                rec_sum(i, 0, nar, dep, 0, &ar, &p_sum) / F::from(2 * k + 1).unwrap()
            }
        }).fold(zero, |m, s| m + s);
        let weight = (numer / denom).abs() * two;
        w[i] = weight;
        w[n-i-1] = weight;
    }
    w
}