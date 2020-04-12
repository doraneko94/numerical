use num_traits::float::Float;

pub fn legendre_core<F: Float>(n: usize) -> Vec<F> {
    let n2 = n / 2;
    let npow = 2usize.pow(n as u32);
    let one: F = One::one();
    let mut core: Vec<F> = vec![one / F::from(npow).unwrap(); n2+1];
    let mut p1: isize = 1;
    let mut p2: usize = (1..2*(n-n2)+1).map(|i| i).fold(1, |m, i| m * i);
    let mut p3: usize = 1;
    let mut p4: usize = (1..n-n2+1).map(|i| i).fold(1, |m, i| m * i);
    let mut p5: usize = (1..n-2*n2+1).map(|i| i).fold(1, |m, i| m * i);
    core[n2] = core[n2] * F::from(p2).unwrap() / (F::from(p4).unwrap() * F::from(p5).unwrap());
    for r in 1..n2+1 {
        p1 *= -1;
        p2 *= 2 * (n - n2 + r) * (2 * (n - n2 + r) - 1);
        p3 *= r;
        p4 *= n - n2 + r;
        p5 *= (n - 2 * (n2 - r)) * (n - 2 * (n2 - r) - 1);
        core[r] = core[r] * F::from(p1).unwrap() / F::from(p3).unwrap();
        core[n2 - r] = core[n2 - r] * F::from(p2).unwrap() / (F::from(p4).unwrap() * F::from(p5).unwrap());
    }
    core
}

pub fn pnx<F: Float>(n: usize, x: F, core: &Vec<F>) -> F {
    let n2 = n / 2;
    let one: F = One::one();
    let mut xn = (0..n-n2).fold(one, |m, _| m * x);
    let mut pnx = core[n2] * xn;
    for r in 1..n2+1 {
        xn = xn * x * x;
        pnx = pnx + core[n2-r] * xn;
    }
    pnx
}