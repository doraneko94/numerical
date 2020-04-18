use eom::traits::*;
use ndarray::*;

#[derive(Clone, Copy, Debug)]
pub struct Parabolic {
    pub range: (f64, f64),
    pub c: f64,
    pub dt: f64,
    pub x_part: usize,
    pub dx: f64,
}

impl Default for Parabolic {
    fn default() -> Self {
        Parabolic {
            range: (0.0, 1.0),
            c: 1.0,
            dt: 0.01,
            x_part: 6,
            dx: 1.0 / 6.0,
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
        Parabolic { range, c, dt, x_part, dx }
    }
}

impl ModelSpec for Parabolic {
    type Scalar = f64;
    type Dim = Ix1;
    fn model_size(&self) -> usize {
        self.x_part + 1
    }
}

impl Explicit for Parabolic {
    fn rhs<'a, S>(&mut self, v: &'a mut ArrayBase<S, Ix1>) -> &'a mut ArrayBase<S, Ix1>
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