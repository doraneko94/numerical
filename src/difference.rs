use eom::traits::*;
use ndarray::*;
use ndarray_linalg::*;
use num_traits::One;

/// Core implementation for difference schemes
pub trait Difference: ModelSpec {
    /// calculate recurrence formula from current state
    fn recr<'a, S>(&mut self, x: &'a mut ArrayBase<S, Self::Dim>) -> &'a mut ArrayBase<S, Self::Dim>
    where
        S: DataMut<Elem = Self::Scalar>;
}

#[derive(Debug, Clone)]
pub struct Diff<F: Difference> {
    f: F,
    x: Array<F::Scalar, F::Dim>,
}

impl<A: Scalar, F: Difference<Scalar = A>> TimeStep for Diff<F> {
    type Time = A::Real;

    fn get_dt(&self) -> Self::Time {
        One::one()
    }

    fn set_dt(&mut self, _dt: Self::Time) {
        panic!("dt of Diff solver is fixed to 1!");
    }
}

impl<F: Difference> Scheme for Diff<F> {
    type Core = F;
    fn new(f: F, dt: Self::Time) -> Self {
        if dt != One::one() {
            panic!("dt of Diff solver should be fixed to 1!");
        }
        let x = Array::zeros(f.model_size());
        Self { f, x }
    }
    fn core(&self) -> &Self::Core {
        &self.f
    }
    fn core_mut(&mut self) -> &mut Self::Core {
        &mut self.f
    }
}

impl<F: Difference> Diff<F> {
    pub fn from(f: F) -> Self {
        Self::new(f, One::one())
    }
}

impl<F: Difference> ModelSpec for Diff<F> {
    type Scalar = F::Scalar;
    type Dim = F::Dim;
    fn model_size(&self) -> <Self::Dim as Dimension>::Pattern {
        self.f.model_size()
    }
}

impl<F:Difference> TimeEvolution for Diff<F> {
    fn iterate<'a, S>(&mut self, x: &'a mut ArrayBase<S, F::Dim>) -> &'a mut ArrayBase<S, Self::Dim>
    where
        S: DataMut<Elem = Self::Scalar>,
    {
        self.x.zip_mut_with(x, |buf, x| *buf = *x);
        let fx = self.f.recr(x);

        fx
    }
}