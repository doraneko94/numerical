use eom::traits::*;
use ndarray::*;
use ndarray_linalg::*;
use num_traits::One;

/// Core implementation for implicit schemes
pub trait Implicit: ModelSpec {
    /// calculate simultaneous equations of implicit from current state
    fn simeq<'a, S>(&mut self, x: &'a mut ArrayBase<S, Self::Dim>) -> &'a mut ArrayBase<S, Self::Dim>
    where
        S: DataMut<Elem = Self::Scalar>;
}

#[derive(Debug, Clone)]
pub struct DiffIm<F: Implicit> {
    f: F,
    x: Array<F::Scalar, F::Dim>,
}

impl<A: Scalar, F: Implicit<Scalar = A>> TimeStep for DiffIm<F> {
    type Time = A::Real;

    fn get_dt(&self) -> Self::Time {
        One::one()
    }

    fn set_dt(&mut self, _dt: Self::Time) {
        panic!("dt of DiffIm solver is fixed to 1!");
    }
}

impl<F: Implicit> Scheme for DiffIm<F> {
    type Core = F;
    fn new(f: F, dt: Self::Time) -> Self {
        if dt != One::one() {
            panic!("dt of DiffIm solver should be fixed to 1!");
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

impl<F: Implicit> DiffIm<F> {
    pub fn from(f: F) -> Self {
        Self::new(f, One::one())
    }
}

impl<F: Implicit> ModelSpec for DiffIm<F> {
    type Scalar = F::Scalar;
    type Dim = F::Dim;
    fn model_size(&self) -> <Self::Dim as Dimension>::Pattern {
        self.f.model_size()
    }
}

impl<F:Implicit> TimeEvolution for DiffIm<F> {
    fn iterate<'a, S>(&mut self, x: &'a mut ArrayBase<S, F::Dim>) -> &'a mut ArrayBase<S, Self::Dim>
    where
        S: DataMut<Elem = Self::Scalar>,
    {
        self.x.zip_mut_with(x, |buf, x| *buf = *x);
        let fx = self.f.simeq(x);

        fx
    }
}