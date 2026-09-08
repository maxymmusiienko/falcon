use std::fmt;

pub(crate) struct Complex {
    real: f64,
    imag: f64,
}

impl Complex {
    pub(crate) fn new(real: f64, imag: f64) -> Complex {
        Complex { real, imag }
    }

    pub(crate) fn add(&self, other: &Complex) -> Complex {
        Complex::new(self.real + other.real, self.imag + other.imag)
    }
}

impl fmt::Display for Complex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} + {}i", self.real, self.imag)
    }
}
