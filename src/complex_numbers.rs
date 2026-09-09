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

    pub(crate) fn mult(&self, other: &Complex) -> Complex {
        Complex::new(self.real * other.real - self.imag * other.imag,
                     self.imag * other.real + self.real * other.imag)
    }

    pub(crate) fn to_polar(&self) -> ComplexPolar {
        let norm = (self.real * self.real + self.imag * self.imag).sqrt();
        let phi = self.imag.atan2(self.real);
        ComplexPolar::new(norm, phi) // Змінено порядок на (norm, phi)
    }
}

impl fmt::Display for Complex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} + {}i", self.real, self.imag)
    }
}

pub(crate) struct ComplexPolar {
    norm: f64,
    phi: f64,
}

impl ComplexPolar {
    fn new(norm: f64, phi: f64) -> ComplexPolar {
        ComplexPolar { norm, phi }
    }

    pub(crate) fn to_algebraic(&self) -> Complex {
        Complex::new(self.norm * self.phi.cos(), self.norm * self.phi.sin())
    }

    pub(crate) fn pow(&self, power: f64) -> ComplexPolar {
        ComplexPolar::new(self.norm.powf(power), self.phi * power)
    }
}
