use std::fmt;
use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Complex {
    pub real: f64,
    pub imag: f64,
}

impl Complex {
    pub(crate) fn new(real: f64, imag: f64) -> Complex {
        Complex { real, imag }
    }

    pub(crate) fn to_polar(&self) -> ComplexPolar {
        let norm = (self.real * self.real + self.imag * self.imag).sqrt();
        let phi = self.imag.atan2(self.real);
        ComplexPolar::new(norm, phi)
    }

    pub(crate) fn pow(&self, power: f64) -> Complex {
        let polar = self.to_polar();
        let res = polar.pow(power);
        res.to_algebraic()
    }

    pub(crate) fn exp(&self) -> Complex {
        Complex::new(self.real.exp() * self.imag.cos(),
                     self.real.exp() * self.imag.sin())
    }
}

impl Add for Complex {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Complex::new(self.real + other.real, self.imag + other.imag)
    }
}

impl Sub for Complex {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Complex::new(self.real - other.real, self.imag - other.imag)
    }
}

impl Mul for Complex {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Complex::new(
            self.real * other.real - self.imag * other.imag,
            self.imag * other.real + self.real * other.imag
        )
    }
}

impl fmt::Display for Complex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.imag > 0.0 {
            write!(f, "{:.2} + {:.2}i", self.real, self.imag)
        } else if self.imag < 0.0 {
            write!(f, "{:.2} - {:.2}i", self.real, -self.imag)
        } else {
            write!(f, "{:.2}", self.real)
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct ComplexPolar {
    norm: f64,
    phi: f64,
}

impl ComplexPolar {
    fn new(norm: f64, phi: f64) -> ComplexPolar {
        ComplexPolar { norm, phi }
    }

    fn to_algebraic(&self) -> Complex {
        Complex::new(self.norm * self.phi.cos(), self.norm * self.phi.sin())
    }

    fn pow(&self, power: f64) -> ComplexPolar {
        ComplexPolar::new(self.norm.powf(power), self.phi * power)
    }
}

//todo add unit-testing for this struct
