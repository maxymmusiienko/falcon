use std::fmt;
use crate::math::complex_numbers::Complex;
use crate::math::fft::falcon_fft;
pub(crate) use crate::math::polynomial_fft::PolynomialFFT;

//todo split polynomials in different files
pub(crate) struct Polynomial {
    coefficients: Vec<f64>,
}

impl Polynomial {
    //coef are in increasing order w.r.t power of x:
    //exmpl: 5 + 3x + 4x^2
    //index 0 of coefficients is coef for x^0
    pub(crate) fn new(coefficients: Vec<f64>) -> Polynomial {
        Polynomial { coefficients }
    }

    pub(crate) fn add(&mut self, polynomial: &Polynomial) {
        let self_len = self.coefficients.len();
        let add_len = polynomial.coefficients.len();

        if self_len >= add_len {
            for i in 0..add_len {
                self.coefficients[i] += polynomial.coefficients[i];
            }
        } else {
            for i in 0..self_len {
                self.coefficients[i] += polynomial.coefficients[i];
            }
            let poly_diff = add_len - self_len;
            for i in 0..poly_diff {
                self.coefficients.push(polynomial.coefficients[self_len + i]);
            }
        }
    }

    pub(crate) fn f(&self, arg: f64) -> f64 {
        let mut sum = 0.0;
        for &coef in self.coefficients.iter().rev() {
            sum = sum * arg + coef;
        }
        sum
    }

    pub(crate) fn f_complex(&self, arg: &Complex) -> Complex {
        let mut sum = Complex::new(0.0, 0.0);
        for &coef in self.coefficients.iter().rev() {
            let complex_coef = Complex::new(coef, 0.0);
            sum = sum * (*arg) + complex_coef;
        }
        sum
    }

    pub(crate) fn dft(&self, roots: &Vec<Complex>) -> PolynomialFFT {
        let mut values = Vec::new();
        for root in roots {
            values.push(self.f_complex(root));
        }
        PolynomialFFT::new(values)
    }

    fn to_complex(&self) -> Vec<Complex> {
        let mut values = Vec::new();
        for &coef in self.coefficients.iter() {
            values.push(Complex::new(coef, 0.0));
        }
        values
    }

    pub(crate) fn fft(&self) -> PolynomialFFT {
        let complex_coefs = self.to_complex();
        let fft = falcon_fft(complex_coefs);
        PolynomialFFT::new(fft)
    }
}

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let length = self.coefficients.len();
        let mut builder = String::new();
        for i in (1..length).rev() {
            builder.push_str(&format!("{:.2}x^{} + ", self.coefficients[i], i));
        }
        builder.push_str(&format!("{:.2}", self.coefficients[0]));
        write!(f, "{}", builder)
    }
}
