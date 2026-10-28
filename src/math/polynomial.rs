use std::fmt;
use crate::math::complex_numbers::Complex;
use crate::math::fft::FFTUtil;
pub(crate) use crate::math::polynomial_fft::PolynomialFFT;

//todo split polynomials in different files
pub(crate) struct Polynomial {
    pub(crate) coefficients: Vec<f64>,
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
        let fft = FFTUtil::falcon_fft(complex_coefs);
        PolynomialFFT::new(fft)
    }

    fn split(&self) -> (Polynomial, Polynomial) {
        let n = self.coefficients.len();
        let mut pol_even = Polynomial::new(Vec::with_capacity(n / 2));
        let mut pol_odd = Polynomial::new(Vec::with_capacity(n / 2));

        for i in 0..n / 2 {
            pol_even.coefficients.push(self.coefficients[2 * i]);
            pol_odd.coefficients.push(self.coefficients[2 * i + 1]);
        }

        (pol_even, pol_odd)
    }

    fn modular_reduction(&self, n: usize) -> Polynomial {
        let mut reduced_coeffs = vec![0.0; n];

        let lower_bound = std::cmp::min(n, self.coefficients.len());
        for i in 0..lower_bound {
            reduced_coeffs[i] = self.coefficients[i];
        }

        for i in n..self.coefficients.len() {
            reduced_coeffs[i - n] -= self.coefficients[i];
        }

        Polynomial {coefficients: reduced_coeffs}
    }

    fn mul_by_x(&self) -> Polynomial {
        let mut new_coeffs = Vec::with_capacity(self.coefficients.len() + 1);
        new_coeffs.push(0.0);
        new_coeffs.extend_from_slice(&self.coefficients);

        Polynomial {coefficients: new_coeffs}
    }

    fn mul_by_scalar(&self, scalar: f64) -> Polynomial {
        let new_coeffs = self.coefficients.iter().map(|&c| c * scalar).collect();
        Polynomial {coefficients: new_coeffs}
    }

    pub(crate) fn field_norm(&self) -> Polynomial {
        let n = self.coefficients.len();
        let (f0, f1) : (Polynomial, Polynomial) = self.split();
        let f0_fft = f0.fft();
        let f1_fft = f1.fft();
        let f0_sq_fft = f0_fft.multiply_fft(&f0_fft);
        let f1_sq_fft = f1_fft.multiply_fft(&f1_fft);
        let mut f0_sq = f0_sq_fft.inv_fft();
        let f1_sq = f1_sq_fft.inv_fft();

        let f1_sq_x = f1_sq.mul_by_x();
        let opposite = f1_sq_x.mul_by_scalar(-1.0);
        f0_sq.add(&opposite);
        f0_sq.modular_reduction(n / 2)
    }

    pub(crate) fn get_polynomial_in_square(&self) -> Polynomial {
        let n = self.coefficients.len();
        let mut res : Vec<f64> = Vec::with_capacity(n * 2);
        for i in 0..n {
            res.push(0.0);
            res.push(self.coefficients[i]);
        }
        Polynomial::new(res)
    }

    pub(crate) fn get_polynomial_adjoint(&self) -> Polynomial {
        let n = self.coefficients.len();
        let mut res = self.coefficients.clone();
        for i in 1..n {
            res[i] *= -1.0;
        }
        Polynomial::new(res)
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
