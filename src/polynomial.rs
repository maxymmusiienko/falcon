use std::fmt;
use crate::complex_numbers::Complex;

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

pub struct PolynomialFFT {
    //values of polynomial in complex roots in rev clockwise order?
    //todo ensure the correctness of the representation
    values: Vec<Complex>,
}

impl PolynomialFFT {
    pub(crate) fn new(values: Vec<Complex>) -> PolynomialFFT {
        PolynomialFFT { values }
    }

    pub(crate) fn add(&self, polynomial: &PolynomialFFT) -> PolynomialFFT {
        if self.values.len() != polynomial.values.len() {
            panic!("PolynomialFFT must have the same length");
        }
        let mut res = Vec::new();
        for i in 0..polynomial.values.len() {
            res.push(self.values[i] + polynomial.values[i]);
        }
        PolynomialFFT::new(res)
    }

    pub(crate) fn multiply(&self, polynomial: &PolynomialFFT) -> PolynomialFFT {
        if self.values.len() != polynomial.values.len() {
            panic!("PolynomialFFT must have the same length");
        }
        let mut res = Vec::new();
        for i in 0..polynomial.values.len() {
            res.push(self.values[i] * polynomial.values[i]);
        }
        PolynomialFFT::new(res)
    }

    //pub(crate) fn inv_fft(&self) -> Polynomial {
        
    //}
}

impl fmt::Display for PolynomialFFT {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut res = String::from("(");
        for i in 0..self.values.len() - 1 {
            res.push_str(self.values[i].to_string().as_str());
            res.push_str(", " );
        }
        res.push_str(self.values.last().unwrap().to_string().as_str());
        res.push_str(")");
        write!(f, "{}", res)
    }
}
