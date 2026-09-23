use std::fmt;
use crate::complex_numbers::Complex;
use crate::falcon_config_scripts::{find_complex_positive_image_roots, find_complex_roots};

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
        let fft = PolynomialFFT::falcon_fft(complex_coefs);
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

pub struct PolynomialFFT {
    //values of polynomial in complex roots in counter clockwise order?
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

    fn get_float_coefs(vec: &Vec<Complex>) -> Vec<f64> {
        let mut res = Vec::new();
        for i in 0..vec.len() {
            res.push(vec[i].real);
        }
        res
    }

    pub(crate) fn inv_fft(&self) -> Polynomial {
        let complex_coefs = self.values.clone();
        let n = complex_coefs.len() as f64;

        let inv_fft_result = Self::falcon_inv_fft(complex_coefs);

        let float_coefs = PolynomialFFT::get_float_coefs(&inv_fft_result);
        Polynomial::new(float_coefs)
    }

    fn falcon_fft(mut pol: Vec<Complex>) -> Vec<Complex> {
        let n = pol.len();

        let psi_angle = std::f64::consts::PI / n as f64;
        let psi = Complex::new(0.0, psi_angle).exp();
        let mut current_psi = Complex::new(1.0, 0.0);

        for j in 0..n {
            pol[j] = pol[j] * current_psi;
            current_psi = current_psi * psi;
        }

        fft_routine(pol, false)
    }

    fn falcon_inv_fft(pol: Vec<Complex>) -> Vec<Complex> {
        let n = pol.len();

        let mut res = fft_routine(pol, true);

        let psi_angle = -std::f64::consts::PI / n as f64;
        let psi = Complex::new(0.0, psi_angle).exp();
        let mut current_psi = Complex::new(1.0, 0.0);

        let inv_n = 1.0 / n as f64;

        for j in 0..n {
            res[j] = res[j] * current_psi * inv_n;
            current_psi = current_psi * psi;
        }

        res
    }

    pub(crate) fn splitfft(&self) -> (PolynomialFFT, PolynomialFFT) {
        let n = self.values.len();
        let mut pol_even = PolynomialFFT::new(Vec::with_capacity(n / 2));
        let mut pol_odd = PolynomialFFT::new(Vec::with_capacity(n / 2));
        let roots = find_complex_positive_image_roots(n);
        for i in 0..n / 2 {
            let even_arg = (self.values[i] + self.values[i + n / 2]) * 0.5;
            pol_even.values.push(even_arg);
            let odd_arg = (self.values[i] - self.values[i + n / 2]) * 0.5 * roots[i].conj();
            pol_odd.values.push(odd_arg);
        }
        (pol_even, pol_odd)
    }

    pub(crate) fn mergefft(even_polynomial: &PolynomialFFT, odd_polynomial: &PolynomialFFT) -> PolynomialFFT {
        let n = even_polynomial.values.len() * 2;
        //todo take the roots from cache
        let roots = find_complex_positive_image_roots(n);
        let mut merged_values = vec![Complex::new(0.0, 0.0); n];

        for i in 0..n / 2 {
            let zeta = roots[i];
            let term = zeta * odd_polynomial.values[i];
            merged_values[i] = even_polynomial.values[i] + term;
            merged_values[i + n / 2] = even_polynomial.values[i] - term;
        }

        PolynomialFFT::new(merged_values)
    }
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

fn fft_routine(pol: Vec<Complex>, inv: bool) -> Vec<Complex> {
    let n = pol.len();
    if n == 1 {
        return pol;
    }

    let angle = if inv {
        -2.0 * std::f64::consts::PI / n as f64
    } else {
        2.0 * std::f64::consts::PI / n as f64
    };

    let complex_arg = Complex::new(0.0, angle);
    let w = Complex::exp(&complex_arg);

    let mut pol_even: Vec<Complex> = Vec::with_capacity(n / 2);
    let mut pol_odd: Vec<Complex> = Vec::with_capacity(n / 2);

    for i in 0..n / 2 {
        pol_even.push(pol[2 * i]);
        pol_odd.push(pol[2 * i + 1]);
    }

    let y_even = fft_routine(pol_even, inv);
    let y_odd = fft_routine(pol_odd, inv);

    let mut y: Vec<Complex> = vec![Complex::new(0.0, 0.0); n];
    let mut w_j = Complex::new(1.0, 0.0);

    for j in 0..n / 2 {
        let term = w_j * y_odd[j];

        y[j] = y_even[j] + term;
        y[j + n / 2] = y_even[j] - term;

        w_j = w_j * w;
    }

    y
}
