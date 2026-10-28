use crate::math::complex_numbers::Complex;
use crate::math::fft::FFTUtil;
use crate::math::polynomial::Polynomial;
use crate::utils::falcon_config_scripts::find_complex_positive_image_roots;
use std::fmt;

pub(crate) struct PolynomialFFT {
    //values of polynomial in complex roots in counter clockwise order?
    //todo ensure the correctness of the representation
    values: Vec<Complex>,
}

impl PolynomialFFT {
    pub(crate) fn new(values: Vec<Complex>) -> PolynomialFFT {
        PolynomialFFT { values }
    }

    pub(crate) fn add_fft(&self, polynomial: &PolynomialFFT) -> PolynomialFFT {
        if self.values.len() != polynomial.values.len() {
            panic!("PolynomialFFT must have the same length");
        }
        let mut res = Vec::new();
        for i in 0..polynomial.values.len() {
            res.push(self.values[i] + polynomial.values[i]);
        }
        PolynomialFFT::new(res)
    }

    pub(crate) fn multiply_fft(&self, polynomial: &PolynomialFFT) -> PolynomialFFT {
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

        let inv_fft_result = FFTUtil::falcon_inv_fft(complex_coefs);

        let float_coefs = PolynomialFFT::get_float_coefs(&inv_fft_result);
        Polynomial::new(float_coefs)
    }

    pub(crate) fn splitfft(&self) -> (PolynomialFFT, PolynomialFFT) {
        let n = self.values.len();
        let mut pol_even = PolynomialFFT::new(Vec::with_capacity(n / 2));
        let mut pol_odd = PolynomialFFT::new(Vec::with_capacity(n / 2));
        //todo precompute it, use context
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
        //todo precompute it, use context
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
