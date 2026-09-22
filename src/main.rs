mod complex_numbers;
mod polynomial;
mod falcon_config_scripts;

use polynomial::Polynomial;
use complex_numbers::Complex;
use crate::polynomial::PolynomialFFT;

const n: u16 = 512;
const q: u16 = 12289;

fn main() {
    let pol : Polynomial = Polynomial::new(vec![5.0, 3.0, 1.0, 8.0]);
    let res = pol.fft();
    println!("{}", res);
    
    let res_inv = res.inv_fft();
    println!("{}", res_inv);
}
