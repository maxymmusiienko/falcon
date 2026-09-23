mod complex_numbers;
mod polynomial;
mod falcon_config_scripts;

use polynomial::Polynomial;
use complex_numbers::Complex;
use crate::falcon_config_scripts::find_complex_roots;
use crate::polynomial::PolynomialFFT;

const n: u16 = 512;
const q: u16 = 12289;

fn main() {
    let pol : Polynomial = Polynomial::new(vec![5.0, 3.0, 1.0, 8.0]);
    println!("{}", pol);
    let res = pol.fft();
    println!("{}", res);

    let (p1, p2) = res.splitfft();
    println!("{}", p1);
    println!("{}", p2);

    let p1_int = p1.inv_fft();
    let p2_int = p2.inv_fft();
    println!("{}", p1_int);
    println!("{}", p2_int);

    let mergedfft = PolynomialFFT::mergefft(&p1, &p2);
    println!("{}", mergedfft);

    let mergedint = mergedfft.inv_fft();
    println!("{}", mergedint);
}
