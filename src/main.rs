mod complex_numbers;
mod polynomial;
mod falcon_config_scripts;

use polynomial::Polynomial;
use complex_numbers::Complex;

const n: u16 = 512;
const q: u16 = 12289;

fn main() {
    let complex_roots = falcon_config_scripts::find_complex_roots(4);
    for root in &complex_roots {
        println!("{}", root);
    }
    let polynomial1 = Polynomial::new(vec![1.0, 2.0, -1.0, 0.0]);
    let polynomial2 = Polynomial::new(vec![2.0, -1.0, 1.0, 1.0]);

    let fft_polynomial1 = polynomial1.fft(&complex_roots);
    let fft_polynomial2 = polynomial2.fft(&complex_roots);

    println!("{}", fft_polynomial1);
    println!("{}", fft_polynomial2);

    let fft_sum = fft_polynomial1.add(&fft_polynomial2);
    println!("{}", fft_sum);
    let fft_mult = fft_polynomial1.multiply(&fft_polynomial2);
    println!("{}", fft_mult);
}
