mod complex_numbers;
mod polynomial;

use polynomial::Polynomial;
use complex_numbers::Complex;

const n: u16 = 512;
const q: u16 = 12289;

fn main() {
    let mut polynomial = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
    println!("first poly is f1 = {}", polynomial);
    let polynomial2 = Polynomial::new(vec![1.0, 2.0, 3.0]);
    let polynomial3 = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
    let polynomial4 = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    println!("first poly is f2 = {}", polynomial2);
    println!("first poly is f3 = {}", polynomial3);
    println!("first poly is f4 = {}", polynomial4);
    polynomial.add(&polynomial4);
    println!("f1 + f4 = {}", polynomial);
    let val = polynomial.f(2.0);
    println!("value of {} in arg = 2.0 = {}", polynomial, val);
}
