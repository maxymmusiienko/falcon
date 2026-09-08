use std::fmt;

const n: u16 = 512;
const q: u16 = 12289;

struct Polynomial {
    coefficients: Vec<f64>,
}

impl Polynomial {
    fn new(coefficients: Vec<f64>) -> Polynomial {
        Polynomial { coefficients }
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

fn main() {
    let polynomial = Polynomial::new(vec![1.0, 2.0, 3.0, 4.0]);
    println!("{}", polynomial);
}
