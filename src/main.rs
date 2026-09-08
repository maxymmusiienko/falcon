use std::fmt;

const n: u16 = 512;
const q: u16 = 12289;
struct Polynomial {
    coefficients: Vec<f64>,
}

impl Polynomial {
    //coef are in increasing order w.r.t power of x:
    //exmpl: 5 + 3x + 4x^2
    //index 0 of coefficients is coef for x^0
    fn new(coefficients: Vec<f64>) -> Polynomial {
        Polynomial { coefficients }
    }

    fn add(&mut self, polynomial: &Polynomial) {
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

    fn f(&self, arg: f64) -> f64 {
        let mut sum = self.coefficients[0];
        for i in 1..self.coefficients.len() {
            sum += self.coefficients[i] * arg.powf(i as f64);
        }
        sum
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
