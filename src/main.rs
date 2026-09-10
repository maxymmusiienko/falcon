mod complex_numbers;
mod polynomial;
mod falcon_config_scripts;

use polynomial::Polynomial;
use complex_numbers::Complex;

const n: u16 = 512;
const q: u16 = 12289;

fn main() {
    let complex_roots = falcon_config_scripts::find_complex_roots(2);
    for root in complex_roots {
        println!("{}", root);
    }
}
