use crate::complex_numbers::Complex;

 pub(crate) fn find_complex_roots(power: usize) -> Vec<Complex>{
     //find complex roots for mod polynomial of the type x^n + 1, where n = 2k
     //todo add error handling
    let mut complex_roots = Vec::new();
    let pi_f64: f64 = std::f64::consts::PI;
    for i in 0..power {
        let root = Complex::new(0.0,
                                ((2.0 * i as f64 + 1.0) * pi_f64) / power as f64);
        complex_roots.push(root.exp());
    }
    complex_roots
}
