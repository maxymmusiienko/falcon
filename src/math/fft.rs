use crate::math::complex_numbers::Complex;

pub(crate) fn falcon_fft(mut pol: Vec<Complex>) -> Vec<Complex> {
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

pub(crate) fn falcon_inv_fft(pol: Vec<Complex>) -> Vec<Complex> {
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
