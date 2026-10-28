use crate::crypto::crypto_utils::CryptoUtil;
use crate::math::fft::FFTUtil;
use crate::math::polynomial::Polynomial;
use crate::utils::config::FalconDegree;

pub(crate) struct NtruSolved {
    pub(crate) F : Polynomial,
    pub(crate) G : Polynomial,
}

impl NtruSolved {
    pub(crate) fn new(F : Polynomial, G : Polynomial) -> Self {
        NtruSolved {F,G}
    }
}

pub(crate) fn ntru_solve(f : Polynomial, g : Polynomial, n : u32) -> Option<NtruSolved> {
    if n == 1 {
        let f_z = f.coefficients[0] as i64;
        let g_z = g.coefficients[0] as i64;
        let egcd_res = CryptoUtil::egcd(f_z, g_z);
        let gcd = egcd_res.0;
        if gcd != 1 {
            return None;
        }
        let q = FalconDegree::q as i64;
        let F = Polynomial::new(vec![(egcd_res.2 * q) as f64]);
        let G = Polynomial::new(vec![(egcd_res.1 * q) as f64]);
        let ntru_solved = NtruSolved::new(F, G);
        Some(ntru_solved)
    } else {
        let f_prime = f.field_norm();
        let g_prime = g.field_norm();
        let ntru_solved = ntru_solve(f_prime, g_prime, n / 2);
        if let Some(solved) = ntru_solved {
            let f_prime_sq: Polynomial = solved.F.get_polynomial_in_square();
            let g_prime_sq: Polynomial = solved.G.get_polynomial_in_square();
            let f_adj: Polynomial = f.get_polynomial_adjoint();
            let g_adj: Polynomial = g.get_polynomial_adjoint();
            let F = FFTUtil::multiply_two_polynomials(&f_prime_sq, &g_adj);
            let G = FFTUtil::multiply_two_polynomials(&g_prime_sq, &f_adj);
            //todo Reduce function
            Some(NtruSolved::new(F, G))
        } else {
            return None;
        }
    }
}
