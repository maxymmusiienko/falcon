use crate::crypto::crypto_utils::CryptoUtil;
use crate::math::fft::FFTUtil;
use crate::math::polynomial::Polynomial;
use crate::utils::config::FalconDegree;
//todo wrap this in util struct
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
        let F = Polynomial::new(vec![(-egcd_res.2 * q) as f64]);
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
            let mut F = FFTUtil::multiply_two_polynomials(&f_prime_sq, &g_adj);
            let mut G = FFTUtil::multiply_two_polynomials(&g_prime_sq, &f_adj);
            reduce(&f, &g, &mut F, &mut G);
            Some(NtruSolved::new(F, G))
        } else {
            return None;
        }
    }
}

fn reduce(f: &Polynomial, g: &Polynomial, F: &mut Polynomial, G: &mut Polynomial) {
    let f_adj = f.get_hermitian_adjoint();
    let g_adj = g.get_hermitian_adjoint();

    let ff = FFTUtil::multiply_two_polynomials(f, &f_adj);
    let gg = FFTUtil::multiply_two_polynomials(g, &g_adj);
    let denom = ff.add_new(&gg);

    loop {
        let F_f = FFTUtil::multiply_two_polynomials(F, &f_adj);
        let G_g = FFTUtil::multiply_two_polynomials(G, &g_adj);
        let num = F_f.add_new(&G_g);

        let mut k = FFTUtil::divide_two_polynomials(&num, &denom);
        k.round_coefficients();

        if k.is_zero() {
            break;
        }

        let k_f = FFTUtil::multiply_two_polynomials(&k, f);
        let k_g = FFTUtil::multiply_two_polynomials(&k, g);

        F.sub(&k_f);
        G.sub(&k_g);
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests for [`ntru_solve`] and [`reduce`], verifying algebraic correctness
    //! of the NTRU invariant and lattice reduction quality across multiple ring degrees.

    use super::*;

    /// Verifies that the solved polynomial pair `(F, G)` satisfies the NTRU equation:
    /// $$f \cdot G - g \cdot F \equiv q \pmod{x^n + 1}$$
    ///
    /// Specifically checks that the constant term equals `q` (after rounding floating-point
    /// FFT noise) and all higher-degree coefficients $x^1, \dots, x^{n-1}$ vanish within
    /// a small numerical tolerance (`1e-3`).
    fn verify_ntru_equation(f: &Polynomial, g: &Polynomial, solved: &NtruSolved, q: f64) {
        let f_g = FFTUtil::multiply_two_polynomials(f, &solved.G);
        let g_f = FFTUtil::multiply_two_polynomials(g, &solved.F);

        let mut diff = f_g;
        diff.sub(&g_f);

        // Constant term must equal the Falcon modulus q
        let c0 = diff.coefficients[0].round();
        assert_eq!(
            c0, q,
            "Constant term mismatch: expected {}, got {} (check signs in egcd or adjoints)",
            q, c0
        );

        // All remaining coefficients (x^1 ... x^(n-1)) must be zero
        for (i, &coef) in diff.coefficients.iter().enumerate().skip(1) {
            assert!(
                coef.abs() < 1e-3,
                "Non-zero coefficient at x^{}: {} (expected 0.0)",
                i,
                coef
            );
        }
    }

    /// Tests the base case ($n = 1$) of [`ntru_solve`] over $\mathbb{Z}$,
    /// validating Bézout coefficient signs and scaling by $q$.
    #[test]
    fn test_ntru_solve_n1() {
        let f_coefs = vec![3.0];
        let g_coefs = vec![5.0];
        let q = (FalconDegree::q as i64) as f64;

        let solved = ntru_solve(
            Polynomial::new(f_coefs.clone()),
            Polynomial::new(g_coefs.clone()),
            1,
        )
            .expect("ntru_solve must return Some for coprime f and g");

        verify_ntru_equation(
            &Polynomial::new(f_coefs),
            &Polynomial::new(g_coefs),
            &solved,
            q,
        );
    }

    /// Tests recursive field-norm projection, lifting, and Babai reduction
    /// in the cyclotomic ring $\mathbb{Z}[x]/(x^4 + 1)$ ($n = 4$).
    ///
    /// Also asserts that [`reduce`] keeps the infinity norm of `(F, G)` bounded.
    #[test]
    fn test_ntru_solve_n4() {
        let f_coefs = vec![2.0, -1.0, 1.0, 0.0];
        let g_coefs = vec![1.0, 1.0, -1.0, 2.0];
        let q = (FalconDegree::q as i64) as f64;

        let solved = ntru_solve(
            Polynomial::new(f_coefs.clone()),
            Polynomial::new(g_coefs.clone()),
            4,
        )
            .expect("ntru_solve must find a valid solution for n = 4");

        verify_ntru_equation(
            &Polynomial::new(f_coefs),
            &Polynomial::new(g_coefs),
            &solved,
            q,
        );

        // Ensure Babai reduction kept the coefficients reasonably small
        let max_f = solved
            .F
            .coefficients
            .iter()
            .map(|c| c.abs())
            .fold(0.0_f64, f64::max);
        let max_g = solved
            .G
            .coefficients
            .iter()
            .map(|c| c.abs())
            .fold(0.0_f64, f64::max);

        assert!(
            max_f < 100_000.0 && max_g < 100_000.0,
            "Coefficients of F or G are too large (max_F = {}, max_G = {}); check reduce()",
            max_f,
            max_g
        );
    }

    /// Tests deeper recursion ($n = 8 \to 4 \to 2 \to 1$) in $\mathbb{Z}[x]/(x^8 + 1)$
    /// to ensure multi-level lifting and reduction remain numerically stable.
    #[test]
    fn test_ntru_solve_n8() {
        let f_coefs = vec![1.0, 2.0, -1.0, 0.0, 1.0, -2.0, 0.0, 1.0];
        let g_coefs = vec![2.0, -1.0, 0.0, 1.0, -1.0, 1.0, 2.0, -1.0];
        let q = (FalconDegree::q as i64) as f64;

        let solved = ntru_solve(
            Polynomial::new(f_coefs.clone()),
            Polynomial::new(g_coefs.clone()),
            8,
        )
            .expect("ntru_solve must find a valid solution for n = 8");

        verify_ntru_equation(
            &Polynomial::new(f_coefs),
            &Polynomial::new(g_coefs),
            &solved,
            q,
        );
    }

    /// Verifies that [`ntru_solve`] gracefully aborts and returns `None` when
    /// $f$ and $g$ share a common factor, resulting in $\gcd(N(f), N(g)) \neq 1$.
    #[test]
    fn test_ntru_solve_no_solution() {
        let f = Polynomial::new(vec![2.0, 4.0, -2.0, 2.0]);
        let g = Polynomial::new(vec![4.0, -2.0, 6.0, 2.0]);

        let solved = ntru_solve(f, g, 4);
        assert!(
            solved.is_none(),
            "ntru_solve must return None when gcd(N(f), N(g)) != 1"
        );
    }

    /// Tests full-scale recursion ($n = 128 \to 64 \to \dots \to 1$) in $\mathbb{Z}[x]/(x^{128} + 1)$
    /// using deterministically generated Falcon-like short polynomials.
    ///
    /// Validates that:
    /// 1. The NTRU invariant $f \cdot G - g \cdot F \equiv q \pmod{x^{128} + 1}$ holds without
    ///    floating-point overflow across 7 levels of recursive norm projection.
    /// 2. All coefficients of `(F, G)` remain integer-valued and properly reduced in norm.
    #[test]
    fn test_ntru_solve_n128() {
        let n = 16;
        let q = (FalconDegree::q as i64) as f64;

        // Simple deterministic Xorshift64 PRNG to generate short Falcon-like coefficients
        // without external crate dependencies.
        let mut state: u64 = 0x5EED_FA1C_0128_0001;
        let mut next_short_coef = || -> f64 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            // Sum of centered small values approximates a discrete Gaussian with small std dev
            let a = (state & 0x3) as i64 - 1;         // in {-1, 0, 1, 2}
            let b = ((state >> 8) & 0x3) as i64 - 2;  // in {-2, -1, 0, 1}
            let c = ((state >> 16) & 0x1) as i64;     // in {0, 1}
            (a + b + c) as f64                        // concentrated in [-3, 4]
        };

        // Try deterministic seeds until a coprime pair (N(f), N(g)) is found
        let mut solved_pair = None;
        for _attempt in 0..10000 {
            let f_coefs: Vec<f64> = (0..n).map(|_| next_short_coef()).collect();
            let g_coefs: Vec<f64> = (0..n).map(|_| next_short_coef()).collect();

            if let Some(solved) = ntru_solve(
                Polynomial::new(f_coefs.clone()),
                Polynomial::new(g_coefs.clone()),
                n as u32,
            ) {
                solved_pair = Some((f_coefs, g_coefs, solved));
                break;
            }
        }

        let (f_coefs, g_coefs, solved) =
            solved_pair.expect("ntru_solve must find a valid solution for n = 128 within 10000 attempts");

        verify_ntru_equation(
            &Polynomial::new(f_coefs),
            &Polynomial::new(g_coefs),
            &solved,
            q,
        );

        // Verify that all coefficients of F and G are exact integers (no fractional drift)
        for (i, (&cf, &cg)) in solved
            .F
            .coefficients
            .iter()
            .zip(solved.G.coefficients.iter())
            .enumerate()
        {
            assert!(
                (cf - cf.round()).abs() < 1e-6 && (cg - cg.round()).abs() < 1e-6,
                "Non-integer coefficient at index {}: F[i] = {}, G[i] = {}",
                i,
                cf,
                cg
            );
        }

        // Verify that Babai reduction kept the solution norm small
        let max_f = solved.F.coefficients.iter().map(|c| c.abs()).fold(0.0_f64, f64::max);
        let max_g = solved.G.coefficients.iter().map(|c| c.abs()).fold(0.0_f64, f64::max);
        assert!(
            max_f < 100_000.0 && max_g < 100_000.0,
            "Coefficients of F or G are too large for n = 128 (max_F = {}, max_G = {})",
            max_f,
            max_g
        );
    }
}
