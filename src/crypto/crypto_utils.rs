pub(crate) struct CryptoUtil {
}

impl CryptoUtil {
    /// Returns (gcd, u, v) such that u*a + v*b = gcd in Z
    pub(crate) fn egcd(a: i64, b: i64) -> (i64, i64, i64) {
        if a == 0 {
            (b, 0, 1)
        } else {
            let (gcd, u1, v1) = Self::egcd(b % a, a);
            (gcd, v1 - (b / a) * u1, u1)
        }
    }
}
