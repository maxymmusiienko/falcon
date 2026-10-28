use crate::math::complex_numbers::Complex;
use crate::utils::falcon_config_scripts;

pub enum FalconDegree {
    N512 = 512,
    N1024 = 1024,
    Test4 = 4,
    Test8 = 8,
    q = 12289,
}

pub struct FalconContext {
    pub(crate) degree: usize,
    pub(crate) roots: Vec<Complex>,
    pub(crate) pos_roots: Vec<Complex>,
}

impl FalconContext {
    pub fn new(degree_enum: FalconDegree) -> Self {
        let degree = degree_enum as usize;
        let roots = falcon_config_scripts::find_complex_roots(degree);
        let pos_roots = falcon_config_scripts::find_complex_positive_image_roots(degree);

        FalconContext { degree, roots, pos_roots }
    }
}
