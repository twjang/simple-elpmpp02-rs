use include_bytes_plus::include_bytes;

pub static DATA_COEFF_MAIN_0: [u64; 6138] = include_bytes!("src/bins/coeff_main_0.bin" as u64le);
pub static DATA_COEFF_MAIN_1: [u64; 5508] = include_bytes!("src/bins/coeff_main_1.bin" as u64le);
pub static DATA_COEFF_MAIN_2: [u64; 4224] = include_bytes!("src/bins/coeff_main_2.bin" as u64le);
pub static DATA_COEFF_PERTURB_0_0: [u64; 67884] = include_bytes!("src/bins/coeff_perturb_0_0.bin" as u64le);
pub static DATA_COEFF_PERTURB_0_1: [u64; 7194] = include_bytes!("src/bins/coeff_perturb_0_1.bin" as u64le);
pub static DATA_COEFF_PERTURB_0_2: [u64; 1314] = include_bytes!("src/bins/coeff_perturb_0_2.bin" as u64le);
pub static DATA_COEFF_PERTURB_0_3: [u64; 12] = include_bytes!("src/bins/coeff_perturb_0_3.bin" as u64le);
pub static DATA_COEFF_PERTURB_1_0: [u64; 38772] = include_bytes!("src/bins/coeff_perturb_1_0.bin" as u64le);
pub static DATA_COEFF_PERTURB_1_1: [u64; 3096] = include_bytes!("src/bins/coeff_perturb_1_1.bin" as u64le);
pub static DATA_COEFF_PERTURB_1_2: [u64; 312] = include_bytes!("src/bins/coeff_perturb_1_2.bin" as u64le);
pub static DATA_COEFF_PERTURB_2_0: [u64; 72690] = include_bytes!("src/bins/coeff_perturb_2_0.bin" as u64le);
pub static DATA_COEFF_PERTURB_2_1: [u64; 6990] = include_bytes!("src/bins/coeff_perturb_2_1.bin" as u64le);
pub static DATA_COEFF_PERTURB_2_2: [u64; 1260] = include_bytes!("src/bins/coeff_perturb_2_2.bin" as u64le);
pub static DATA_COEFF_PERTURB_2_3: [u64; 12] = include_bytes!("src/bins/coeff_perturb_2_3.bin" as u64le);
pub const CONST_A405: f64 = 384747.9613701725;
pub const CONST_AELP: f64 = 384747.980674318;
pub const CONST_W: [f64;15] = [
   3.81034392032190893929e+00,  8.39968473020743294910e+03, -3.31919929752746043178e-05,
   3.20170955004737532337e-08, -1.53637455543611967059e-10,  1.45478934807000870322e+00,
   7.09933054823066100880e+01, -1.85504743616622859046e-04, -2.18394018929412647762e-07,
   1.03270162213142248579e-09,  2.18243867555731885233e+00, -3.37814274330531461032e+01,
   3.08293019817553899642e-05,  3.69670431846021155713e-08, -1.73854186045879597399e-10
];
pub const CONST_P_POLY: [f64;5] = [
   1.01803910000000003586e-05,  4.70204390000000004823e-07, -5.41736700000000015858e-10,
  -2.50794800000000008912e-12,  4.63485999999999974455e-15
];
pub const CONST_Q_POLY: [f64;5] = [
  -1.13469001999999999472e-04,  1.23726739999999989995e-07,  1.26541700000000003609e-09,
  -1.37180800000000004512e-12, -3.20334000000000010177e-15
];
pub const CONST_SC: f64 = 36525.0;