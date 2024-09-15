#![feature(portable_simd)]

use core::f64;
use std::simd::{f64x4, num::SimdFloat};

mod coeff;

const RAD: f64 = 648000.0 / f64::consts::PI; 
const J2000: f64 = 2451545.0;


fn u64_as_f64_slice(data: &[u64]) -> &[f64] {
    unsafe {
        std::slice::from_raw_parts(
            data.as_ptr() as *const f64,
            data.len(),
        )
    }
}

fn accumulate_main(data_u64: &[u64], tj: f64, tol: f64) -> (f64, f64) {
    let data = u64_as_f64_slice(data_u64);
    let n_blocks = data.len() / 6;
    let tj2 = tj * tj;
    let tj3 = tj * tj2;
    let tj4 = tj * tj3;

    let tj_vec = f64x4::from_array([tj, tj2, tj3, tj4]);
    let tj_1_vec = f64x4::from_array([1.0, tj, tj2, tj3]);
    let k_vec = f64x4::from_array([1.0, 2.0, 3.0, 4.0]);

    let mut res: f64 = 0.0;
    let mut res_p: f64 = 0.0;

    let mut blk_start_idx= 0;
    let mut blk_end_idx = n_blocks;

    loop {
        let blk_mid_idx = (blk_start_idx + blk_end_idx) / 2;
        if data[blk_mid_idx * 6].abs() > tol {
            blk_end_idx = blk_mid_idx;
        } else {
            blk_start_idx = blk_mid_idx;
        }
        if blk_end_idx - blk_start_idx < 2 {
            break
        }
    }

    let mut idx = blk_start_idx * 6;

    loop {
        let x = data[idx]; // cmpb[n];

        let mut y = data[idx+1]; // fmpb[0, n];
        let mut yp = 0.0;

        let fmpb_1_5_vec = f64x4::from_slice(&data[idx+2..idx+6]);
        y += (fmpb_1_5_vec * tj_vec).reduce_sum();
        yp += (k_vec * fmpb_1_5_vec * tj_1_vec).reduce_sum();

        res += x * y.sin();
        res_p += x * yp * y.cos();
        
        idx += 6;
        if idx >= data.len() {
            break
        }
    }

    (res, res_p)
}

fn accumulate_perturb(data_u64: &[u64], tj: f64, it: i32, tol: f64) -> (f64, f64) {
    let data = u64_as_f64_slice(data_u64);
    let n_blocks = data.len() / 6;

    let tj2 = tj * tj;
    let tj3 = tj * tj2;
    let tj4 = tj * tj3;

    let tj_vec = f64x4::from_array([tj, tj2, tj3, tj4]);
    let tj_1_vec = f64x4::from_array([1.0, tj, tj2, tj3]);
    let k_vec = f64x4::from_array([1.0, 2.0, 3.0, 4.0]);

    let tit = match it {
        1 => tj,
        2 => tj2,
        3 => tj3,
        _ => 1.0,
    };

    let tit_1 = match it {
        1 => 0.0,
        2 => tj,
        3 => tj2,
        _ => 1.0,
    };


    let mut res: f64 = 0.0;
    let mut res_p: f64 = 0.0;


    let mut blk_start_idx= 0;
    let mut blk_end_idx = n_blocks;

    loop {
        let blk_mid_idx = (blk_start_idx + blk_end_idx) / 2;
        if data[blk_mid_idx * 6].abs() > tol {
            blk_end_idx = blk_mid_idx;
        } else {
            blk_start_idx = blk_mid_idx;
        }
        if blk_end_idx - blk_start_idx < 2 {
            break
        }
    }

    let mut idx = blk_start_idx * 6;

    loop {
        let x = data[idx]; // cper[n];

        let mut y: f64 = data[idx+1]; // fper[0, n];
        let xp: f64 = (it as f64) * x * tit_1;

        let fper_1_5_vec = f64x4::from_slice(&data[idx+2..idx+6]);
        y += (fper_1_5_vec * tj_vec).reduce_sum();
        let yp = (k_vec * fper_1_5_vec * tj_1_vec).reduce_sum();

        res += x * tit * y.sin();
        res_p += xp * y.sin() + x * tit * yp * y.cos();
        idx += 6;

        if idx >= data.len() {
            break
        }
    }

    (res, res_p)
}


/// Calculate geocentric spherical coordinates of the Moon. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `cmpb` or `cper` are smaller than `tol`.
///
/// # Return 
///  * out.0: longitude  (rad)
///  * out.1: latitude   (rad)
///  * out.2: distance   (km)
///  * out.3: longitude' (rad/day)
///  * out.4: latitude'  (rad/day)
///  * out.5: distance'  (km/day)
pub fn spherical(t:f64, tol:f64) -> (f64, f64, f64, f64, f64, f64) {
    let tj = (t - J2000) / coeff::CONST_SC;
    let mut lng: f64 = 0.0; 
    let mut lat: f64 = 0.0; 
    let mut dst: f64 = 0.0; 
    let mut v_lng: f64 = 0.0; 
    let mut v_lat: f64 = 0.0; 
    let mut v_dst: f64 = 0.0; 

    // Main problem
    let src = &coeff::DATA_COEFF_MAIN_0;
    let (d, dv) = accumulate_main(src, tj, tol);
    lng += d; v_lng += dv;

    let src = &coeff::DATA_COEFF_MAIN_1;
    let (d, dv) = accumulate_main(src, tj, tol);
    lat += d; v_lat += dv;

    let src = &coeff::DATA_COEFF_MAIN_2;
    let (d, dv) = accumulate_main(src, tj, tol);
    dst += d; v_dst += dv;


    // Perturbation for lng
    let src = &coeff::DATA_COEFF_PERTURB_0_0;
    let (d, dv) = accumulate_perturb(src, tj, 0, tol);
    lng += d; v_lng += dv;

    let src = &coeff::DATA_COEFF_PERTURB_0_1;
    let (d, dv) = accumulate_perturb(src, tj, 1, tol);
    lng += d; v_lng += dv;

    let src = &coeff::DATA_COEFF_PERTURB_0_2;
    let (d, dv) = accumulate_perturb(src, tj, 2, tol);
    lng += d; v_lng += dv;

    let src = &coeff::DATA_COEFF_PERTURB_0_3;
    let (d, dv) = accumulate_perturb(src, tj, 3, tol);
    lng += d; v_lng += dv;

    // Perturbation for lat
    let src = &coeff::DATA_COEFF_PERTURB_1_0;
    let (d, dv) = accumulate_perturb(src, tj, 0, tol);
    lat += d; v_lat += dv;

    let src = &coeff::DATA_COEFF_PERTURB_1_1;
    let (d, dv) = accumulate_perturb(src, tj, 1, tol);
    lat += d; v_lat += dv;

    let src = &coeff::DATA_COEFF_PERTURB_1_2;
    let (d, dv) = accumulate_perturb(src, tj, 2, tol);
    lat += d; v_lat += dv;


    // Perturbation for dist
    let src = &coeff::DATA_COEFF_PERTURB_2_0;
    let (d, dv) = accumulate_perturb(src, tj, 0, tol);
    dst += d; v_dst += dv;

    let src = &coeff::DATA_COEFF_PERTURB_2_1;
    let (d, dv) = accumulate_perturb(src, tj, 1, tol);
    dst += d; v_dst += dv;

    let src = &coeff::DATA_COEFF_PERTURB_2_2;
    let (d, dv) = accumulate_perturb(src, tj, 2, tol);
    dst += d; v_dst += dv;

    let src = &coeff::DATA_COEFF_PERTURB_2_3;
    let (d, dv) = accumulate_perturb(src, tj, 3, tol);
    dst += d; v_dst += dv;

    (lng, lat, dst, v_lng, v_lat, v_dst)
}


/// Calculate geocentric rectangular coordinates of the Moon. The frame is dynamical
/// equinox and ecliptic J2000. The time scale of `t` is TDB. To speed up, we ignore
/// terms whose `cmpb` or `cper` are smaller than `tol`.
///
/// # Return 
///  * out.0: X  (km)
///  * out.1: Y  (km)
///  * out.2: Z  (km)
///  * out.3: X' (km/day)
///  * out.4: Y' (km/day)
///  * out.5: Z' (km/day)
pub fn cartesian(t: f64, tol: f64) -> (f64, f64, f64, f64, f64, f64) {
    let mut v = spherical(t, tol);

    let tj = (t - J2000) / coeff::CONST_SC;
    let tj2 = tj * tj;
    let tj3 = tj2 * tj;
    let tj4 = tj3 * tj;

    v.0 = v.0 / RAD  
        + coeff::CONST_W[0] 
        + coeff::CONST_W[1] * tj
        + coeff::CONST_W[2] * tj2
        + coeff::CONST_W[3] * tj3
        + coeff::CONST_W[4] * tj4;
    v.1 = v.1 / RAD;
    v.2 = v.2 * coeff::CONST_A405 / coeff::CONST_AELP;
    v.3 = v.3 / RAD
        + coeff::CONST_W[1] 
        + 2.0 * coeff::CONST_W[2] * tj
        + 3.0 * coeff::CONST_W[3] * tj2
        + 4.0 * coeff::CONST_W[4] * tj3;
    v.4 = v.4 / RAD;

    let mut xyz: (f64, f64, f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    let clamb = v.0.cos();
    let slamb = v.0.sin();
    let cbeta = v.1.cos();
    let sbeta = v.1.sin();
    let cw = v.2 * cbeta;
    let sw = v.2 * sbeta;

    let x1 = cw * clamb;
    let x2 = cw * slamb;
    let x3 = sw;
    let xp1    = (v.5*cbeta-v.4*sw)*clamb-v.3*x2;
    let xp2    = (v.5*cbeta-v.4*sw)*slamb+v.3*x1;
    let xp3    = v.5*sbeta+v.4*cw;
    
    let pw     = (
        coeff::CONST_P_POLY[0]
        +coeff::CONST_P_POLY[1]*tj
        +coeff::CONST_P_POLY[2]*tj2
        +coeff::CONST_P_POLY[3]*tj3
        +coeff::CONST_P_POLY[4]*tj4)*tj;

    let qw     = (
        coeff::CONST_Q_POLY[0]
        +coeff::CONST_Q_POLY[1]*tj
        +coeff::CONST_Q_POLY[2]*tj2
        +coeff::CONST_Q_POLY[3]*tj3
        +coeff::CONST_Q_POLY[4]*tj4)*tj;

    let ra     = 2.0*((1.0-pw*pw-qw*qw).sqrt());
    let pwqw   = 2.0*pw*qw;
    let pw2    = 1.0-2.0*pw*pw;
    let qw2    = 1.0-2.0*qw*qw;
    let pwra   = pw*ra;
    let qwra   = qw*ra;
    
    xyz.0 = pw2*x1+pwqw*x2+pwra*x3;
    xyz.1 = pwqw*x1+qw2*x2-qwra*x3;
    xyz.2 = -pwra*x1+qwra*x2+(pw2+qw2-1.0)*x3;
    
    let ppw     = coeff::CONST_P_POLY[0] + (
        2.0*coeff::CONST_P_POLY[1]
        +3.0*coeff::CONST_P_POLY[2]*tj
        +4.0*coeff::CONST_P_POLY[3]*tj2
        +5.0*coeff::CONST_P_POLY[4]*tj3)*tj;
    let qpw    = coeff::CONST_Q_POLY[0] + (
        2.0*coeff::CONST_Q_POLY[1]
        +3.0*coeff::CONST_Q_POLY[2]*tj
        +4.0*coeff::CONST_Q_POLY[3]*tj2
        +5.0*coeff::CONST_Q_POLY[4]*tj3)*tj;

    let ppw2   = -4.0*pw*ppw;
    let qpw2   = -4.0*qw*qpw;
    let ppwqpw = 2.0*(ppw*qw+pw*qpw);
    let rap    = (ppw2+qpw2)/ra;
    let ppwra  = ppw*ra+pw*rap;
    let qpwra  = qpw*ra+qw*rap;
    
    xyz.3 = (pw2*xp1+pwqw*xp2+pwra*xp3 + ppw2*x1+ppwqpw*x2+ppwra*x3)/coeff::CONST_SC;
    xyz.4 = (pwqw*xp1+qw2*xp2-qwra*xp3 + ppwqpw*x1+qpw2*x2-qpwra*x3)/coeff::CONST_SC;
    xyz.5 = (-pwra*xp1+qwra*xp2+(pw2+qw2-1.0)*xp3 - ppwra*x1+qpwra*x2+(ppw2+qpw2)*x3)/coeff::CONST_SC;

    xyz
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cartesian() {
        let tj0=2444239.5;
        let dtj=2000.0;
        for n in 1..6 {
            let t = tj0 + ((n-1) as f64) * dtj;
            println!("JD {t}: {:#?}", cartesian(t, 0.0));
        }
    }

    fn is_matched(
        a: &(f64, f64, f64, f64, f64, f64),
        b: &(f64, f64, f64, f64, f64, f64),
    ) -> bool {
        let tol = 1e-4;
        if (a.0 - b.0).abs() > tol { return false; }
        if (a.1 - b.1).abs() > tol { return false; }
        if (a.2 - b.2).abs() > tol { return false; }
        if (a.3 - b.3).abs() > tol { return false; }
        if (a.4 - b.4).abs() > tol { return false; }
        if (a.5 - b.5).abs() > tol { return false; }

        return true;
    }
    #[test]
    fn check_values() {
        let cases = vec![
            (2444239.5, 
                (43890.2824005, 381188.7274523, -31633.3816524, -87516.1974842, 13707.6644351, 2754.2212424),
            ),
            (2446239.5, 
                (-313664.5964499, 212007.2667385, 33744.7512039, -47315.9128069, -75710.8750091, -1475.6286887),
            ),
            (2448239.5, 
                (-273220.0606714, -296859.7682229, -34604.3569962, 60542.3275903, -58162.3167425, 2270.8869053),
            ),
            (2450239.5, 
                (171613.1427993, -318097.3375025, 31293.5482404, 83266.7799036, 42585.8302849, -1695.8261107),
            ),
            (2452239.5, 
                (396530.0063512, 47487.9224886, -36085.3090343, -12664.2869360, 83512.7571912, 1507.3675566),
            ),
        ];

        let mut failcnt = 0;
        for (t, expected) in cases.iter() {
            let actual = cartesian(*t, 0.0);

            let mut failed = false;
            if !is_matched(expected, &actual) {
                failed = true;
            }
            if !failed {
                println!("[t={t}] PASSED")
            } else {
                failcnt += 1; 
                println!("[t={t}] FAILED")
            }
        }

        assert_eq!(failcnt, 0);
    }

}