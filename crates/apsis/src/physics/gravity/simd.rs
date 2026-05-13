//! AVX2 leaf-pair Plummer kernel for the Barnes-Hut walk.
//!
//! Vectorises the leaf-pair phase of the two-phase walk
//! ([`bh_process_lists`](super::engine::bh_process_lists)) via
//! `std::arch::x86_64` AVX2 + FMA intrinsics. Scattered leaf-mate body
//! indices are loaded with `_mm256_i32gather_pd` against the SoA
//! [`BodyArrays`](crate::domain::body_arrays::BodyArrays) field vectors.
//! The accepted-node phase remains scalar.

#![allow(dead_code)]

use std::arch::x86_64::*;

use crate::domain::body_arrays::BodyArrays;
use crate::math::Vec3;

use super::kernel::G;

// ── Tier 0: saxpy hardware sanity microbench ────────────────────────────────── //

/// Scalar baseline for the saxpy hardware sanity check.
///
/// `y[i] += a * x[i]` over `x.len().min(y.len())` lanes, no SIMD.
#[inline(never)]
pub(crate) fn saxpy_scalar(a: f64, x: &[f64], y: &mut [f64]) {
    let n = x.len().min(y.len());
    for i in 0..n {
        y[i] += a * x[i];
    }
}

/// AVX2 saxpy via `_mm256_fmadd_pd` (4 doubles per FMA instruction).
///
/// Caller must ensure AVX2 is available (`is_x86_feature_detected!("avx2")`).
/// Tier 0 gate predicts `t_scalar / t_avx2 ≥ 2.5×` on Zen 4.
#[target_feature(enable = "avx2,fma")]
pub(crate) unsafe fn saxpy_avx2(a: f64, x: &[f64], y: &mut [f64]) {
    let n = x.len().min(y.len());
    let n_chunks = n / 4;

    unsafe {
        let a_vec = _mm256_set1_pd(a);

        for chunk in 0..n_chunks {
            let off = chunk * 4;
            let xv = _mm256_loadu_pd(x.as_ptr().add(off));
            let yv = _mm256_loadu_pd(y.as_ptr().add(off));
            let res = _mm256_fmadd_pd(a_vec, xv, yv);
            _mm256_storeu_pd(y.as_mut_ptr().add(off), res);
        }
    }

    for i in (n_chunks * 4)..n {
        y[i] += a * x[i];
    }
}

// ── Tier 2a: kernel-isolated Plummer leaf-pair microbench ──────────────────── //

/// Per-interaction tuple format for the kernel-isolated microbench:
/// `(dx, dy, dz, eps2, other_mass)` — already pre-computed (no gather).
/// Measures pure arithmetic throughput, not load throughput.
#[inline]
pub(crate) fn plummer_kernel_scalar_micro(
    target_acc: &mut Vec3,
    interactions: &[(f64, f64, f64, f64, f64)],
) {
    for &(dx, dy, dz, eps2, other_mass) in interactions {
        let r_sq = dx * dx + dy * dy + dz * dz;
        let inv_r = (r_sq + eps2).sqrt().recip();
        let inv_r3 = inv_r * inv_r * inv_r;
        let fac = G * other_mass * inv_r3;
        target_acc.x += dx * fac;
        target_acc.y += dy * fac;
        target_acc.z += dz * fac;
    }
}

/// AVX2 Plummer kernel microbench, 4 interactions per chunk via aligned
/// SIMD ops on pre-laid-out arrays.
///
/// Caller must pass arrays of length divisible by 4 with
/// 32-byte alignment (AVX2 `_mm256_loadu_pd` accepts unaligned but
/// tighter alignment is faster). The microbench harness arranges this.
#[target_feature(enable = "avx2,fma")]
pub(crate) unsafe fn plummer_kernel_avx2_micro(
    target_acc: &mut Vec3,
    dx: &[f64],
    dy: &[f64],
    dz: &[f64],
    eps2: &[f64],
    other_mass: &[f64],
) {
    let n = dx.len();
    let n_chunks = n / 4;
    debug_assert_eq!(n % 4, 0, "microbench requires N divisible by 4");
    debug_assert_eq!(dy.len(), n);
    debug_assert_eq!(dz.len(), n);
    debug_assert_eq!(eps2.len(), n);
    debug_assert_eq!(other_mass.len(), n);

    unsafe {
        let g_vec = _mm256_set1_pd(G);
        let one = _mm256_set1_pd(1.0);
        let mut ax_acc = _mm256_setzero_pd();
        let mut ay_acc = _mm256_setzero_pd();
        let mut az_acc = _mm256_setzero_pd();

        for chunk in 0..n_chunks {
            let off = chunk * 4;
            let dx_v = _mm256_loadu_pd(dx.as_ptr().add(off));
            let dy_v = _mm256_loadu_pd(dy.as_ptr().add(off));
            let dz_v = _mm256_loadu_pd(dz.as_ptr().add(off));
            let eps2_v = _mm256_loadu_pd(eps2.as_ptr().add(off));
            let mass_v = _mm256_loadu_pd(other_mass.as_ptr().add(off));

            let dx_sq = _mm256_mul_pd(dx_v, dx_v);
            let dy_sq = _mm256_mul_pd(dy_v, dy_v);
            let dz_sq = _mm256_mul_pd(dz_v, dz_v);
            let r_sq = _mm256_add_pd(_mm256_add_pd(dx_sq, dy_sq), dz_sq);

            let r_sq_eff = _mm256_add_pd(r_sq, eps2_v);
            let r_eff = _mm256_sqrt_pd(r_sq_eff);
            let inv_r = _mm256_div_pd(one, r_eff);

            let inv_r2 = _mm256_mul_pd(inv_r, inv_r);
            let inv_r3 = _mm256_mul_pd(inv_r2, inv_r);

            let g_mass = _mm256_mul_pd(g_vec, mass_v);
            let fac = _mm256_mul_pd(g_mass, inv_r3);

            ax_acc = _mm256_fmadd_pd(dx_v, fac, ax_acc);
            ay_acc = _mm256_fmadd_pd(dy_v, fac, ay_acc);
            az_acc = _mm256_fmadd_pd(dz_v, fac, az_acc);
        }

        target_acc.x += horizontal_sum_pd(ax_acc);
        target_acc.y += horizontal_sum_pd(ay_acc);
        target_acc.z += horizontal_sum_pd(az_acc);
    }
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn horizontal_sum_pd(v: __m256d) -> f64 {
    let mut buf = [0.0_f64; 4];
    unsafe { _mm256_storeu_pd(buf.as_mut_ptr(), v) };
    buf[0] + buf[1] + buf[2] + buf[3]
}

// ── Production kernel: AVX2 leaf-pair (gather-based) ───────────────────────── //

/// AVX2-vectorised leaf-pair Plummer kernel, processing the
/// `leaf_body_indices` from the two-phase walk's interaction list.
///
/// Loads scattered fields via `_mm256_i32gather_pd` (no Morton ordering
/// in production — leaf-mate indices are arbitrary). Returns the
/// leaf-pair contribution to `(a, phi)`; the caller composes with the
/// scalar accepted-node contribution.
///
/// Caller must ensure AVX2+FMA available.
#[target_feature(enable = "avx2,fma")]
pub(crate) unsafe fn process_leafpair_avx2(
    body_pos_x: f64,
    body_pos_y: f64,
    body_pos_z: f64,
    body_softening: f64,
    arrays: &BodyArrays,
    leaf_body_indices: &[u32],
) -> (Vec3, f64) {
    let n = leaf_body_indices.len();
    let n_chunks = n / 4;

    let (mut a, mut phi) = unsafe {
        let bpx = _mm256_set1_pd(body_pos_x);
        let bpy = _mm256_set1_pd(body_pos_y);
        let bpz = _mm256_set1_pd(body_pos_z);
        let body_soft_sq = _mm256_set1_pd(body_softening * body_softening);
        let g_vec = _mm256_set1_pd(G);
        let half = _mm256_set1_pd(0.5);
        let one = _mm256_set1_pd(1.0);
        let neg_g = _mm256_set1_pd(-G);

        let mut ax_acc = _mm256_setzero_pd();
        let mut ay_acc = _mm256_setzero_pd();
        let mut az_acc = _mm256_setzero_pd();
        let mut phi_acc = _mm256_setzero_pd();

        for chunk in 0..n_chunks {
            let off = chunk * 4;
            let idx = _mm_loadu_si128(leaf_body_indices.as_ptr().add(off) as *const __m128i);

            let opx = _mm256_i32gather_pd::<8>(arrays.pos_x.as_ptr(), idx);
            let opy = _mm256_i32gather_pd::<8>(arrays.pos_y.as_ptr(), idx);
            let opz = _mm256_i32gather_pd::<8>(arrays.pos_z.as_ptr(), idx);
            let omass = _mm256_i32gather_pd::<8>(arrays.mass.as_ptr(), idx);
            let osoft = _mm256_i32gather_pd::<8>(arrays.softening.as_ptr(), idx);

            let dx = _mm256_sub_pd(opx, bpx);
            let dy = _mm256_sub_pd(opy, bpy);
            let dz = _mm256_sub_pd(opz, bpz);

            let r_sq = _mm256_add_pd(
                _mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)),
                _mm256_mul_pd(dz, dz),
            );

            let other_soft_sq = _mm256_mul_pd(osoft, osoft);
            let pair_eps2 = _mm256_mul_pd(_mm256_add_pd(body_soft_sq, other_soft_sq), half);

            let r_sq_eff = _mm256_add_pd(r_sq, pair_eps2);
            let r_eff = _mm256_sqrt_pd(r_sq_eff);
            let inv_r = _mm256_div_pd(one, r_eff);

            let inv_r2 = _mm256_mul_pd(inv_r, inv_r);
            let inv_r3 = _mm256_mul_pd(inv_r2, inv_r);

            let g_mass = _mm256_mul_pd(g_vec, omass);
            let fac = _mm256_mul_pd(g_mass, inv_r3);

            ax_acc = _mm256_fmadd_pd(dx, fac, ax_acc);
            ay_acc = _mm256_fmadd_pd(dy, fac, ay_acc);
            az_acc = _mm256_fmadd_pd(dz, fac, az_acc);

            let neg_g_mass = _mm256_mul_pd(neg_g, omass);
            phi_acc = _mm256_fmadd_pd(neg_g_mass, inv_r, phi_acc);
        }

        (
            Vec3::new(
                horizontal_sum_pd(ax_acc),
                horizontal_sum_pd(ay_acc),
                horizontal_sum_pd(az_acc),
            ),
            horizontal_sum_pd(phi_acc),
        )
    };

    for &raw_bi in &leaf_body_indices[n_chunks * 4..] {
        let bi = raw_bi as usize;
        let other_mass = arrays.mass[bi];
        let dx_s = arrays.pos_x[bi] - body_pos_x;
        let dy_s = arrays.pos_y[bi] - body_pos_y;
        let dz_s = arrays.pos_z[bi] - body_pos_z;
        let other_soft = arrays.softening[bi];
        let pair_eps2 = 0.5 * (body_softening * body_softening + other_soft * other_soft);
        let r_sq = dx_s * dx_s + dy_s * dy_s + dz_s * dz_s;
        let inv_r = (r_sq + pair_eps2).sqrt().recip();
        let inv_r3 = inv_r * inv_r * inv_r;
        let fac = G * other_mass * inv_r3;
        a.x += dx_s * fac;
        a.y += dy_s * fac;
        a.z += dz_s * fac;
        phi += -G * other_mass * inv_r;
    }

    (a, phi)
}
