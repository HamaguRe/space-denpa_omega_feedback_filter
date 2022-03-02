//! 相補フィルタ

use super::DT;
use super::quat;
use super::quat::{Vector3, Quaternion};

const PI: f64 = std::f64::consts::PI;

/// 標準重力
pub const STANDARD_GRAVITY: f64 = 9.80665;

/// 基準座標系上における加速度計測値
pub const ACC_R: [f64; 3] = [0.0, 0.0, STANDARD_GRAVITY];

/// 基準座標系上における地磁気計測値
pub const MAG_R: [f64; 3] = [0.0, 1.0, 0.0];

pub struct AttitudeFilter {
    pub q: Quaternion<f64>,
    k: f64,
}

impl AttitudeFilter {
    pub fn new(cutoff: f64) -> Self {
        let tau = 1.0 / (2.0 * PI * cutoff);  // 時定数
        Self {
            q: (1.0, [0.0; 3]),
            k: tau / (tau + DT),
        }
    }

    pub fn predict(&mut self, gyr: Vector3<f64>) {
        // 積分（q = q + 0.5*Δt*q*ω）
        let tmp0 = quat::scale_vec(self.q.0, gyr);
        let dot = quat::dot_vec(self.q.1, gyr);
        let cross = quat::cross_vec(self.q.1, gyr);
        let tmp1 = (-dot, quat::add_vec(tmp0, cross));
        self.q = quat::scale_add(0.5 * DT, tmp1, self.q);
        // 正規化
        self.q = quat::normalize(self.q);
    }

    pub fn filtering(&mut self, acc: Vector3<f64>, mag: Vector3<f64>) {
        // accとmagから姿勢を計算
        let q_g = quat::rotate_a_to_b(acc, ACC_R);
        let mag_b_to_r = quat::hadamard_vec(quat::vector_rotation(q_g, mag), [1.0, 1.0, 0.0]);
        let q_e = quat::rotate_a_to_b(mag_b_to_r, MAG_R);
        let q_gm = quat::mul(q_e, q_g);
        
        self.q = quat::normalize( quat::lerp(self.q, q_gm, 1.0 - self.k) );
    }
}
