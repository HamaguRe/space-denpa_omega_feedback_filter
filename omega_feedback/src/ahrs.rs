//! 姿勢推定フィルタ

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
    pub q: Quaternion<f64>,      // 姿勢推定値
    pub q_lpf: Quaternion<f64>,  // 一つ前のLPF出力
    gyr_correct: Vector3<f64>,   // 補正角速度（角速度バイアスの推定値を含む）
    coef_lpf_term1: f64,         // LPF第一項の係数
    coef_lpf_term2: f64,         // LPF第二項の係数
    coef_gyr_c: f64,             // 補正角速度を計算するときのパラメータ
    coef_integ: f64,             // 補正角速度の積分係数
    pub gyr_integ: Vector3<f64>, // 補正角速度の積分項
}

impl AttitudeFilter {
    /// * cutoff: ローパスフィルタのカットオフ周波数[Hz]
    /// * alpha : 基準姿勢に収束するまでの時間[s]
    /// * beta  : 補正角速度の積分係数
    pub fn new(cutoff: f64, alpha: f64, beta: f64) -> Self {
        let k = 1.0 / (2.0 * PI * DT * cutoff + 1.0);
        Self {
            q: (1.0, [0.0; 3]),
            q_lpf: (1.0, [0.0; 3]),
            gyr_correct: [0.0; 3],
            coef_lpf_term1: k,
            coef_lpf_term2: 1.0 - k,
            coef_gyr_c: 2.0 / alpha,
            coef_integ: beta,
            gyr_integ: [0.0; 3],
        }
    }

    /// 予測ステップ
    /// 
    /// * gyr: 機体上で計測した角速度[rad/s]
    pub fn predict(&mut self, gyr: Vector3<f64>) {
        let omega = quat::add(gyr, self.gyr_correct);

        // 積分（q[n+1] = q[n] + Δt/2 *q[n]*ω[n]）
        let tmp0 = quat::scale(self.q.0, omega);
        let dot = quat::dot(self.q.1, omega);
        let cross = quat::cross(self.q.1, omega);
        let tmp1 = (-dot, quat::add(tmp0, cross));
        self.q = quat::scale_add(0.5 * DT, tmp1, self.q);
        // 正規化
        self.q = quat::normalize(self.q);
    }

    /// 補正ステップ
    /// 
    /// * acc: 機体上のセンサで計測した加速度[m/s^2]
    /// * mag: 機体上のセンサで計測した地磁気（方向だけわかれば良いので単位不問）
    pub fn correct(&mut self, acc: Vector3<f64>, mag: Vector3<f64>) {
        // accとmagから姿勢q_gmを計算
        let q_g = quat::rotate_a_to_b(acc, ACC_R).unwrap();
        let mag_b2r = quat::hadamard(quat::point_rotation(q_g, mag), [1.0, 1.0, 0.0]);
        let q_e = quat::rotate_a_to_b_shortest(mag_b2r, MAG_R, 1.0).unwrap();
        let q_gm = quat::mul(q_e, q_g);

        // LPFを通す（q_gmの符号をq_lpfに合わせる）
        let term2_coef = if quat::dot(q_gm, self.q_lpf).is_sign_negative() {
            -self.coef_lpf_term2
        } else {
            self.coef_lpf_term2
        };
        let term1 = quat::scale(self.coef_lpf_term1, self.q_lpf);
        let term2 = quat::scale(term2_coef, q_gm);
        self.q_lpf = quat::normalize( quat::add(term1, term2) );

        // qからq_lpfに到達するための角速度を計算
        let term1 = quat::scale(self.q.0, self.q_lpf.1);
        let term2 = quat::scale(self.q_lpf.0, self.q.1);
        let term3 = quat::cross(self.q_lpf.1, self.q.1);
        self.gyr_correct = quat::scale(self.coef_gyr_c, quat::add(quat::sub(term1, term2), term3));
        // 符号を合わせる（LPFの前で合わせてるから、シミュレーションなら無くても平気。実環境では必要だと思う）
        if quat::dot(self.q, self.q_lpf).is_sign_negative() {
            self.gyr_correct = quat::negate(self.gyr_correct);
        }

        // 積分項を更新
        self.gyr_integ = quat::scale_add(DT, self.gyr_correct, self.gyr_integ);

        // 積分項の値を補正角速度に反映
        self.gyr_correct = quat::scale_add(self.coef_integ, self.gyr_integ, self.gyr_correct);
    }
}
