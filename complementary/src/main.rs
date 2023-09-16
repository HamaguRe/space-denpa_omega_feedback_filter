//! 相補フィルタにより姿勢推定を行う
//! 
//! 推定してるのは四元数のみ．

use std::fs;
use std::io::{Write, BufWriter};
use std::mem::MaybeUninit;

use rand::distributions::{Distribution, Normal};
use quaternion_core as quat;

mod ahrs;

const DT: f64 = 0.02;
const SIM_TIME: f64 = 30.0;
const N: usize = (SIM_TIME / DT) as usize + 1;

/// 角速度センサのノイズ分散
const GYR_VAR: f64 = 0.0001;

/// 加速度センサのノイズ分散
const ACC_VAR: f64 = 0.01;

/// 地磁気センサのノイズ分散
const MAG_VAR: f64 = 0.01;

fn main() {
    // CSVファイルにデータ保存（同一ファイルが存在したら上書き）
    let mut file = BufWriter::new( fs::File::create("result.csv").unwrap() );

    // 標準正規分布の乱数を生成
    let randn = Normal::new(0.0, 1.0);  // 平均値:0，標準偏差:1

    // 姿勢推定フィルタ
    let mut filter = ahrs::AttitudeFilter::new(1.0);

    let mut q = (1.0, [0.0; 3]);
    let gyr_bias = [-0.02, 0.01, 0.05];

    // ---- Loop start ---- //
    let gyr = [0.1; 3];
    for t in 0..N {
        // 積分（q = q + 0.5*Δt*q*ω）
        q = {
            let tmp0 = quat::scale(q.0, gyr);
            let dot = quat::dot(q.1, gyr);
            let cross = quat::cross(q.1, gyr);
            let tmp1 = (-dot, quat::add(tmp0, cross));
            quat::scale_add(0.5 * DT, tmp1, q)
        };
        q = quat::normalize(q);

        // 計測値生成
        let mut acc_b = quat::frame_rotation(q, ahrs::ACC_R);
        let mut mag_b = quat::frame_rotation(q, ahrs::MAG_R);
        acc_b = add_noise(&randn, ACC_VAR, acc_b);
        mag_b = add_noise(&randn, MAG_VAR, mag_b);

        // 推定
        let gyr_noisy = add_noise(&randn, GYR_VAR, gyr);
        filter.predict( quat::add(gyr_noisy, gyr_bias) );
        filter.filtering(acc_b, mag_b);

        // ---------- データ書き込み ---------- //
        // 時刻
        file.write( format!("{:.3},", t as f64 * DT ).as_bytes() ).unwrap();
        // オイラー角の真値
        let ypr_true = quat::to_euler_angles(quat::RotationType::Intrinsic, quat::RotationSequence::ZYX, q);
        for i in 0..3 {
            file.write( format!("{:.7},", ypr_true[i] ).as_bytes() ).unwrap();
        }
        // オイラー角の推定値
        let ypr_hat = quat::to_euler_angles(quat::RotationType::Intrinsic, quat::RotationSequence::ZYX, filter.q);
        for i in 0..3 {
            file.write( format!("{:.7},", ypr_hat[i] ).as_bytes() ).unwrap();
        }
        // 角速度バイアスの真値
        for i in 0..3 {
            file.write( format!("{:.7},", gyr_bias[i] ).as_bytes() ).unwrap();
        }
        // 四元数の真値
        file.write( format!("{:.7},", q.0 ).as_bytes() ).unwrap();
        for i in 0..3 {
            file.write( format!("{:.7},", q.1[i] ).as_bytes() ).unwrap();
        }
        // 四元数の推定値
        file.write( format!("{:.7},", filter.q.0 ).as_bytes() ).unwrap();
        for i in 0..2 {
            file.write( format!("{:.7},", filter.q.1[i] ).as_bytes() ).unwrap();
        }

        file.write( format!("{:.7}\n", filter.q.1[2] ).as_bytes() ).unwrap();
        // ------------------------------------ //
    }
}

/// ベクトルxにノイズを加える．
fn add_noise(randn: &rand::distributions::Normal, variance: f64, x: quat::Vector3<f64>) -> quat::Vector3<f64> {
    let mut noisy: quat::Vector3<f64> = unsafe {MaybeUninit::uninit().assume_init()};

    let tmp = variance.sqrt();
    for i in 0..3 {
        noisy[i] = x[i] + randn.sample(&mut rand::thread_rng()) * tmp;
    }
    noisy
}