//! EKFを用いた姿勢推定アルゴリズムを実装
//! 
//! うまく推定できないことが多い。
//! 何回か実行していると綺麗に推定できるときがある。
//! 
//! シミュレーション開始直後は角速度バイアスの推定値が3rad/s程度まで跳ね上がる。

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

    // 誤差共分散行列の初期値
    let mut p = [[0.0; 7]; 7];
    for i in 0..7 {
        p[i][i] = 1.0;
    }

    // 入力行列
    let mut g = [[0.0; 3]; 7];
    for i in 0..3 {
        g[i+4][i] = DT;
    }

    // システムノイズの共分散行列
    let mut q = [[0.0; 3]; 3];
    for i in 0..3 {
        q[i][i] = GYR_VAR;
    }
    // 観測ノイズの分散（共分散行列の対角成分）
    let mut r = [0.0; 6];
    for i in 0..3 {
        r[i]   = ACC_VAR;  // 加速度
        r[i+3] = MAG_VAR;  // 地磁気
    }

    let mut filter = ahrs::ExUdFilter::new(p, g, q, r);
    filter.x[0] = 1.0;  // 初期値は恒等四元数

    let mut x = filter.x;
    let mut y_true = ahrs::calc_h(x);
    let mut y = y_true;

    let gyr_bias = [-0.02, 0.01, 0.05];

    // ---- Loop start ---- //
    let gyr = [0.1; 3];
    for t in 0..N {
        x = ahrs::calc_f(x, gyr);
        y_true = ahrs::calc_h(x);
        for i in 0..y.len() {
            y[i] = y_true[i] + randn.sample(&mut rand::thread_rng()) * r[i].sqrt();
        }

        // 推定
        let gyr_noisy = add_noise(&randn, GYR_VAR, gyr);
        filter.predict( quat::add_vec(gyr_noisy, gyr_bias) );
        filter.filtering(&y);

        // ---------- データ書き込み ---------- //
        // 時刻
        file.write( format!("{:.3},", t as f64 * DT ).as_bytes() ).unwrap();
        // オイラー角の真値
        let q = (x[0], [x[1], x[2], x[3]]);
        let ypr_true = quat::to_euler_angles(q);
        for i in 0..3 {
            file.write( format!("{:.7},", ypr_true[i] ).as_bytes() ).unwrap();
        }
        // オイラー角の推定値
        let q_hat = (filter.x[0], [filter.x[1], filter.x[2], filter.x[3]]);
        let ypr_hat = quat::to_euler_angles(q_hat);
        for i in 0..3 {
            file.write( format!("{:.7},", ypr_hat[i] ).as_bytes() ).unwrap();
        }
        // 角速度バイアスの真値
        for i in 0..3 {
            file.write( format!("{:.7},", gyr_bias[i] ).as_bytes() ).unwrap();
        }
        // 角速度バイアスの推定値
        let gyr_bias_hat = [filter.x[4], filter.x[5], filter.x[6]];
        for i in 0..3 {
            file.write( format!("{:.7},", gyr_bias_hat[i] ).as_bytes() ).unwrap();
        }
        // 四元数の真値
        file.write( format!("{:.7},", q.0 ).as_bytes() ).unwrap();
        for i in 0..3 {
            file.write( format!("{:.7},", q.1[i] ).as_bytes() ).unwrap();
        }
        // 四元数の推定値
        file.write( format!("{:.7},", q_hat.0 ).as_bytes() ).unwrap();
        for i in 0..2 {
            file.write( format!("{:.7},", q_hat.1[i] ).as_bytes() ).unwrap();
        }
        file.write( format!("{:.7}\n", q_hat.1[2] ).as_bytes() ).unwrap();
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