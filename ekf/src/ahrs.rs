//! AHRS

use std::mem::MaybeUninit;

use super::quat;
use super::quat::Vector3;
use super::DT;

const FRAC_DT_2: f64 = DT / 2.0;

/// ジャイロバイアスの係数（対角成分）
const BETA_DIAG: Vector3<f64> = [0.0; 3];  // 一定と仮定したので0

/// 標準重力
pub const STANDARD_GRAVITY: f64 = 9.80665;

/// 基準座標系上における加速度計測値
pub const ACC_R: [f64; 3] = [0.0, 0.0, STANDARD_GRAVITY];

/// 基準座標系上における地磁気計測値
pub const MAG_R: [f64; 3] = [0.0, 1.0, 0.0];

// --------- 入出力数 --------- //
// 以下の状態空間モデルを考える．
// x[k+1] = F*x[k] + G*w[k],
// y[k]   = H*x[k] + v[k].
// このとき，行列F, G, Hのサイズはそれぞれ
// F: N×N, G: N×M, H: P×N
// となる．

/// 状態変数の個数
const SYS_N: usize = 7;

/// 入力数
const SYS_M: usize = 3;

/// 出力数
const SYS_P: usize = 6;
// ---------------------------- //

// -- ベクトル・行列型の定義 -- //
// Vector○: X次元ベクトル
type VectorN<T>  = [T; SYS_N];
type VectorP<T>  = [T; SYS_P];
// Matrix○x□: ○行□列行列
type MatrixNxN<T>  = [[T; SYS_N]; SYS_N];
type MatrixNxM<T>  = [[T; SYS_M]; SYS_N];
type MatrixPxN<T>  = [[T; SYS_N]; SYS_P];
type MatrixMxM<T>  = [[T; SYS_M]; SYS_M];
type MatrixNxNM<T> = [[T; SYS_N + SYS_M]; SYS_N];  // N×(N+M)
// ---------------------------- //


/// U-D分解フィルタ
#[allow(non_snake_case)]
pub struct ExUdFilter {
    pub x: VectorN<f64>,    // 状態変数
    pub U: MatrixNxN<f64>,  // U-D分解した共分散行列
    F: MatrixNxN<f64>,  // システム行列
    G: MatrixNxM<f64>,  // 入力行列
    H: MatrixPxN<f64>,  // 出力行列
    R: VectorP<f64>,    // 観測ノイズの共分散行列の対角成分
}

impl ExUdFilter {
    /// 誤差共分散行列の初期値を零行列にするとU-D分解に失敗するので，
    /// スカラー行列にするのが無難．
    /// 
    /// 状態変数は全て零で初期化する．状態変数xはパブリックメンバにしているので，
    /// 零以外にしたい場合は構造体を作った後に適宜アクセスして書き換えること．
    #[allow(non_snake_case)]
    pub fn new(
        P: MatrixNxN<f64>,  // 共分散行列の初期値
        G: MatrixNxM<f64>,  // 入力行列
        Q: MatrixMxM<f64>,  // システムノイズの共分散行列
        R: VectorP<f64>     // 観測ノイズの共分散行列の対角成分
    ) -> Self {
        // システムノイズの分散Qが単位行列で無い場合には，
        // QをQ=C*C^Tと分解し，Gを改めてG*Cとおく．
        let c = cholesky_decomp(Q);
        let mut gc: MatrixNxM<f64> = unsafe {MaybeUninit::uninit().assume_init()};
        for i in 0..SYS_N {
            for j in 0..SYS_M {
                let mut sum = 0.0;
                for k in 0..SYS_M {
                    sum += G[i][k] * c[k][j];
                }
                gc[i][j] = sum;
            }
        }
        Self {
            x: [0.0; SYS_N],
            U: ud_decomp(P),
            F: [[0.0; SYS_N]; SYS_N],
            G: gc,
            H: [[0.0; SYS_N]; SYS_P],
            R: R,
        }
    }

    /// 予測ステップ
    /// 
    /// x(k+1) = F * x(k)
    /// P_bar = F*P*F^T + G*Q*Q^T
    pub fn predict(&mut self, gyr: Vector3<f64>) {
        // Working array
        let mut qq: VectorN<f64>    = unsafe {MaybeUninit::uninit().assume_init()};
        let mut z:  VectorN<f64>    = unsafe {MaybeUninit::uninit().assume_init()};
        let mut w:  MatrixNxNM<f64> = unsafe {MaybeUninit::uninit().assume_init()};

        // 線形化
        self.calc_jacobian_f(gyr);
        self.x = calc_f(self.x, gyr);

        // qqとwの左NxN要素を初期化
        for j in (1..SYS_N).rev() {
            qq[j] = self.U[j][j];
            for i in 0..SYS_N {
                let mut sum = self.F[i][j];
                for k in 0..j {
                    sum += self.F[i][k] * self.U[k][j];
                }
                w[i][j] = sum;
            }
        }
        qq[0] = self.U[0][0];
        // wの右NxR要素を初期化
        for i in 0..SYS_N {
            for j in 0..SYS_M {
                w[i][j + SYS_N] = self.G[i][j];
            }
            w[i][0] = self.F[i][0];
        }
        // --- ここまででw, qq, self.xを計算

        for j in (1..SYS_N).rev() {
            let mut sum = 0.0;
            for k in 0..SYS_N {
                z[k] = w[j][k] * qq[k];
                sum += z[k] * w[j][k];
            }
            for k in SYS_N..(SYS_N + SYS_M) {
                sum += w[j][k] * w[j][k];
            }
            self.U[j][j] = sum;
            let u_recip = self.U[j][j].recip();
            for i in 0..j {
                sum = 0.0;
                for k in 0..SYS_N {
                    sum += w[i][k] * z[k];
                }
                for k in SYS_N..(SYS_N + SYS_M) {
                    sum += w[i][k] * w[j][k];
                }

                sum *= u_recip;
                for k in 0..(SYS_N + SYS_M) {
                    w[i][k] -= sum * w[j][k];
                }
                self.U[i][j] = sum;
            }
        }
        let mut sum = 0.0;
        for k in 0..SYS_N {
            sum += qq[k] * (w[0][k] * w[0][k]);  // qqには更新前のUの対角要素が入っている
        }
        for k in SYS_N..(SYS_N + SYS_M) {
            sum += w[0][k] * w[0][k];
        }
        self.U[0][0] = sum;
    }

    /// フィルタリングステップ
    pub fn filtering(&mut self, y: &VectorP<f64>) {
        // Working array
        let mut ff: VectorN<f64> = unsafe {MaybeUninit::uninit().assume_init()};  // U^T H^T
        let mut gg: VectorN<f64> = unsafe {MaybeUninit::uninit().assume_init()};  // D U^T H^T

        self.calc_jacobian_h();
        let yhat = calc_h(self.x);

        // 出力の数だけループ
        for l in 0..SYS_P {
            // y_diff := y - h(x)
            let mut y_diff = y[l] - yhat[l];
            
            for j in (1..SYS_N).rev() {
                ff[j] = self.H[l][j];
                for k in 0..j {
                    ff[j] += self.U[k][j] * self.H[l][k];
                }
                gg[j] = self.U[j][j] * ff[j];
            }
            ff[0] = self.H[l][0];
            gg[0] = self.U[0][0] * ff[0];
            // --- ここまででy_diff, ff, ggを計算

            let mut alpha = self.R[l] + ff[0] * gg[0];  // 式 8.46
            let mut gamma = alpha.recip();
            self.U[0][0] = self.R[l] * gamma * self.U[0][0];  // 式 8.46
            for j in 1..SYS_N {
                let mut beta = alpha;
                alpha += ff[j] * gg[j];  // 式 8.47
                let lambda = ff[j] * gamma;  // 式　8.49
                gamma = alpha.recip();
                self.U[j][j] = beta * self.U[j][j] * gamma;  // 式 8.48
                for i in 0..j {
                    beta = self.U[i][j];
                    self.U[i][j] -= lambda * gg[i];  // 式 8.50
                    gg[i] +=  beta * gg[j];  // 式 8.51
                }
            }

            y_diff *= gamma;
            for j in 0..SYS_N {
                self.x[j] += gg[j] * y_diff;
            }
        }

        // 最後に正規化
        let norm_recip = quat::norm((self.x[0], [self.x[1], self.x[2], self.x[3]])).recip();
        for i in 0..4 {
            self.x[i] *= norm_recip;
        }
    }

    fn calc_jacobian_f(&mut self, gyr: Vector3<f64>) {
        // df1/dx1
        let [r1, r2, r3] = quat::scale_vec(FRAC_DT_2, quat::sub_vec(gyr, [self.x[4], self.x[5], self.x[6]]));
        self.F[0][0] = 1.0;
        self.F[0][1] = -r1;
        self.F[0][2] = -r2;
        self.F[0][3] = -r3;
        self.F[1][0] =  r1;
        self.F[1][1] = 1.0;
        self.F[1][2] =  r3;
        self.F[1][3] = -r2;
        self.F[2][0] =  r2;
        self.F[2][1] = -r3;
        self.F[2][2] = 1.0;
        self.F[2][3] =  r1;
        self.F[3][0] =  r3;
        self.F[3][1] =  r2;
        self.F[3][2] = -r1;
        self.F[3][3] = 1.0;

        // df1/dx2
        let (q0, [q1, q2, q3]) = quat::scale(FRAC_DT_2, (self.x[0], [self.x[1], self.x[2], self.x[3]]));
        self.F[0][4] =  q1;
        self.F[0][5] =  q2;
        self.F[0][6] =  q3;
        self.F[1][4] = -q0;
        self.F[1][5] =  q3;
        self.F[1][6] = -q2;
        self.F[2][4] = -q3;
        self.F[2][5] = -q0;
        self.F[2][6] =  q1;
        self.F[3][4] =  q2;
        self.F[3][5] = -q1;
        self.F[3][6] = -q0;

        // df2/dx2
        for i in 0..3 {
            self.F[i+4][i+4] = 1.0 + DT * BETA_DIAG[i];
        }
    }

    fn calc_jacobian_h(&mut self) {
        let q0 = self.x[0];
        let q1 = self.x[1];
        let q2 = self.x[2];
        let q3 = self.x[3];

        // dh1/dx1
        let r1 = ACC_R[0];
        let r2 = ACC_R[1];
        let r3 = ACC_R[2];
        self.H[0][0] =  r2*q3 - r3*q2;
        self.H[0][1] =  r3*q3 + r2*q2;
        self.H[0][2] = -q0*r3 - 2.0 * r1*q2 + r2*q1;
        self.H[0][3] =  q0*r2 + r3*q1 - 2.0 * r1*q3;
        self.H[1][0] =  r3*q1 - r1*q3;
        self.H[1][1] =  q0*r3 + q2*r1 - 2.0 * q1*r2;
        self.H[1][2] =  r1*q1 + r3*q3;
        self.H[1][3] = -q0*r1 - 2.0 * r2*q3 + r3*q2;
        self.H[2][0] =  r1*q2 - r2*q1;
        self.H[2][1] = -q0*r2 - 2.0 * q1*r3 + r1*q3;
        self.H[2][2] =  q0*r1 + r2*q3 - 2.0 * r3*q2;
        self.H[2][3] =  r2*q2 + r1*q1;
        for i in 0..3 {
            for j in 0..4 {
                self.H[i][j] *= 2.0;
            }
        }

        // dh2/dx1
        let r1 = MAG_R[0];
        let r2 = MAG_R[1];
        let r3 = MAG_R[2];
        self.H[3][0] = r2 * q3 - r3 * q2;
        self.H[3][1] = r3 * q3 + r2 * q2;
        self.H[3][2] = -q0 * r3 - 2.0*r1*q2 + r2*q1;
        self.H[3][3] = q0 * r2 + r3 * q1 - 2.0 * r1 * q3;
        self.H[4][0] = r3 * q1 - r1 * q3;
        self.H[4][1] = q0 * r3 + q2 * r1 - 2.0 * q1 * r2;
        self.H[4][2] = r1 * q1 + r3 * q3;
        self.H[4][3] = -q0 * r1 - 2.0 * r2 * q3 + r3 * q2;
        self.H[5][0] = r1 * q2 - r2 * q1;
        self.H[5][1] = -q0 * r2 - 2.0 * q1 * r3 + r1 * q3;
        self.H[5][2] = q0 * r1 + r2 * q3 - 2.0 * r3 * q2;
        self.H[5][3] = r2 * q2 + r1 * q1;
        for i in 0..3 {
            for j in 0..4 {
                self.H[i+3][j] *= 2.0;
            }
        }
    }
}


/// U-D分解（P = U * D * U^T）
/// 
/// * Pをn×n非負正定値対称行列とする．
/// * Uは対角成分を1とするn×n上三角行列．
/// * Dはn×n対角行列．
/// 
/// 返り値は，対角成分をDとし，それ以外の要素をUとした上三角行列．
fn ud_decomp(mut p: MatrixNxN<f64>) -> MatrixNxN<f64> {
    let mut ud: MatrixNxN<f64> = unsafe {MaybeUninit::uninit().assume_init()};

    for k in (1..SYS_N).rev() {  // n-1, n-2, ..., 1
        ud[k][k] = p[k][k];
        let ud_recip = ud[k][k].recip();
        for j in 0..k {
            ud[j][k] = p[j][k] * ud_recip;
            ud[k][j] = 0.0;  // 対角を除いた下三角成分を0埋め

            let tmp = ud[j][k] * ud[k][k];
            for i in 0..=j {  // 両側閉区間
                p[i][j] -= ud[i][k] * tmp;
            }
        }
    }
    ud[0][0] = p[0][0];  // pを書き換えてるから，d[0]の代入はこの位置じゃないとダメ

    ud
}

/// コレスキー分解（P = U * U^T）
/// 
/// * Pをn×n非負正定値対称行列とする．
/// * Uは対角要素が非負の値をとるn×n上三角行列．
fn cholesky_decomp(mut p: MatrixMxM<f64>) -> MatrixMxM<f64> {
    let mut u: MatrixMxM<f64> = unsafe {MaybeUninit::uninit().assume_init()};

    for k in (1..SYS_M).rev() {
        u[k][k] = p[k][k].sqrt();
        let u_recip = u[k][k].recip();
        for j in 0..k {
            u[j][k] = p[j][k] * u_recip;
            u[k][j] = 0.0;  // 対角を除いた下三角成分を0埋め
            for i in 0..=j {
                p[i][j] -= u[i][k] * u[j][k];
            }
        }
    }
    u[0][0] = p[0][0].sqrt();

    u
}

pub fn calc_f(x: VectorN<f64>, gyr: Vector3<f64>) -> VectorN<f64> {
    let mut q = (x[0], [x[1], x[2], x[3]]);
    
    // 積分（q = q + 0.5*Δt*q*ω）
    let gyr = [gyr[0] - x[4], gyr[1] - x[5], gyr[2] - x[6]];
    let tmp0 = quat::scale_vec(q.0, gyr);
    let dot = quat::dot_vec(q.1, gyr);
    let cross = quat::cross_vec(q.1, gyr);
    let tmp1 = (-dot, quat::add_vec(tmp0, cross));
    q = quat::scale_add(FRAC_DT_2, tmp1, q);

    q = quat::normalize(q);

    [
        q.0,
        q.1[0],
        q.1[1],
        q.1[2],
        x[4] + DT * BETA_DIAG[0] * x[4],
        x[5] + DT * BETA_DIAG[1] * x[5],
        x[6] + DT * BETA_DIAG[2] * x[6],
    ]
}

pub fn calc_h(x: VectorN<f64>) -> VectorP<f64> {
    let q = (x[0], [x[1], x[2], x[3]]);
    let h1 = quat::frame_rotation(q, ACC_R);
    let h2 = quat::frame_rotation(q, MAG_R);
    let mut yhat: VectorP<f64> = unsafe {MaybeUninit::uninit().assume_init()};
    for i in 0..3 {
        yhat[i] = h1[i];
        yhat[i + 3] = h2[i];
    }
    yhat
}