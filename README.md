# 姿勢推定フィルタの比較シミュレーション

宇宙電波実験室の記事内で掲載したシミュレーションを行うプログラムです。
実装する際の参考にしてください。

該当記事：[角速度の形でフィードバックを行う姿勢推定フィルタ（宇宙電波実験室）](https://space-denpa.jp/2022/03/01/omega-feedback-filter/)

## 内容

各ディレクトリにRustで記述したシミュレーション用プログラムが入っています。

* /omega_feedback：記事中で説明したフィルタ
* /complementary：相補フィルタ
* /ekf：拡張カルマンフィルタ

## 使い方

シミュレーションを実行する際は、各ディレクトリに移動して以下のコマンドを実行してください。

```
$ cargo run && python3 data_plot.py
```
