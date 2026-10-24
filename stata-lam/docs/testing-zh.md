# 在 Stata 中测试 LAM

建议先运行中国态度数据示例，确认结果和图能出来，再运行完整测试。先下载或克隆 [LAM 仓库](https://github.com/DonghuiWangDemography/lam)，把下文的 `/path/to/lam` 替换为本机仓库根目录。该目录应当包含 `stata-lam` 文件夹和原始数据。示例和测试文件会执行 `clear all`，因此请先保存 Stata 内存中尚未保存的数据。

在 Stata 中打开 Do-file Editor，把下面一个代码块完整复制进去，再点击 Do 运行。也可以逐行粘贴到 Command 窗口。不要把 Markdown 的代码围栏一起复制进去。

## 第一步，运行原来的数据

```stata
cd "/path/to/lam/stata-lam"
do examples/china.do
```

这个示例自动加载原始数据、估计模型、计算预测、绘图并保存结果。原始数据文件保持不变。它支持原始数据位于项目的 `lam-main/lam-china.dta`，或在 GitHub 仓库根目录直接存放为 `lam-china.dta`。

运行结束后应该看到一张 Americans' attitudes toward China 趋势图。结果位于 `stata-lam/results` 文件夹。

| 文件 | 内容 |
| --- | --- |
| `china_trend.png` | 潜在态度趋势及 95% 置信区间 |
| `china_trend.csv` | 年份、估计值、标准误和置信区间 |
| `china_trend.dta` | 同一趋势表的 Stata 数据文件 |
| `china_cutpoints.csv` | 各调查的阈值估计与标准误 |
| `china_lam.ster` | 保存的 Stata 估计结果 |

接着运行下面的代码检查结果。`assert` 成功时不打印消息，最后的文字表示这几项检查都通过。

```stata
assert e(converged) == 1
assert e(N) == 108036
assert e(n_surveys) == 13
assert e(n_periods) == 40
assert e(k) == 89
assert abs(e(ll) - (-155574.5182877433)) < 0.0001
display %18.4f e(ll)
display "CHINA QUICK CHECK PASSED"
```

对数似然应显示约 `-155574.5183`。系数表标题中可能使用科学计数法，这是显示精度问题，以上 `display` 命令会显示更多位数。小的跨平台浮点差异是允许的。

```stata
matrix list e(trend), format(%9.4f)
matrix list e(cutpoints), format(%9.4f)
ereturn list
```

`e(trend)` 的五列是时期、估计值、标准误、置信区间下界和上界。1974 年固定为零，用来确定尺度。它的标准误和区间也为零，这是归一化，不是对总体态度没有不确定性。其他年份的估计是相对于 1974 年的潜在态度均值。

`e(b)` 保存所有自由参数，`e(V)` 保存协方差矩阵。当前示例的置信区间使用观测信息矩阵，不是 bootstrap 区间。

## 第二步，自己调用一次命令

如果希望直接调用估计器，可以运行下面这段代码。它支持原始数据在子目录或仓库根目录的两种结构。

```stata
cd "/path/to/lam/stata-lam"
adopath ++ "."
local china "../lam-main/lam-china.dta"
capture confirm file "`china'"
if _rc local china "../lam-china.dta"
use "`china'", clear
lam k, survey(q) time(s) nolog
predict double fitted_mean, xb
predict double fitted_mean_se, stdp
predict double probability_category2, pr outcome(2)
lam_plot
```

三项输入分别是回答类别 `k`、调查序列 `q` 和年份 `s`。以后换成自己的数据时，替换这三个变量名即可。`fitted_mean` 是该年份的总体潜在态度均值，不是每个受访者的个人态度分数。`probability_category2` 是该调查、该年份回答代码为 2 的预测概率。

重新加载数据后再运行整段代码，可以避免重复创建同名预测变量。仅重新运行 `predict` 时，如果变量已经存在，Stata 会提示变量已定义，这不表示模型失败。

## 第三步，运行不依赖原数据的模拟示例

```stata
cd "/path/to/lam/stata-lam"
do examples/synthetic.do
```

示例自动生成 14,000 条记录、3 个调查和 10 个时期，然后拟合模型。预计对数似然约为 `-15409.43571`。示例内还会检查每条记录的类别概率之和是否为 1，以及频数输入和个体记录输入是否得到相同的对数似然。

这个示例结束时没有独立的总通过标记。应完整运行到 `end of do-file`，且没有中途停止的 `assertion is false`。主要输出是 `results/synthetic_trend.png` 和 `results/synthetic_trend.csv`。CSV 还包含真实均值 `truth` 和估计误差 `error`，便于比较。这里重新生成的数据不是逐字复现旧模拟文件的随机抽样。

## 第四步，运行完整自动测试

完整测试需要把工作目录切回仓库根目录，也就是 `stata-lam` 的上一级。这和运行示例时的工作目录不同。请把整个代码块放在同一个 do-file 中运行。

```stata
cd "/path/to/lam"
do stata-lam/tests/test_lam.do
do stata-lam/tests/test_china.do
do stata-lam/tests/test_derivatives.do
do stata-lam/tests/test_separation.do
do stata-lam/tests/test_postestimation.do
display "ALL FIVE LAM TEST SUITES COMPLETED"
```

每个测试文件都应到达自己的最终通过标记。

| 测试文件 | 最终通过标记 |
| --- | --- |
| `test_lam.do` | `ALL CORE TESTS PASSED` |
| `test_china.do` | `CHINA REPLICATION PASSED` |
| `test_derivatives.do` | `DERIVATIVE AND TAIL TESTS PASSED` |
| `test_separation.do` | `SEPARATION AND PERIOD RECODING TESTS PASSED` |
| `test_postestimation.do` | `ALL POSTESTIMATION TESTS PASSED` |

测试会故意输入不连通调查、非法权重、分离数据或非法绘图选项，因此中途出现错误文字可能是预期行为。程序会捕获这些错误并检查错误码。判断是否成功应看五个最终标记是否都出现，而不是要求输出中完全没有红色文字。如果在某个标记前停止，先查看该测试日志最后几行。

日志保存在 `stata-lam/tests` 中，每个文件对应同名 `.log`。中国数据测试会独立重新估计原始似然，并比较系数、完整协方差矩阵及对数似然。其他测试覆盖官方 `oprobit` 对照、导数、频数等价性、重叠关系、识别、预测和作图。

如果更习惯 Mac 终端，也可以使用已经提供的运行脚本。

```sh
cd "/path/to/lam"
sh stata-lam/tests/run_tests.sh
```

最终应显示 `All available LAM checks passed.`。如果终端提示跳过中国数据测试，说明没有找到原始数据，这时不能把该次运行视为中国数据复现也已经通过。Stata 批处理的进程退出码本身不足以判断全部断言是否运行，脚本会额外检查通过标记。

## 安装为日常命令

确认示例成功后，可以从本地文件夹安装。

```stata
net install lam, from("/path/to/lam/stata-lam") replace
which lam
help lam
```

`which lam` 显示 Stata 实际使用的文件位置。如果修改代码后仍在使用之前加载的版本，先保存数据并重新打开 Stata，再执行测试。

也可以直接从 GitHub 安装 0.1.0 版本。程序位于仓库的 `stata-lam` 目录，原来的论文代码和数据保留在仓库根目录。

```stata
net install lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
```

`net install` 安装程序和帮助文件。若需要附加说明及独立模拟示例，可在希望保存文件的目录运行下面的命令。

```stata
net get lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
do synthetic.do
```

`net get` 会把示例文件直接复制到当前目录，因此此时使用 `do synthetic.do`。中国态度原始数据不属于安装包，完整复现测试应在包含原始数据的项目或仓库副本中运行。
