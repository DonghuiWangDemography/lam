# LAM 中文使用说明

这是 LAM 0.1.0，论文主模型的通用 Stata 实现。不同调查可以有不同数量的回答类别和不同阈值。同一个时期的调查共享潜在态度均值，利用调查之间重叠的时期把结果放到同一尺度上。

最基本的调用只需要一行。

```stata
lam response, survey(source) time(year)
```

其中 `response` 是数值型有序回答，`source` 是调查序列标识，`year` 是数值型时期变量。变量名可以换成自己的数据中的名称。

## 先运行已有数据

下载或克隆 [LAM 仓库](https://github.com/DonghuiWangDemography/lam)，在 Stata 中把工作目录设为其中的 `stata-lam` 文件夹。先保存尚未保存的数据，再运行下面的代码。

```stata
adopath ++ "."
capture mkdir results
local china "../lam-main/lam-china.dta"
capture confirm file "`china'"
if _rc local china "../lam-china.dta"
use "`china'", clear
lam k, survey(q) time(s)
lam_plot
graph export "results/china_trend.png", replace width(1800)
```

也可以直接运行完整示例。示例会自动建立结果文件夹、保存估计结果、导出趋势表和阈值表，并生成 PNG 图片。

```stata
do examples/china.do
```

以上代码和完整示例都会先检查 `../lam-main/lam-china.dta`，再检查 `../lam-china.dta`，因此支持数据在子目录或仓库根目录的两种结构。示例读取原数据，不会修改原文件。

第一次使用时，可以按 [中文测试教程](testing-zh.md) 操作，也有对应的 [English testing guide](testing-en.md)。

如果只想先检查命令是否能运行，可以使用不依赖外部数据的模拟示例。

```stata
do examples/synthetic.do
```

## 安装为常用命令

可以直接从 GitHub 安装。

```stata
net install lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
help lam
```

上面的命令安装估计程序和帮助文件。如果还需要把说明文档和示例复制到当前工作目录，可以另外运行下面的命令。

```stata
net get lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
do synthetic.do
```

`net get` 会把这些附加文件直接复制到当前目录，不保留 `examples` 和 `docs` 子目录。通过这种方式取得的模拟示例使用 `do synthetic.do` 运行。

如果已经下载整个仓库，也可以从本地安装。将示例路径替换为本机仓库目录。

```stata
net install lam, from("/path/to/lam/stata-lam") replace
```

以后打开 Stata 就可以直接调用 `lam`。这个版本使用 ado 和 Mata 源代码，不需要安装第三方包，也不需要针对不同电脑编译 plugin。源码以 Stata 16 及以上版本为目标，实际运行验证使用 macOS 上的 StataNow 19.5。更早的 Stata 版本和其他操作系统尚未进行运行测试。

## 数据需要满足的条件

- 所有调查中的回答方向一致。数值越大应当表示同一个方向的态度。反向题需要先反向编码。
- 每个调查至少有两个实际观测到的回答类别。编码不必连续，但必须具有明确的数值顺序。
- 调查之间必须通过重叠时期连接起来。A 与 B 重叠，B 与 C 重叠，就可以连接 A、B、C，不要求三个调查在同一年同时存在。
- 同一调查的测量方式和阈值在时间上保持稳定。重叠时期的调查应当代表同一个目标总体并测量同一概念。

程序可以检查数据结构和重叠关系，但不能替代研究者对测量一致性和样本代表性的判断。如果某个名义回答类别在一个调查的整个估计样本中都没有出现，程序只使用实际观测到的类别，不能单独识别这个空类别的边界。

## 查看估计结果

```stata
matrix list e(trend)
matrix list e(cutpoints)
ereturn list
```

`e(trend)` 的五列分别是时期、均值估计、标准误、置信区间下界和上界。默认把最早的时期均值设为零。这个时期的标准误和区间也为零，这是尺度归一化的结果。

如果要换一个基准时期，可以明确指定它。

```stata
lam k, survey(q) time(s) baseperiod(1989) nolog
```

`e(cutpoints)` 保存调查编号、阈值下方的原始类别编码、阈值估计和标准误。`e(survey1)`、`e(survey2)` 等保存内部调查编号对应的原始调查值。`e(categories)` 保存每个调查中的类别顺序。

程序同时提供标准的 `e(b)` 和 `e(V)`，可以保存估计结果或计算线性组合。均值系数的名称使用排序后的时期序号，而不是直接使用年份。比如默认基准时期为第一个时期时，`mu:t2` 表示第二个观测时期相对于基准时期的均值。

## 预测与作图

```stata
predict double attitude, xb
predict double attitude_se, stdp
predict double p_category2, pr outcome(2)
lam_plot, title("Latent attitude trend")
graph export "results/trend.png", replace width(1800)
```

`attitude` 是对应时期的总体潜在态度均值，不是每个受访者个人态度的后验估计。`p_category2` 是该记录所属调查和时期下，原始回答代码为 2 的预测概率。预测所需的时期、调查或类别没有出现在拟合模型中时，对应预测为缺失值。

## 使用频数表

如果已经有调查、时期、回答类别及该类别人数，就不需要再把数据扩展成每人一行。

```stata
lam category [fw=count], survey(source) time(year)
```

`count` 必须是非负整数人数。频数输入与逐人展开数据使用相同的似然，得到相同的估计和观测信息矩阵标准误。频数权重不是复杂抽样权重，不能直接把调查权重当作人数填入。

默认置信区间来自观测信息矩阵。对个体记录使用 Stata 的 `bootstrap` 前缀时，需要保留相应抽样结构，并确认各次重复的类别支持和系数维度一致。直接对汇总频数表的行进行重抽样，不等于对受访者进行重抽样。

使用 `bootstrap` 后，`predict, stdp` 和 `lam_plot` 会读取更新后的 `e(V)`。此时图中的区间使用 bootstrap 标准误和正态近似，不是百分位数区间。`e(trend)` 和 `e(cutpoints)` 是原始拟合时保存的表，外部前缀不会自动更新它们。可以在运行 `lam_plot` 后，用 `matrix trend = r(trend)` 取得当前图中使用的结果表。

程序也会检查是否存在有限的最大似然估计。如果某些时期只有极端类别等响应结构导致分离，程序会明确报错，避免把优化器停止迭代当作有效估计。

## 当前版本范围

当前版本实现论文主模型，即调查特定、时间不变的阈值和自由的时期均值。它不包含协变量、复杂抽样权重、聚类稳健标准误或随时间变化的阈值，也不会插补未观测时期。

维护者为 [Donghui Wang](https://github.com/DonghuiWangDemography)。可复现的问题可以提交到 [仓库 Issues](https://github.com/DonghuiWangDemography/lam/issues)。0.1.0 版本没有附带软件许可证。
