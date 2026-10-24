# LAM 的 Stata 实现

LAM 整合回答量表不同的调查序列，估计同一潜在态度随时间变化的趋势。本仓库保留原论文的复现材料，并提供可直接调用的通用 Stata 命令。

[English README](README.md) · [完整中文说明](stata-lam/docs/quickstart-zh.md) · [Stata 测试教程](stata-lam/docs/testing-zh.md)

## 安装

在 Stata 中运行下面的命令。

```stata
net install lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
help lam
```

程序使用 ado 和 Mata 源代码，无需第三方包。源码以 Stata 16 及以上版本为目标，实际运行验证使用 StataNow 19.5 的 macOS 版本。其他版本和操作系统尚未逐一运行验证。

## 运行完整示例

安装后，在希望保存结果的工作目录中运行下面的代码。示例会载入新数据，运行前请先保存 Stata 中尚未保存的数据。

```stata
net get lam, from("https://raw.githubusercontent.com/DonghuiWangDemography/lam/main/stata-lam") replace
do synthetic.do
```

它自动生成三个相互连接的调查序列、估计模型、检查预测概率和频数输入，并在 `results` 中生成 PNG 趋势图及结果表。

要运行原论文中的中国态度数据，请下载或克隆整个仓库，再把下面的路径换成本机仓库位置。

```stata
cd "/path/to/lam/stata-lam"
do examples/china.do
matrix list e(trend), format(%9.4f)
```

示例读取仓库根目录的 `lam-china.dta`，不修改原始文件。结果应包括 108,036 个观测、13 个调查和 40 个时期，对数似然约为 -155574.5183。趋势图和结果表保存在 `stata-lam/results`。

## 用自己的数据估计

```stata
lam response, survey(source) time(year) nolog
matrix list e(trend)
matrix list e(cutpoints)
predict double attitude, xb
predict double probability, pr outcome(2)
lam_plot
graph export "trend.png", replace width(1800)
```

`response` 是数值型有序回答，`source` 是调查序列标识，`year` 是数值型时期。调查标识可以是数值或字符串。响应代码和时期不必连续，但所有调查的回答方向需要一致。

如果已经汇总为人数，可以直接使用频数。

```stata
lam response [fw=count], survey(source) time(year)
```

`count` 是非负整数人数。频数输入与逐人记录使用相同的似然。它不是复杂抽样权重，不能直接填入调查的 sampling weight。

## 方法和结果解释

各调查有各自的有序阈值，阈值在调查内保持时间不变。同一时期的调查共享潜在态度均值。模型将潜在误差方差固定为一，并默认将最早观测时期的均值固定为零。可以用 `baseperiod(#)` 更换基准时期。

调查必须通过重叠时期连成一个网络，间接连接也可以。调查还需要测量同一概念，并在重叠时期代表同一目标总体。这些实质假设需要研究者判断。程序会检查网络连接，并拒绝不存在有限最大似然估计的分离响应结构。

`predict, xb` 返回某时期的总体潜在态度均值，不是受访者个人态度分数。基准时期的零标准误来自尺度归一化。默认标准误和置信区间使用观测信息矩阵，bootstrap 属于另外的推断步骤。直接对汇总频数表的行重抽样，不等于对受访者重抽样。

当前版本实现论文固定阈值主模型，尚未包含协变量、复杂抽样权重、聚类稳健标准误、logit 链接、时间变化阈值或未观测时期插补。

## 测试和计算改进

新版与原程序的 89 个自由参数及完整协方差矩阵进行了比较。系数最大绝对差异约为 0.000008832，协方差元素最大绝对差异约为 0.0000003143。

计算时将中国数据的 108,036 条记录汇总为 483 个有正频数的调查、时期和回答类别组合，再直接计算解析梯度及 Hessian。频数保留每条记录对似然的贡献，从而减少重复计算。

在仓库根目录运行以下五组测试。

```stata
cd "/path/to/lam"
do stata-lam/tests/test_lam.do
do stata-lam/tests/test_china.do
do stata-lam/tests/test_derivatives.do
do stata-lam/tests/test_separation.do
do stata-lam/tests/test_postestimation.do
```

部分测试会故意触发错误，再核对错误码。判断通过与否，应检查每组最后的 `PASSED` 标记。详细步骤见 [中文测试教程](stata-lam/docs/testing-zh.md)。

原来的 `lam-ml.do`、`lam-sim.do` 和数据保留在仓库根目录。旧脚本中有机器特定路径，日常使用请运行 `stata-lam/examples` 中的示例。安装包本身不包含原始实证数据。

## 论文与反馈

Wang, Donghui, Yu Xie, and Junming Huang. 2023. *Trend Analysis with Pooled Data from Different Survey Series. The Latent Attitude Method.* Sociological Methodology. [DOI 10.1177/00811750231193641](https://doi.org/10.1177/00811750231193641).

软件问题可以提交到 [GitHub Issues](https://github.com/DonghuiWangDemography/lam/issues)，附上 Stata 版本、运行命令、错误信息及可复现的小例子。方法与旧代码中的差异见 [方法审查](stata-lam/docs/method-audit.md)，版本说明见 [release notes](stata-lam/docs/release-notes.md)。此版本未附加明确的软件许可证。
