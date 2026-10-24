# SeawaterTracker：单独整理的论文核心代码

本文件夹对应你提供的 **SeawaterTracker_portable.rar**。它整理的是测序数据来源分析流程，核心文件可以直接查看和运行。

## 文件与功能对应

| 文件名 | 完成的功能 |
| --- | --- |
| `source_tracking.py` | 来源比例采样、独立重复、均值、标准差、相对标准差 |
| `primer_trimming.py` | 去除双端测序引物，调用 Cutadapt |
| `data_preparation.py` | ASV 表读取、特征对齐、按来源类别聚合 |
| `dada2_denoise.R` | DADA2 过滤、去噪、拼接、去嵌合体、输出 ASV 表 |
| `run_source_tracking.py` | 读取数据与参数，运行来源比例分析并保存结果 |

函数和变量沿用已规范的英文命名，参数集中在 `configs/source_tracking.json`。本包可独立使用。

## 最快运行方式

在本文件夹中打开 PowerShell，使用 Python 3.11：

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run_source_tracking.py --sources examples/source_counts.tsv --sinks examples/sink_counts.tsv --config examples/demo_parameters.json --output results/demo
```

这个小示例使用合成数据，便于检查安装和运行。结果保存在 `results/demo/`：

- `mean_proportions.tsv`：来源比例，0–1。
- `std_proportions.tsv`：独立重复结果的标准差，0–1。
- `relative_std_percent.tsv`：相对标准差，百分数。
- `run_metadata.json`：本次参数与输入文件校验值。

附带的研究数据在 `data/sequencing/`：21 类来源、3 个待测样本、78,705 个对齐后的特征。其运行命令和 Python 函数调用示例见 [英文说明](README.md)。

## 已完成的核对

已重新核对压缩包 SHA-256，并将本包的核心代码与先前验证版本对应起来。测试覆盖固定种子的采样结果、独立重复、均值/标准差/RSD、重采样、数据对齐及命令行运行。R 脚本的变量改名已通过解析结构对比。

本次保持原计算规则，软件中已有的算法和参数差异仍记录在 [方法说明](docs/methods.md)。完整原始测序流程和论文全部图表尚未重新运行；恢复出的自定义采样器也未被认定为与官方 SourceTracker 实现完全相同。

具体记录见 `validation/archive_provenance.json`、`validation/numerical_comparison.json` 和 `validation/release_checks.json`。
