# 结果文件

- `steps_000_200_eight_methods.csv`：201 行宽表，每种方法的能量和标准误各一列；包含 θ、理想模拟器值和未校正的 N=2 保留率。
- `steps_000_200_long.csv`：201 × 8 = 1608 行，适合绘图和分组分析。
- `0-200逐步结果.md`：201 步八方法的可阅读表格。
- `method_summary.csv`：每种方法的初始值、末值、最低值、对应步数和 θ。
- `offline_methods_comparison.csv`：完整中间量，包括 IBU 10/50/100 轮、投影分子/分母和有效 bootstrap 数。
- `offline_methods_comparison.png/pdf`：两套端点校准下的曲线和采样误差条；每 10 步绘一根误差条。
- `integrity_check.json` 与 `offline_methods_validation.json`：数据完整性、参数一致性及约束优化的数值检查。

最低值未经独立重复测量；采样标准误不包含门误差、设备漂移、读出模型及自适应选择误差。端点校准曲线用于敏感性分析。
