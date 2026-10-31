# 离线复现与核验

1. 在 Python 环境中安装 `requirements.lock.txt` 中的依赖。该文件记录本次离线分析的软件版本，不是完整的全平台环境锁定文件。
2. 运行 `python code/offline_methods_comparison.py`：从 `data/experiment_state.json` 的原始 counts 和开始/结束对照，离线重算三种 REM + N=2 曲线，包含 200 次 bootstrap；同时更新 `results/offline_methods_comparison.csv`、图和验证记录。该脚本不访问量子云。
3. 运行 `python scripts/derive_results.py`：检查每一步五基 counts、布局、shots、参数及 Raw/N2/模拟器数值，与保存的状态核对后重新生成八方法宽表、长表、汇总表和逐步 Markdown。
4. 比较 `results/integrity_check.json` 和 `results/offline_methods_validation.json`。校验只是确认数据与实现一致，不构成对纠错模型真实准确性的证明。

IBU 固定均匀先验，并输出 10、50、100 轮；主要表格和图使用 50 轮。约束最小二乘对每个测量基独立拟合。采样误差处理、端点校准假设及已知局限详见 `protocols/experiment.md`。

原始线路快照为去除个人身份和本机凭据路径后的版本，源文件原始哈希记录于实验状态的 manifest；该快照不是已配置认证的云端可运行脚本。实验记录中的任务没有逐任务时间戳。原始状态记录零个失败尝试，不排除未记录的云端失败；若有未记录失败，其实际 shots 未知。

`SHA256SUMS.txt` 是整理时的文件快照；重算结果后如内容改变，请重新生成文件哈希再核对。实验没有被重新提交。
