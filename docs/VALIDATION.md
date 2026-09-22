**验证记录：Algorithm 1–6**

环境：SageMath 10.7，Python 3.12。日期：2026-09-22。

48 项测试通过，覆盖采样概率质量、格基变换、Gaussian 支持界、局部同余、留数误差区间、满秩真子格、例外单位、紧凑表示和记录复核。端到端求解期间禁用 `pari_bnf`，随后才使用 Sage 的 S-unit 群接口作独立对照。

射线类群检查另覆盖 Q(i) 的模 (3)：仅根单位不能生成全部剩余类单位，加入 1+i 后才满射。这区分了普通类群生成与射线类群生成两个前提。0.2.0 安装包已构建并通过从 wheel 导入的检查。

**已完成的端到端实例**

| 数域 / 目标 | 设置 | 本次结果 |
| --- | --- | --- |
| Q(sqrt(2))，普通单位群 | smooth_bound=59，generation_radius=2，seed=2026 | 20 条关系、5,452 次采样；自由生成元 sqrt(2)-1，根单位 -1 |
| Q(sqrt(-5))，S 为 2 上方的素理想 | smooth_bound=97，generation_radius=2，seed=11 | 25 条关系、5,244 次采样；自由生成元 -2，根单位 -1 |
| Q(i)，S 为 3 上方的素理想 | smooth_bound=59，generation_radius=2，seed=7 | 自动测试核对自由生成元在 Sage 基下的指数为 ±1 |
| Q(i)，S 为 2 和 3 上方的素理想 | smooth_bound=59，generation_radius=2，seed=31；模显式设为 2 上方素理想 | 完成例外单位路径；生成元 -3i、-2187(i+1)，在 Sage 自由基中的矩阵行列式为 -1 |

前两个实例的精确输出已保存为 [real_units.json](../examples/results/real_units.json) 和 [class2_sunits.json](../examples/results/class2_sunits.json)。它们的完整生成结论在 GRH 下成立，默认参考随机游走参数没有混合时间保证。

实二次域实例在达到满秩后，索引区间依次仍包含约 10^13、52–56、6，继续收集后才缩小到 {1}。这检查了“满秩”与“完整生成”的区别。类数为 2 的实例中，生成元 -2 在相应非主素理想上的赋值是 2，吻合类群约束。

**复核保存结果**

```bash
PYTHONPATH=src sage -python examples/verify_sunits.py examples/results/real_units.json
PYTHONPATH=src sage -python examples/verify_sunits.py examples/results/class2_sunits.json
```

复核器重新求留数区间、工作生成集的格基、额外赋值的整数核和最终子格指数。它忽略 JSON 中保存的完成标记和体积声明。测试会把一个记录中的最终生成元平方，并保留满秩及原来的成功标记；复核器必须得到指数 2 并拒绝该记录。

**认证混合界的局部集成检查**

对 Q(sqrt(2)) 的一次 Algorithm 3 运行，显式谱校准得到 B=38、N=114，专用 Sample 使用 epsilon=1/2199023255552。22 次候选后得到精确关系，运行报告为 `verified_spectral_bound`，组合论成功概率条件也通过。该校准使用了额外类群／单位群预处理，不作为默认 Algorithm 6 的依赖。

**尚未宣称的范围**

没有完成通用渐近位复杂度认证；生成半径候选值在部分例子中只能保持未验证。完整 S-unit 求解的端到端验证集中在二次数域，采样器另有混合三次域及全复四次域测试。上述关系次数与耗时不构成高次数或大判别式的性能保证。
