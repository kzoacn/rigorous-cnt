**rigorous-cnt：计算数论论文 Algorithm 1–6 的参考实现**

[English overview](README.en.md)

论文：[Markdownxiv 2609.00011v1](https://markdownxiv.github.io/abs/2609.00011v1/)（英文，Kaiyi Zhang）。[Markdown 源稿](paper/manuscript.md) · [发布记录](docs/PUBLICATION.md)

对应 [Rigorous methods for computational number theory，arXiv:2512.01588v2](https://arxiv.org/abs/2512.01588v2)。使用 SageMath 中的普通 Python，已在 SageMath 10.7 / Python 3.12 验证。

已实现理想采样、平滑关系、S-unit 关系、例外单位、BKP 后处理、完整生成判定和任意 S 的限制。输出采用紧凑幂乘积；小结果可以通过带高度界的模运算恢复，避免巨大中间展开。

**保证与前提**

- 理想关系、整数变换和支持集合使用精确运算。
- Algorithm 6 的完整生成结论在 GRH 下由独立的 zeta 留数区间和子格指数检查得到。
- 默认求解链不调用已有类群／单位群／S-unit 群求解器；Sage 对照只用于测试。
- 默认参考采样参数尚无通用混合时间保证，完整论文位复杂度也未认证。采样不足时会继续收集；预算耗尽会报错，不把部分生成集标成完整群。
- 可显式启用 Algorithm 2 的谱校准。它使用额外的类群、单位群预处理，需与核心求解成本分开理解。

详见 [Algorithm 3–6 数学说明](docs/ALGORITHMS_3_6.md)和 [Algorithm 1–2 数学说明](docs/MATHEMATICS.md)。

**直接运行**

在项目根目录执行：

```bash
# 实二次域的普通单位群
PYTHONPATH=src sage -python examples/demo_sunits.py --field real

# Q(i) 中 3 上方素理想对应的 S-unit 群
PYTHONPATH=src sage -python examples/demo_sunits.py --field imaginary --primes 3

# 显式覆盖例外单位分支
PYTHONPATH=src sage -python examples/demo_sunits.py --field imaginary --primes 2,3 --modulus-prime 2

# 保留 Algorithm 1–2 的采样演示
PYTHONPATH=src sage -python examples/demo.py --calibrate

# 全部测试
PYTHONPATH=src sage -python -m unittest discover -s tests -v

# 重新核验保存的结果，不信任记录中的完成标记
PYTHONPATH=src sage -python examples/verify_sunits.py examples/results/real_units.json
```

S-unit 示例会显示收集到的关系数、当前秩和子格指数区间。首次达到满秩后，可能仍需要更多关系。随机重试次数依数域、因子基和种子变化，例子通常需要数十秒到数分钟。

也可以在 Sage 环境中安装：

```bash
sage -pip install -e . --no-build-isolation
```

**完整 S-unit 求解接口**

```python
from sage.all import QQ, PolynomialRing, NumberField
from rigorous_cnt import NumberFieldContext, RandomBits, SUnitSettings, algorithm6

R = PolynomialRing(QQ, 't')
t = R.gen()
K = NumberField(t*t - 2, 'a')
ctx = NumberFieldContext(K)

settings = SUnitSettings(smooth_bound=59, generation_radius=2)
result = algorithm6(
    ctx, S=[], settings=settings, rng=RandomBits(seed=2026),
    progress=lambda state: print(state),
)

assert result.verify()
print(result.rank)
print(result.torsion_order, result.torsion_generator)
print(result.materialize_generators())
print(result.report)
```

`S=[]` 求普通单位群。一般 S 用素理想列表表示，例如 `S=K.primes_above(3)`。内部会扩张工作集合，认证其完整 S-unit 群后，通过额外赋值的饱和整数核收缩到请求的 S。

`result.generators` 是自由部分的紧凑生成元，`torsion_generator` 单独生成根单位部分。大结果可保持紧凑表示；`materialize_generators(max_bits=...)` 只有在认证高度预算允许时才恢复显式数域元素。`result.to_dict()` 导出精确的基础元素、整数指数、赋值与验证区间，演示脚本支持 `--output result.json`。

`verify_sunit_record(data)` 会忽略保存的完成标记，重新计算留数、工作格完整性、额外赋值的整数核及最终子格指数。项目保存了两个实际运行记录：[实二次域单位](examples/results/real_units.json)、[类数为 2 的数域 S-unit](examples/results/class2_sunits.json)。它们可以快速复核，无需重新执行随机关系搜索。

**单独运行 Algorithm 3–5**

```python
from rigorous_cnt import (
    RelationEngine, FactorBase, prime_ideals_up_to,
    algorithm3, algorithm4, algorithm5,
)

engine = RelationEngine(ctx, settings)
T = FactorBase(ctx, [
    P for P in prime_ideals_up_to(ctx, 59)
    if P not in engine.exceptional_primes
])

relation = algorithm3(engine, K.ideal(1), T, rng=RandomBits(2))
assert relation.verify()  # (alpha) = input_ideal * product(P**v)

observation = algorithm4(engine, T, rng=RandomBits(3))
print(observation.valuations)

# 若 engine.modulus 非平凡，可对其中每个素因子调用：
# exceptional = algorithm5(engine, q, T, rng=RandomBits(4))
```

Algorithm 3 的模由留数近似按照论文分支选择。`modulus_override` 用于显式指定平方自由模，报告会标记此选择；它也便于测试 Algorithm 5。

**采样参数与输出正确性分开报告**

`SUnitSettings` 记录因子基界、随机游走步数、生成半径候选值和资源预算。默认参数用于参考实验，报告中的 `sampling_mixing_guarantee` 通常为 `not_established`。结果仍须通过独立的完整生成判据才能返回。

生成半径和射线类群生成条件可以在完成后由实际生成集进一步检查，对应的布尔值会写入报告。即使这些检查通过，也不把它们当作整个程序的渐近复杂度证明。

若需要认证的采样混合界，可明确传入参数工厂：

```python
from rigorous_cnt import calibrate_walk

settings = SUnitSettings(
    require_mixing_bound=True,
    walk_factory=lambda context, ray, epsilon:
        calibrate_walk(context, ray, epsilon),
)
```

这个工厂沿用较重的谱校准预处理，确实会计算已有类群／单位群数据；默认 Algorithm 6 不采用它。它仅提供采样参数，返回的 S-unit 关系和完整生成证据仍由本项目计算。

**实现模块**

| 文件 | 内容 |
| --- | --- |
| `sampling.py`、`geometry.py` | Algorithm 1–2 与精确均匀采样 |
| `relations.py` | 专用 Sample、Algorithm 3–5 |
| `residue.py` | Belabas–Friedman 显式留数误差界 |
| `postprocess.py` | 两遍 BKP 与区间 Gram 行列式 |
| `sunits.py` | Algorithm 6、类群覆盖证据、完整性判定、任意 S |
| `compact.py` | 紧凑幂乘积、误差补偿及模恢复 |
| `factor_base.py` | 因子基和精确平滑性检测 |
| `calibration.py` | 可选的独立谱校准 |
| `verification.py` | 重新核验导出记录，拒绝满秩但不完整的结果 |

Algorithm 1–2 的详细 API 见 [采样使用说明](docs/ALGORITHMS_1_2.md)。目前采用极大整环和绝对数域表示，非极大整环扩展尚未实现；高次数 S-unit 求解尚无性能承诺。
