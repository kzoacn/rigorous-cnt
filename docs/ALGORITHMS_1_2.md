**rigorous-cnt：论文 Algorithm 1–2 的参考实现**

本项目实现 [Rigorous methods for computational number theory，arXiv:2512.01588v2](https://arxiv.org/abs/2512.01588v2) 中的受约束格点均匀采样与理想采样。运行环境为 SageMath 中的 Python，已在 SageMath 10.7 / Python 3.12 验证。

支持实、复及混合签名的绝对数域；整数与分式理想；有理平移；有限模；实嵌入符号限制；每个无限位上的有理实/复缩放；以及调用方提供的有限指数子群成员 oracle。

**可运行范围与保证**

- Algorithm 1 实现精确均匀采样，输出保留为符号形式 `x * alpha`。理想成员关系、局部同余、严格符号条件和几何边界用精确运算验证。
- Algorithm 2 实现随机素理想乘积、Lemma 2.22 离散高斯、有理扰动和 Algorithm 1 调用。返回精确的数域元素 `beta`，同时保留中间数据与运行报告。
- Algorithm 2 可用 `calibrate_walk(...)` 为小数域自动选取并检查 `B,N`。它显式计算类群、单位群和有限个字符的特征值界，结合 Appendix A.1 的 Gaussian 尾界验证混合误差。这个独立预处理支持完整 Arakelov 射线类群；一般子群仍需外部混合依据。
- **论文 Corollary 6.5 的高效通用参数界尚未具体化。** `calibrate_walk` 是较重的直接校准，不是对隐含常数的任意赋值。仅指定未经校准的 `B,N` 时，默认拒绝声明混合保证；显式关闭要求可做流程实验。另有 `WalkParameters.from_prime_bound(...)`，只从体积上界自动选择 `N`。
- **没有声称达到论文完整的位复杂度界。** 当前先用整数近似矩阵产生 LLL 候选，并验证原始精确格基的长度界；必要时回退到精确 LLL 和块 HKZ 枚举。枚举后端不是论文复杂度证明所用 Kannan 子程序的完整复现。

所有概率论断以独立均匀随机位为模型。固定种子支持实验重现。资源上限触发时抛出异常，不返回伪造样本；算法的无中断分布与截断执行需区分。数学细节、实现差异和未完成的理论保证见 [MATHEMATICS.md](MATHEMATICS.md)。

**运行**

在项目根目录执行，无需修改现有 Sage 环境：

```bash
PYTHONPATH=src sage -python examples/demo.py --calibrate
PYTHONPATH=src sage -python examples/demo.py --field imaginary --calibrate
PYTHONPATH=src sage -python examples/demo.py --field mixed --calibrate
PYTHONPATH=src sage -python -m unittest discover -s tests -v
```

也可在 Sage 环境中安装该包：

```bash
sage -pip install -e . --no-build-isolation
```

**Algorithm 1 的使用**

```python
from sage.all import QQ, PolynomialRing, NumberField
from rigorous_cnt import NumberFieldContext, RayConditions, RandomBits, Algorithm1

R = PolynomialRing(QQ, 't')
t = R.gen()
K = NumberField(t*t - 2, 'a')
ctx = NumberFieldContext(K)

# 实嵌入按生成元的像从小到大排序；这里要求第 0 个实嵌入为负。
ray = RayConditions(ctx, modulus=K.ideal(2), real_places=[0], tau=-1)
sampler = Algorithm1(ctx, K.ideal(3), ray=ray, block_size=2)
sample = sampler.sample(RandomBits(seed=7))
assert sampler.verify(sample.element)
print(sample.element)       # 精确 alpha
print(sample.point(ctx))    # Algorithm 1 的 x*alpha
print(sample.report)
```

同一个 `Algorithm1` 对象可以重复采样，重用格基和网格准备数据。一般缩放通过 `x` 传入，每个无限位一个数：实位用精确有理数，复位用 `(实部, 虚部)`。复位仅选择共轭对中虚部为正的嵌入；完整坐标按“实位，复位的实部/虚部”排列。

`gamma` 可以是分式元素。若其在有限模的某个素因子处不整，则相应采样集合为空，抛出 `EmptySupport`。同余是局部赋值条件，不能用“差属于整模理想”替代所有分式情形。

**Algorithm 2：自动校准后运行**

```python
from rigorous_cnt import algorithm2, calibrate_walk

# 显式执行较重的参数预处理；结果可重用。
walk = calibrate_walk(ctx, ray, epsilon=QQ(1)/100)
result = algorithm2(ctx, K.ideal(3), walk, ray=ray, rng=RandomBits(seed=7))
assert result.report['mixing_guarantee'] == 'verified_spectral_bound'
print(walk.prime_bound, walk.steps, walk.mixing_l1_bound)
```

校准器从经过 `proof=True` 认证的 PARI 数据构造射线单位格和字符，使用有向球计算筛选所有相关低频字符并上界特征值；它有显式群大小、频率数量和素理想界预算，超出预算抛出 `ResourceLimit`。校准证据绑定数域上下文、模、步数与范数界，不能移用到其他参数。

**Algorithm 2：未经校准的流程实验**

```python
from rigorous_cnt import algorithm2, WalkParameters

result = algorithm2(
    ctx, K.ideal(3),
    walk=WalkParameters.from_prime_bound(ctx, ray, prime_bound=19),
    ray=ray, epsilon=QQ(1)/100,
    rng=RandomBits(seed=7),
    require_mixing_bound=False,  # 明确接受 B 尚无混合保证
)
assert result.element in K.ideal(3)
assert ray.contains(result.element)
assert K.ideal(result.element) == result.output_ideal() * result.input_ideal
print(result.report['mixing_guarantee'])  # not_established
```

`result.output_ideal()` 是 `(beta) * input_ideal^(-1)`，属于输入理想类的逆类。输出是候选理想，不保证每一个候选都是平滑或近素；目标理想族的检测由应用层处理。

若有外部推导的参数，可以提供 `WalkParameters(..., mixing_l1_bound=..., justification=...)`，误差不得超过 `epsilon/2`。运行报告标为 `conditional_on_supplied_bound`，表示以调用方所给的数学依据为前提，程序不会仅凭一段说明文字把它当作已验证证明。

一般子群还需提供 `subgroup_index` 和 `subgroup_oracle(ideal)`；其指数、同态性质以及对应的有限指数闭子群是调用方的前提。程序检查抽中的素理想满足该 oracle，并检查 tau 对应的条件。

**验证与实现结构**

测试包括小网格全部提议点的概率质量穷举、分裂/惰性/分歧素理想的均匀质量检查、格变换与 HKZ 条件、Gaussian 支持界、精确球端点、有理扰动归一化、分式 CRT、复嵌入、重放、资源失败及未证明混合界的拒绝行为。测试还禁止核心采样路径调用 `pari_bnf`；`calibrate_walk` 中的类群/单位群计算属于显式的独立预处理。

| 文件 | 对应内容 |
| --- | --- |
| `arithmetic.py` | 数域上下文、嵌入、局部射线条件、分式 CRT |
| `geometry.py` | Proposition 8.8、实区间与复圆盘的网格采样 |
| `gaussian.py`、`numerics.py` | Lemma 2.22、精确随机决定、带误差界的超越函数 |
| `lattice.py` | 精确 LLL、HKZ 和所需短基的构造与检查 |
| `primes.py` | Lemma 5.4 / Appendix A.2 |
| `parameters.py` | 有理网格参数、随机游走参数与未证明前提 |
| `calibration.py` | 小数域的显式字符与 Gaussian 尾项校准 |
| `sampling.py` | Algorithm 1–2 的完整流程及结果对象 |

本文档介绍 Algorithm 1–2；Algorithm 3–6 的用法见[主文档](../README.md)。当前不含非极大整环扩展。相对数域应先转为绝对数域表示。基础数域、极大整环和精确代数数运算复用 SageMath；预处理成本与论文核心算法成本分开理解。
