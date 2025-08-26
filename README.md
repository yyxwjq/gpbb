# GPBB - Generalized Potential Bond-length Balancing
# GPBB - 通用势能键长平衡算法

[English](#english) | [中文](#chinese)

---

<a name="english"></a>
## English

### Overview

GPBB is a bond length adjustment tool for machine learning dataset conversion. It optimizes interatomic distances in crystal and molecular structures to match target bond length distributions, with advanced molecule detection and analysis capabilities.

### Key Features

- **Intelligent Bond Length Optimization**: Adaptive step-size algorithm with state tracking
- **Unified Molecule Analysis**: Advanced detection of molecules and single atom adsorbates
- **Molecule Protection**: Optional rigid molecule protection to preserve molecular structures
- **Parallel Processing**: Multi-core parallel processing for large datasets
- **Flexible Configuration**: Full control via YAML configuration files
- **Comprehensive Logging**: Detailed optimization process logging and analysis output

### Installation

#### Dependencies

```bash
pip install numpy>=1.19.0 ase>=3.20.0 pyyaml>=5.3.0
```

#### Quick Install

```bash
# Clone repository
git clone https://github.com/yourusername/gpbb.git
cd gpbb

# Install package
pip install -e .
```

### Usage

GPBB provides two main commands: `adjust` for bond length optimization and `detect` for molecule analysis.

#### 1. Bond Length Adjustment

**Prepare Configuration**

Create `config.yaml`:

```yaml
# Required parameters
filename: 'input.xyz'          # Input structure file
elements:                      # Element mapping
  'Pd': 'Cu'
  'Au': 'Cu'
scale_factors:                 # Bond scale factors (use ORIGINAL elements)
  'Pd-Pd': 0.927              # NOT 'Cu-Cu'!
  'Au-Au': 0.885
  'C-O': 1.000

# Optimization parameters
tolerance: 0.05               # Tolerance (Å)
confidence_level: 0.90        # Target confidence
steps: 2000                   # Max steps

# Molecule detection (NEW)
min_molecule_size: 2          # Set to 1 to detect single atoms
output_molecule_analysis: false  # Save detailed analysis
molecule_analysis_dir: "molecule_analysis"  # Output directory
```

**Run Adjustment**

```bash
# Use default config.yaml
gpbb

# Or explicitly
gpbb adjust

# Specify config file
gpbb adjust -c my_config.yaml

# Show help
gpbb adjust --help
```

#### 2. Molecule Detection

**Basic Detection**

```bash
# Detect molecules only (min size = 2 atoms)
gpbb detect structure.xyz

# Include single atom adsorbates
gpbb detect structure.xyz -mins 1

# Save analysis to directory
gpbb detect structure.xyz -mins 1 -o analysis_results/
```

**Advanced Detection**

```bash
# Custom parameters
gpbb detect structure.xyz -t 1.8 -e C O H N -mins 1 -maxs 50

# Use config file settings
gpbb detect structure.xyz -c config.yaml -o results/

# Verbose output
gpbb detect structure.xyz -mins 1 -v
```

**Detection Parameters**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-t, --threshold` | Detection threshold (Å) | 1.5 |
| `-e, --elements` | Molecular elements | ['C','H','O','N','S','P'] |
| `-mins` | Min molecule size (1=single atoms) | 2 |
| `-maxs` | Max molecule size | 50 |
| `-o, --output` | Output directory | None |
| `-v, --verbose` | Verbose output | False |

#### 3. Check Results

**Adjustment Results:**
- Output structures: `*_traj/` directory
- Detailed logs: `gpbb.log`
- Molecule analysis: `molecule_analysis/` (if enabled)

**Detection Results:**
- Console output with species statistics
- JSON files: `output_dir/{index}.out` (if `-o` specified)

### Configuration Parameters

#### Core Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `filename` | Input structure file path | Required |
| `elements` | Element replacement mapping | Required |
| `scale_factors` | Bond length scale factors | Required |
| `tolerance` | Bond length tolerance (Å) | 0.10 |
| `confidence_level` | Target confidence | 0.80 |

#### Molecule Detection (NEW)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `enable_molecule_protection` | Enable molecule protection | true |
| `molecule_detection_threshold` | Detection threshold (Å) | 1.5 |
| `min_molecule_size` | Min atoms per molecule (1=single atoms) | 2 |
| `max_molecule_size` | Max atoms per molecule | 20 |
| `molecular_elements` | Elements that can form molecules | ['C','H','O','N','S','P'] |
| `output_molecule_analysis` | Save detailed analysis | false |
| `molecule_analysis_dir` | Analysis output directory | "molecule_analysis" |

#### Optimization Control

| Parameter | Description | Default |
|-----------|-------------|---------|
| `use_adaptive_step` | Enable adaptive step size | true |
| `initial_step_size` | Initial step size | 0.002 |
| `steps` | Maximum optimization steps | 2000 |
| `early_stop_no_improvement` | No improvement stop threshold | 500 |
| `evaluation_distance_cutoff` | Only optimize bonds within distance (Å) | 6.0 |

### Molecule Detection Capabilities

#### Supported Species

**Known Molecules:**
- H2O (Water)
- H3O (Hydronium ion) - NEW
- CO, CO2 (Carbon oxides)
- CH4 (Methane)
- NH3 (Ammonia)
- O2, H2 (Diatomic gases)

**General Categories:**
- Hydrocarbons (C-H containing)
- Oxides (O containing)
- Single Atom Adsorbates (when `min_molecule_size: 1`)

#### Analysis Output

Each detected species includes:
- Atomic indices and composition
- Chemical formula and type identification
- Center of mass coordinates
- Molecular radius and mass
- Classification (molecule vs. single atom)

### Algorithm Principles

#### BL_adjust Algorithm

1. **Initialization**: Calculate distance matrices and detect molecules
2. **State Tracking**: Monitor optimization progress with best state preservation
3. **Iterative Optimization**:
   - Calculate current bond errors
   - Adaptively adjust step size
   - Update positions with molecule protection
   - Track convergence metrics
4. **Convergence**: Based on confidence (fraction of bonds within tolerance)

#### Key Innovations

- **Best State Tracking**: Always returns the best configuration found
- **Unified Molecule Analysis**: Single system for detection and protection
- **Single Atom Detection**: Identifies isolated adsorbate species
- **Enhanced Logging**: Comprehensive progress tracking
- **Flexible Output**: Optional detailed molecule analysis files

### Examples

#### Basic Bond Length Adjustment

```yaml
# config.yaml
filename: 'structures.xyz'
elements:
  'Pd': 'Cu'
scale_factors:
  'Pd-Pd': 0.927    # Use original element names!
tolerance: 0.05
confidence_level: 0.90
```

```bash
gpbb adjust -c config.yaml
```

#### With Single Atom Detection

```yaml
# Enable single atom detection and analysis output
min_molecule_size: 1
output_molecule_analysis: true
molecule_analysis_dir: "single_atom_analysis"
enable_molecule_protection: true
```

#### Standalone Molecule Detection

```bash
# Comprehensive analysis including single atoms
gpbb detect my_structures.xyz -mins 1 -o detailed_analysis/

# Quick molecule-only scan
gpbb detect my_structures.xyz -mins 2 -v
```

#### High-Precision Optimization

```yaml
tolerance: 0.03
confidence_level: 0.95
steps: 5000
convergence_check_interval: 25
output_molecule_analysis: true
```

### Troubleshooting

#### Common Issues

1. **Slow Convergence**:
   - Decrease `confidence_level` to 0.80 or lower
   - Increase `tolerance` to 0.08-0.10
   - Reduce `evaluation_distance_cutoff`

2. **Inaccurate Results**:
   - Decrease `tolerance` to 0.03-0.05
   - Increase `steps` to 3000-5000
   - Verify `scale_factors` use original element names

3. **Memory Issues**:
   - Reduce `num_cores`
   - Process structures in smaller batches
   - Increase `evaluation_distance_cutoff` to limit bonds

4. **Molecule Detection Issues**:
   - Adjust `molecule_detection_threshold` (1.0-2.0 Å)
   - Check `molecular_elements` list
   - Use `-v` flag for debugging

### Performance Tips

1. **Parallel Processing**: Set `num_cores` to available CPU cores
2. **Distance Optimization**: Use reasonable `evaluation_distance_cutoff` (4-8 Å)
3. **Logging**: Use `INFO` for production, `DEBUG` for troubleshooting
4. **Memory Management**: For large datasets, process in batches

---

<a name="chinese"></a>
## 中文

### 概述

GPBB 是一个用于机器学习数据集转换的键长调整工具。它通过优化原子间距离来调整晶体和分子结构，使其符合目标键长分布，具备先进的分子检测和分析功能。

### 主要特性

- **智能键长优化**：带状态跟踪的自适应步长算法
- **统一分子分析**：先进的分子和单原子吸附物种检测
- **分子保护**：可选的分子刚体保护，保持分子内部结构
- **并行处理**：支持多核并行处理大型数据集
- **灵活配置**：通过 YAML 配置文件控制所有参数
- **全面日志**：详细的优化过程日志和分析输出

### 安装

#### 依赖项

```bash
pip install numpy>=1.19.0 ase>=3.20.0 pyyaml>=5.3.0
```

#### 快速安装

```bash
# 克隆代码库
git clone https://github.com/yourusername/gpbb.git
cd gpbb

# 安装包
pip install -e .
```

### 使用方法

GPBB 提供两个主要命令：`adjust` 用于键长优化，`detect` 用于分子分析。

#### 1. 键长调整

**准备配置文件**

创建 `config.yaml`：

```yaml
# 必需参数
filename: 'input.xyz'          # 输入结构文件
elements:                      # 元素映射
  'Pd': 'Cu'
  'Au': 'Cu'
scale_factors:                 # 键长缩放因子（使用原始元素名！）
  'Pd-Pd': 0.927              # 不是 'Cu-Cu'！
  'Au-Au': 0.885
  'C-O': 1.000

# 优化参数
tolerance: 0.05               # 容差 (Å)
confidence_level: 0.90        # 目标置信度
steps: 2000                   # 最大步数

# 分子检测（新功能）
min_molecule_size: 2          # 设为1可检测单原子
output_molecule_analysis: false  # 保存详细分析
molecule_analysis_dir: "molecule_analysis"  # 输出目录
```

**运行调整**

```bash
# 使用默认配置文件
gpbb

# 或显式指定
gpbb adjust

# 指定配置文件
gpbb adjust -c my_config.yaml

# 显示帮助
gpbb adjust --help
```

#### 2. 分子检测

**基础检测**

```bash
# 仅检测分子（最小尺寸 = 2个原子）
gpbb detect structure.xyz

# 包含单原子吸附物种
gpbb detect structure.xyz -mins 1

# 保存分析到目录
gpbb detect structure.xyz -mins 1 -o analysis_results/
```

**高级检测**

```bash
# 自定义参数
gpbb detect structure.xyz -t 1.8 -e C O H N -mins 1 -maxs 50

# 使用配置文件设置
gpbb detect structure.xyz -c config.yaml -o results/

# 详细输出
gpbb detect structure.xyz -mins 1 -v
```

**检测参数**

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `-t, --threshold` | 检测阈值 (Å) | 1.5 |
| `-e, --elements` | 分子元素 | ['C','H','O','N','S','P'] |
| `-mins` | 最小分子尺寸（1=单原子） | 2 |
| `-maxs` | 最大分子尺寸 | 50 |
| `-o, --output` | 输出目录 | 无 |
| `-v, --verbose` | 详细输出 | False |

#### 3. 查看结果

**调整结果：**
- 输出结构：`*_traj/` 目录
- 详细日志：`gpbb.log`
- 分子分析：`molecule_analysis/`（如果启用）

**检测结果：**
- 控制台输出物种统计
- JSON文件：`output_dir/{index}.out`（如果指定了 `-o`）

### 配置参数详解

#### 核心参数

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `filename` | 输入结构文件路径 | 必需 |
| `elements` | 元素替换映射 | 必需 |
| `scale_factors` | 键长缩放因子 | 必需 |
| `tolerance` | 键长容差 (Å) | 0.10 |
| `confidence_level` | 目标置信度 | 0.80 |

#### 分子检测（新功能）

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `enable_molecule_protection` | 启用分子保护 | true |
| `molecule_detection_threshold` | 检测阈值 (Å) | 1.5 |
| `min_molecule_size` | 最小分子原子数（1=单原子） | 2 |
| `max_molecule_size` | 最大分子原子数 | 20 |
| `molecular_elements` | 可形成分子的元素 | ['C','H','O','N','S','P'] |
| `output_molecule_analysis` | 保存详细分析 | false |
| `molecule_analysis_dir` | 分析输出目录 | "molecule_analysis" |

#### 优化控制

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `use_adaptive_step` | 启用自适应步长 | true |
| `initial_step_size` | 初始步长 | 0.002 |
| `steps` | 最大优化步数 | 2000 |
| `early_stop_no_improvement` | 无改善停止阈值 | 500 |
| `evaluation_distance_cutoff` | 仅优化此距离内的键 (Å) | 6.0 |

### 分子检测能力

#### 支持的物种

**已知分子：**
- H2O（水）
- H3O（水合氢离子）- 新增
- CO, CO2（碳氧化物）
- CH4（甲烷）
- NH3（氨）
- O2, H2（双原子气体）

**通用类别：**
- 碳氢化合物（含C-H）
- 氧化物（含O）
- 单原子吸附物种（当 `min_molecule_size: 1` 时）

#### 分析输出

每个检测到的物种包括：
- 原子索引和组成
- 化学式和类型识别
- 质心坐标
- 分子半径和质量
- 分类（分子 vs 单原子）

### 算法原理

#### BL_adjust 算法

1. **初始化**：计算距离矩阵并检测分子
2. **状态跟踪**：监控优化进程并保存最佳状态
3. **迭代优化**：
   - 计算当前键误差
   - 自适应调整步长
   - 带分子保护的位置更新
   - 跟踪收敛指标
4. **收敛判定**：基于置信度（容差内键的比例）

#### 关键创新

- **最佳状态跟踪**：始终返回找到的最佳配置
- **统一分子分析**：检测和保护的单一系统
- **单原子检测**：识别孤立的吸附物种
- **增强日志**：全面的进度跟踪
- **灵活输出**：可选的详细分子分析文件

### 示例

#### 基础键长调整

```yaml
# config.yaml
filename: 'structures.xyz'
elements:
  'Pd': 'Cu'
scale_factors:
  'Pd-Pd': 0.927    # 使用原始元素名！
tolerance: 0.05
confidence_level: 0.90
```

```bash
gpbb adjust -c config.yaml
```

#### 启用单原子检测

```yaml
# 启用单原子检测和分析输出
min_molecule_size: 1
output_molecule_analysis: true
molecule_analysis_dir: "single_atom_analysis"
enable_molecule_protection: true
```

#### 独立分子检测

```bash
# 包括单原子的全面分析
gpbb detect my_structures.xyz -mins 1 -o detailed_analysis/

# 快速分子扫描
gpbb detect my_structures.xyz -mins 2 -v
```

#### 高精度优化

```yaml
tolerance: 0.03
confidence_level: 0.95
steps: 5000
convergence_check_interval: 25
output_molecule_analysis: true
```

### 故障排除

#### 常见问题

1. **收敛慢**：
   - 降低 `confidence_level` 至 0.80 或更低
   - 增加 `tolerance` 至 0.08-0.10
   - 减少 `evaluation_distance_cutoff`

2. **结果不准确**：
   - 减小 `tolerance` 至 0.03-0.05
   - 增加 `steps` 至 3000-5000
   - 验证 `scale_factors` 使用原始元素名

3. **内存问题**：
   - 减少 `num_cores`
   - 分批处理结构
   - 增加 `evaluation_distance_cutoff` 限制键数

4. **分子检测问题**：
   - 调整 `molecule_detection_threshold`（1.0-2.0 Å）
   - 检查 `molecular_elements` 列表
   - 使用 `-v` 标志调试

### 性能优化建议

1. **并行处理**：根据可用CPU核心数设置 `num_cores`
2. **距离优化**：使用合理的 `evaluation_distance_cutoff`（4-8 Å）
3. **日志级别**：生产环境使用 `INFO`，故障排除使用 `DEBUG`
4. **内存管理**：大型数据集分批处理

---

## Important Notes / 重要说明

### Scale Factors / 缩放因子

**English**: The `scale_factors` must use the ORIGINAL element names (before replacement), not the target element names after mapping!

**中文**: `scale_factors` 必须使用原始元素名称（替换前的），而不是映射后的目标元素名称！

Example / 示例:
- ✅ Correct / 正确: `'Pd-Pd': 0.927` 
- ❌ Wrong / 错误: `'Cu-Cu': 0.927` (if Pd→Cu mapping)

### New Features in v2.1 / v2.1新功能

**English**:
- Single atom adsorbate detection (`min_molecule_size: 1`)
- H3O+ molecule recognition
- Unified molecule analysis system
- Enhanced output options (`-o directory/`)
- Simplified command interface (`-mins`, `-maxs`)

**中文**:
- 单原子吸附物种检测（`min_molecule_size: 1`）
- H3O+分子识别
- 统一分子分析系统
- 增强输出选项（`-o directory/`）
- 简化命令接口（`-mins`，`-maxs`）

---

## Code Structure / 代码结构

```
gpbb/
├── base.py           # Base class with unified molecule analysis / 统一分子分析的基类
├── bl.py             # BL_adjust algorithm / BL_adjust算法
├── cli.py            # Command line interface / 命令行接口
├── config.yaml       # Configuration file / 配置文件
├── setup.py          # Installation script / 安装脚本
└── README.md         # This document / 本文档
```

## Version History / 版本历史

- **v2.1.0** (2025-08)
  - Unified molecule detection and analysis system / 统一分子检测分析系统
  - Single atom adsorbate detection / 单原子吸附物种检测
  - H3O+ molecule recognition / H3O+分子识别
  - Enhanced CLI with simplified parameters / 简化参数的增强CLI
  - Best state tracking in optimization / 优化中的最佳状态跟踪

- **v2.0.0** (2025-08)
  - Refactored codebase with improved architecture / 改进架构的重构代码
  - Enhanced molecule protection system / 增强分子保护系统
  - Better state management and convergence / 更好的状态管理和收敛
  - Comprehensive logging system / 全面日志系统

- **v1.0.0** (2025-08)
  - Initial release with basic functionality / 基础功能的初始版本
  - Multi-architecture support / 多架构支持
  - Parallel processing / 并行处理

## License / 许可

MIT License

## Contributing / 贡献

Welcome Issues and Pull Requests! / 欢迎提交 Issue 和 Pull Request！

## Contact / 联系方式

Please contact via GitHub Issues. / 请通过 GitHub Issues 联系。