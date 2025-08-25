# GPBB - Generalized Potential Bond-length Balancing
# GPBB - 通用势能键长平衡算法

[English](#english) | [中文](#chinese)

---

<a name="english"></a>
## English

### Overview

GPBB is a bond length adjustment tool for machine learning dataset conversion. It optimizes interatomic distances in crystal and molecular structures to match target bond length distributions.

### Key Features

- **Intelligent Bond Length Optimization**: Adaptive step-size algorithm for bond length optimization
- **Molecule Protection**: Optional rigid molecule protection to preserve molecular structures
- **Parallel Processing**: Multi-core parallel processing for large batches
- **Flexible Configuration**: Full control via YAML configuration files
- **Detailed Logging**: Comprehensive optimization process logging

### Installation

#### Dependencies

```bash
pip install numpy>=1.19.0 ase>=3.20.0 pyyaml>=5.3.0
```

#### Quick Install

```bash
# Clone or download
git clone https://github.com/yourusername/gpbb.git
cd gpbb

# Install package
pip install -e .

# Or use setup script
chmod +x setup.sh
./setup.sh
```

### Usage

#### 1. Prepare Configuration

Copy `config_template.yaml` to `config.yaml` and modify:

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
```

#### 2. Run Program

```bash
# Use default config.yaml
gpbb

# Specify config file
gpbb -c my_config.yaml

# Show help
gpbb --help
```

#### 3. Check Results

- Output structures saved in `*_traj/` directory
- Detailed logs in `gpbb.log`

### Configuration Parameters

#### Core Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `filename` | Input structure file path | Required |
| `elements` | Element replacement mapping | Required |
| `scale_factors` | Bond length scale factors | Required |
| `tolerance` | Bond length tolerance (Å) | 0.05 |
| `confidence_level` | Target confidence | 0.90 |

#### Optimization Control

| Parameter | Description | Default |
|-----------|-------------|---------|
| `use_adaptive_step` | Enable adaptive step size | true |
| `initial_step_size` | Initial step size | 0.002 |
| `steps` | Maximum optimization steps | 2000 |
| `early_stop_no_improvement` | No improvement stop threshold | 500 |

#### Molecule Protection

| Parameter | Description | Default |
|-----------|-------------|---------|
| `enable_molecule_protection` | Enable molecule protection | false |
| `molecule_detection_threshold` | Molecule detection threshold (Å) | 1.5 |
| `molecular_elements` | Elements that can form molecules | ['C','H','O','N','S','P'] |

### Algorithm Principles

#### BL_adjust Algorithm

1. **Initialization**: Calculate original and target distance matrices
2. **Iterative Optimization**:
   - Calculate current errors
   - Adaptively adjust step size
   - Update atomic positions
   - Check convergence
3. **Convergence**: Based on confidence (fraction of bonds within tolerance)

#### Key Innovations

- **Adaptive Step Size**: Dynamic step adjustment based on optimization history
- **Rigid Molecule Protection**: Move detected molecules as rigid bodies
- **Early Stopping**: Auto-stop when no improvement
- **Divergence Detection**: Auto-revert to best state

### Example

#### Basic Usage

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

#### Enable Molecule Protection

```yaml
enable_molecule_protection: true
molecular_elements: ['C', 'O', 'H']
molecule_detection_threshold: 1.5
```

#### Debug Mode

```yaml
log_level: 'DEBUG'
convergence_check_interval: 10
```

### Troubleshooting

#### Common Issues

1. **Slow Convergence**:
   - Decrease `confidence_level`
   - Increase `tolerance`
   - Adjust step size parameters

2. **Inaccurate Results**:
   - Decrease `tolerance`
   - Increase `steps`
   - Verify `scale_factors` are correct

3. **Out of Memory**:
   - Reduce `num_cores`
   - Process structures in batches

### Performance Tips

1. **Parallel Processing**: Set `num_cores` based on CPU cores
2. **Distance Cutoff**: Set reasonable `evaluation_distance_cutoff`
3. **Log Level**: Use `INFO` for production, `DEBUG` for debugging

---

<a name="chinese"></a>
## 中文

### 概述

GPBB 是一个用于机器学习数据集转换的键长调整工具。它通过优化原子间距离来调整晶体和分子结构，使其符合目标键长分布。

### 主要特性

- **智能键长优化**：使用自适应步长算法优化键长
- **分子保护**：可选的分子刚体保护，保持分子内部结构
- **并行处理**：支持多核并行处理大批量结构
- **灵活配置**：通过 YAML 配置文件控制所有参数
- **详细日志**：提供优化过程的详细日志记录

### 安装

#### 依赖项

```bash
pip install numpy>=1.19.0 ase>=3.20.0 pyyaml>=5.3.0
```

#### 快速安装

```bash
# 克隆或下载代码
git clone https://github.com/yourusername/gpbb.git
cd gpbb

# 安装包
pip install -e .

# 或使用安装脚本
chmod +x setup.sh
./setup.sh
```

### 使用方法

#### 1. 准备配置文件

复制 `config_template.yaml` 为 `config.yaml` 并根据需要修改：

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
```

#### 2. 运行程序

```bash
# 使用默认配置文件 (config.yaml)
gpbb

# 指定配置文件
gpbb -c my_config.yaml

# 显示帮助
gpbb --help
```

#### 3. 查看结果

- 输出结构保存在 `*_traj/` 目录
- 详细日志记录在 `gpbb.log`

### 配置参数详解

#### 核心参数

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `filename` | 输入结构文件路径 | 必需 |
| `elements` | 元素替换映射 | 必需 |
| `scale_factors` | 键长缩放因子 | 必需 |
| `tolerance` | 键长容差 (Å) | 0.05 |
| `confidence_level` | 目标置信度 | 0.90 |

#### 优化控制

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `use_adaptive_step` | 启用自适应步长 | true |
| `initial_step_size` | 初始步长 | 0.002 |
| `steps` | 最大优化步数 | 2000 |
| `early_stop_no_improvement` | 无改善停止阈值 | 500 |

#### 分子保护

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `enable_molecule_protection` | 启用分子保护 | false |
| `molecule_detection_threshold` | 分子检测阈值 (Å) | 1.5 |
| `molecular_elements` | 可形成分子的元素 | ['C','H','O','N','S','P'] |

### 算法原理

#### BL_adjust 算法

1. **初始化**：计算原始距离矩阵和目标距离矩阵
2. **迭代优化**：
   - 计算当前误差
   - 自适应调整步长
   - 更新原子位置
   - 检查收敛条件
3. **收敛判定**：基于置信度（容差内键的比例）

#### 关键创新

- **自适应步长**：根据优化历史动态调整步长
- **分子刚体保护**：将检测到的分子作为刚体移动
- **早停机制**：长时间无改善时自动停止
- **发散检测**：自动恢复到最佳状态

### 示例

#### 基础使用

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

#### 启用分子保护

```yaml
enable_molecule_protection: true
molecular_elements: ['C', 'O', 'H']
molecule_detection_threshold: 1.5
```

#### 调试模式

```yaml
log_level: 'DEBUG'
convergence_check_interval: 10
```

### 故障排除

#### 常见问题

1. **收敛慢**：
   - 降低 `confidence_level`
   - 增加 `tolerance`
   - 调整步长参数

2. **结果不准确**：
   - 减小 `tolerance`
   - 增加 `steps`
   - 检查 `scale_factors` 是否正确

3. **内存不足**：
   - 减少 `num_cores`
   - 分批处理结构

### 性能优化建议

1. **并行处理**：根据 CPU 核心数设置 `num_cores`
2. **距离截断**：设置合理的 `evaluation_distance_cutoff`
3. **日志级别**：生产环境使用 `INFO`，调试使用 `DEBUG`

---

## Important Note / 重要说明

**English**: The `scale_factors` must use the ORIGINAL element names (before replacement), not the target element names after mapping!

**中文**: `scale_factors` 必须使用原始元素名称（替换前的），而不是映射后的目标元素名称！

Example / 示例:
- ✅ Correct / 正确: `'Pd-Pd': 0.927` 
- ❌ Wrong / 错误: `'Cu-Cu': 0.927` (if Pd→Cu mapping)

---

## Code Structure / 代码结构

```
gpbb/
├── base.py           # Base class / 基类
├── bl.py            # BL_adjust algorithm / BL_adjust算法
├── gpbb.py          # Main entry / 主程序入口
├── config.yaml      # Configuration / 配置文件
├── setup.py         # Installation / 安装脚本
└── README.md        # This document / 本文档
```

## Version History / 版本历史

- **v1.0.0** (2025-08)
  - Supports multiple architectures / 支持多种体系
  - logging system / 日志系统
  - More configuration options / 更多的配置选项


## License / 许可证

MIT License

## Contributing / 贡献

Welcome Issues and Pull Requests! / 欢迎提交 Issue 和 Pull Request！

## Contact / 联系方式

Please contact via GitHub Issues. / 请通过 GitHub Issues 联系。s