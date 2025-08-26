#!/usr/bin/env python
from setuptools import setup, find_packages
from pathlib import Path

# Read long description
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text() if readme_file.exists() else ""

setup(
    name="gpbb",
    version="1.0.0",
    author="GPBB Development Team",
    description="Generalized Potential Bond-length Balancing for ML dataset conversion",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/gpbb",
    
    packages=find_packages(),
    
    install_requires=[
        "numpy>=1.19.0",
        "ase>=3.20.0",
        "pyyaml>=5.3.0",
    ],
    
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
            "mypy>=0.900",
        ],
    },
    
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    
    python_requires=">=3.7",
    
    package_data={
        "gpbb": ["*.yaml"],
    },
    
    # 关键修改：使用正确的入口点
    entry_points={
        "console_scripts": [
            "gpbb=gpbb.cli:main",  # 指向 gpbb/cli.py 中的 main 函数
        ],
    },
)