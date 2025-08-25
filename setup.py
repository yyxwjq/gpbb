#!/usr/bin/env python
"""
Setup script for GPBB package
"""

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
    
    py_modules=["base", "bl"],
    scripts=["gpbb"],
    
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
        "License :: OSI Approved :: MIT License",
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
        "": ["*.yaml", "*.md"],
    },
    
    entry_points={
        "console_scripts": [
            "gpbb=gpbb:main",
        ],
    },
)
