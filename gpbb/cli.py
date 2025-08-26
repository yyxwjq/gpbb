#!/usr/bin/env python
"""
GPBB - Generalized Potential Bond-length Balancing
Main entry point for bond length adjustment algorithms
"""

import os
import sys
import argparse
import yaml
from pathlib import Path


def print_header(title: str, width: int = 60) -> None:
    """Print formatted header"""
    print("=" * width)
    print(title.center(width))
    print("=" * width)


def validate_config(config_path: str) -> dict:
    """Load and validate configuration file"""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file '{config_path}' not found")
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Validate required fields
    if 'filename' not in config:
        raise ValueError("'filename' not specified in configuration")
    
    if not os.path.exists(config['filename']):
        raise FileNotFoundError(f"Input file '{config['filename']}' not found")
    
    return config


def print_config_summary(config: dict) -> None:
    """Print configuration summary"""
    print(f"Input file: {config['filename']}")
    print(f"Log file: gpbb.log")
    print(f"Cores: {config.get('num_cores', 4)}")
    print(f"Molecule protection: {config.get('enable_molecule_protection', False)}")
    print(f"Adaptive optimization: {config.get('use_adaptive_step', True)}")
    print(f"Target confidence: {config.get('confidence_level', 0.90)}")
    print(f"Tolerance: {config.get('tolerance', 0.05)} Å")


def run_algorithm(algorithm: str, config_file: str) -> None:
    """Run specified algorithm"""
    # Clear previous log
    if os.path.exists('gpbb.log'):
        os.remove('gpbb.log')
    
    if algorithm.lower() == 'bl' or algorithm.lower() == 'bl_adjust':
        from .bl import BLAdjustGPBB
        runner = BLAdjustGPBB(config_file)
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")
    
    runner.run()


def main():
    """Main entry point with argument parsing"""
    parser = argparse.ArgumentParser(
        description='GPBB - Generalized Potential Bond-length Balancing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  gpbb                    # Run with default config.yaml
  gpbb -c my_config.yaml  # Use custom config file
  gpbb --algorithm bl     # Specify algorithm (default: bl)
  
Configuration file must contain:
  - filename: path to input structure file
  - elements: element mapping dictionary
  - scale_factors: bond scale factors
  
See documentation for full configuration options.
        """
    )
    
    parser.add_argument(
        '-c', '--config',
        default='config.yaml',
        help='Configuration file (default: config.yaml)'
    )
    
    parser.add_argument(
        '-a', '--algorithm',
        default='bl',
        choices=['bl', 'bl_adjust'],
        help='Algorithm to use (default: bl)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='GPBB 2.0.0'
    )
    
    args = parser.parse_args()
    
    try:
        # Print header
        print_header("GPBB - Bond Length Adjustment")
        
        # Validate configuration
        config = validate_config(args.config)
        
        # Print configuration summary
        print(f"Configuration: {args.config}")
        print(f"Algorithm: {args.algorithm.upper()}")
        print("-" * 60)
        print_config_summary(config)
        print("=" * 60)
        
        # Run algorithm
        run_algorithm(args.algorithm, args.config)
        
        # Success message
        print("\n" + "=" * 60)
        print("Algorithm completed successfully!")
        print("Results saved in *_traj/ directory")
        print("Check gpbb.log for details")
        print("=" * 60)
        
    except FileNotFoundError as e:
        print(f"\nError: {e}")
        print("\nPlease ensure:")
        print("  1. Configuration file exists")
        print("  2. Input structure file path is correct")
        sys.exit(1)
        
    except ValueError as e:
        print(f"\nConfiguration Error: {e}")
        sys.exit(1)
        
    except KeyboardInterrupt:
        print("\n\nProcess interrupted by user")
        sys.exit(1)
        
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        print("\nCheck gpbb.log for details")
        sys.exit(1)


if __name__ == '__main__':
    main()