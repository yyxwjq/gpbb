#!/usr/bin/env python
"""
GPBB Command Line Interface
Handles both bond length adjustment and molecule detection
Simplified with unified molecule analysis system
"""

import os
import sys
import yaml
import argparse
import numpy as np
from typing import List, Dict, Set
from collections import Counter
from ase.io import read


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
    
    # Validate required fields for adjustment
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
    print(f"Min molecule size: {config.get('min_molecule_size', 2)} ({'Single atoms included' if config.get('min_molecule_size', 2) == 1 else 'Multi-atom only'})")
    print(f"Adaptive optimization: {config.get('use_adaptive_step', True)}")
    print(f"Target confidence: {config.get('confidence_level', 0.90)}")
    print(f"Tolerance: {config.get('tolerance', 0.05)} Å")
    if config.get('output_molecule_analysis', False):
        print(f"Molecule analysis output: {config.get('molecule_analysis_dir', 'molecule_analysis')}/")


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


def run_adjustment(args):
    """Run bond length adjustment"""
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
        if config.get('output_molecule_analysis', False):
            print(f"Molecule analysis saved in {config.get('molecule_analysis_dir', 'molecule_analysis')}/")
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


def run_detection(args):
    """Run molecule detection using unified analysis system"""
    try:
        print_header("GPBB - Molecule Detection")
        
        # Load configuration if provided
        config = {}
        if args.config and os.path.exists(args.config):
            with open(args.config, 'r') as f:
                full_config = yaml.safe_load(f)
                # Extract molecule-related parameters
                config['molecule_detection_threshold'] = full_config.get(
                    'molecule_detection_threshold', args.threshold)
                config['molecular_elements'] = full_config.get(
                    'molecular_elements', args.elements)
                config['min_molecule_size'] = full_config.get(
                    'min_molecule_size', args.mins)
                config['max_molecule_size'] = full_config.get(
                    'max_molecule_size', args.maxs)
                config['enable_molecule_protection'] = True
        else:
            config = {
                'molecule_detection_threshold': args.threshold,
                'molecular_elements': args.elements,
                'min_molecule_size': args.mins,
                'max_molecule_size': args.maxs,
                'enable_molecule_protection': True
            }
        
        # Import unified analyzer
        from .base import MoleculeAnalyzer
        
        # Read structures
        if args.index == ':':
            structures = read(args.input, ':')
        else:
            structures = read(args.input, index=int(args.index))
            
        if not isinstance(structures, list):
            structures = [structures]
        
        print(f"Input file: {args.input}")
        print(f"Structures to analyze: {len(structures)}")
        print(f"Detection threshold: {config['molecule_detection_threshold']} Å")
        print(f"Molecular elements: {config['molecular_elements']}")
        print(f"Size range: {config['min_molecule_size']}-{config['max_molecule_size']} atoms")
        if config['min_molecule_size'] == 1:
            print("Single atom detection: Enabled")
        print("=" * 60)
        
        analyzer = MoleculeAnalyzer(config)
        
        # Process each structure
        all_results = []
        for idx, atoms in enumerate(structures):
            print(f"\nStructure {idx}: {len(atoms)} atoms")
            print("-" * 40)
            
            # Detect and analyze species
            results = analyzer.detect_and_analyze(atoms)
            all_results.append(results)
            
            if not results:
                print("  No species detected")
                continue
            
            # Separate single atoms and multi-atom molecules
            single_atoms = [r for r in results if r['n_atoms'] == 1]
            molecules = [r for r in results if r['n_atoms'] > 1]
            
            # Display molecules
            if molecules:
                by_formula = {}
                for mol in molecules:
                    formula = mol['formula']
                    if formula not in by_formula:
                        by_formula[formula] = []
                    by_formula[formula].append(mol)
                
                print(f"  Found {len(molecules)} molecule(s):")
                for formula in sorted(by_formula.keys()):
                    mols = by_formula[formula]
                    mol_type = mols[0]['type']
                    print(f"    {formula} ({mol_type}): {len(mols)} molecule(s)")
                    
                    if args.verbose and len(mols) <= 3:
                        for mol in mols:
                            indices = mol['indices']
                            if len(indices) <= 10:
                                print(f"      Atoms: {indices}")
                            else:
                                print(f"      Atoms: {indices[:10]}... ({len(indices)} total)")
                            print(f"      Center: [{mol['center'][0]:.2f}, {mol['center'][1]:.2f}, {mol['center'][2]:.2f}]")
                            print(f"      Radius: {mol['radius']:.2f} Å")
            
            # Display single atoms
            if single_atoms:
                by_element = {}
                for atom in single_atoms:
                    element = atom['formula']
                    if element not in by_element:
                        by_element[element] = []
                    by_element[element].append(atom)
                
                print(f"  Found {len(single_atoms)} single atom adsorbate(s):")
                for element in sorted(by_element.keys()):
                    atoms_list = by_element[element]
                    print(f"    {element}: {len(atoms_list)} atom(s)")
                    
                    if args.verbose:
                        for atom_info in atoms_list:
                            atom_idx = atom_info['indices'][0]
                            pos = atom_info['center']
                            print(f"      Atom {atom_idx}: [{pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f}]")
            
            # Save results if requested
            if args.output:
                output_dir = args.output
                analyzer.save_analysis_results(results, output_dir, idx)
                print(f"  Results saved to: {output_dir}/{idx}.out")
        
        # Overall summary for multiple structures
        if len(structures) > 1:
            print("\n" + "=" * 60)
            print("OVERALL SUMMARY")
            print("-" * 60)
            
            total_molecules = sum(len([r for r in results if r['n_atoms'] > 1]) for results in all_results)
            total_single_atoms = sum(len([r for r in results if r['n_atoms'] == 1]) for results in all_results)
            
            print(f"Total molecules across all structures: {total_molecules}")
            if config['min_molecule_size'] == 1:
                print(f"Total single atom adsorbates across all structures: {total_single_atoms}")
            
            # Collect all formulas
            formula_counts = Counter()
            for results in all_results:
                for species in results:
                    formula_counts[species['formula']] += 1
            
            if formula_counts:
                print("\nMost common species:")
                for formula, count in formula_counts.most_common(10):
                    print(f"  {formula}: {count}")
        
        print("=" * 60)
        
    except Exception as e:
        print(f"\nError: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def main():
    """Main entry point with subcommands"""
    parser = argparse.ArgumentParser(
        description='GPBB - Generalized Potential Bond-length Balancing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Commands:
  adjust  - Run bond length adjustment (default)
  detect  - Detect molecules in structures

Examples:
  gpbb                           # Run adjustment with config.yaml
  gpbb adjust -c my_config.yaml  # Custom config file
  gpbb detect structure.xyz      # Detect molecules
  gpbb detect structure.xyz -t 1.8 -e C O  # Custom parameters
  gpbb detect structure.xyz -mins 1  # Include single atoms
  gpbb detect structure.xyz -o results/  # Save to directory
        """
    )
    
    parser.add_argument('--version', action='version', version='GPBB 2.1.0')
    
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Adjust subcommand
    adjust_parser = subparsers.add_parser(
        'adjust',
        help='Run bond length adjustment'
    )
    adjust_parser.add_argument(
        '-c', '--config',
        default='config.yaml',
        help='Configuration file (default: config.yaml)'
    )
    adjust_parser.add_argument(
        '-a', '--algorithm',
        default='bl',
        choices=['bl', 'bl_adjust'],
        help='Algorithm to use (default: bl)'
    )
    adjust_parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    # Detect subcommand
    detect_parser = subparsers.add_parser(
        'detect',
        help='Detect molecules in structures'
    )
    detect_parser.add_argument('input', help='Input structure file')
    detect_parser.add_argument('-c', '--config', help='Configuration file')
    detect_parser.add_argument(
        '-i', '--index',
        default=':',
        help='Structure index (default: all)'
    )
    detect_parser.add_argument(
        '-t', '--threshold',
        type=float,
        default=1.5,
        help='Detection threshold (Å)'
    )
    detect_parser.add_argument(
        '-e', '--elements',
        nargs='+',
        default=['C', 'H', 'O', 'N', 'S', 'P'],
        help='Molecular elements'
    )
    detect_parser.add_argument(
        '-mins', '--mins',
        type=int,
        default=2,
        help='Min molecule size (set to 1 to include single atoms)'
    )
    detect_parser.add_argument(
        '-maxs', '--maxs',
        type=int,
        default=50,
        help='Max molecule size'
    )
    detect_parser.add_argument(
        '-o', '--output',
        help='Output directory (files named as {index}.out)'
    )
    detect_parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args()
    
    # Default to adjust if no command
    if args.command is None:
        args.command = 'adjust'
        args.config = 'config.yaml'
        args.algorithm = 'bl'
        args.verbose = False
    
    # Execute command
    if args.command == 'adjust':
        run_adjustment(args)
    elif args.command == 'detect':
        run_detection(args)


if __name__ == '__main__':
    main()