#!/usr/bin/env python
"""
GPBB Algorithm Base Class - Refactored Version
Provides common interfaces and shared functionality for bond length adjustment
Enhanced with unified molecule detection and analysis
"""

import os
import copy
import shutil
import numpy as np
import yaml
import logging
import json
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Optional, Set, Any
from multiprocessing import Process
from collections import Counter
from ase import Atoms
from ase.io import read, write, Trajectory


class MoleculeAnalyzer:
    """Unified molecule detection and analysis functionality"""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        
    def detect_and_analyze(self, atoms: Atoms) -> List[Dict[str, Any]]:
        """
        Main entry point for molecule detection and analysis
        Returns list of detected species (molecules and single atoms)
        """
        if not self.config.get('enable_molecule_protection', False):
            return []
            
        # Detect connected components
        molecules = self._detect_molecules(atoms)
        
        # Analyze each detected species
        results = []
        for molecule in molecules:
            info = self._analyze_species(molecule, atoms)
            results.append(info)
            
        return results
    
    def _detect_molecules(self, atoms: Atoms) -> List[Set[int]]:
        """Detect molecular units in structure using connectivity"""
        n_atoms = len(atoms)
        positions = atoms.get_positions()
        symbols = atoms.symbols
        cell = atoms.get_cell()
        pbc = atoms.pbc
        
        molecular_elements = set(self.config['molecular_elements'])
        cutoff = self.config['molecule_detection_threshold']
        min_size = self.config.get('min_molecule_size', 2)
        max_size = self.config.get('max_molecule_size', 20)
        
        # Build adjacency list for molecular atoms
        adjacency = [[] for _ in range(n_atoms)]
        
        for i in range(n_atoms):
            if symbols[i] not in molecular_elements:
                continue
                
            for j in range(i + 1, n_atoms):
                if symbols[j] not in molecular_elements:
                    continue
                
                # Calculate distance with PBC
                distance = self._calculate_distance(
                    positions[i], positions[j], cell, pbc
                )
                
                if distance <= cutoff:
                    adjacency[i].append(j)
                    adjacency[j].append(i)
        
        # Find connected components using DFS
        molecules = []
        visited = [False] * n_atoms
        
        for start in range(n_atoms):
            if visited[start] or symbols[start] not in molecular_elements:
                continue
            
            # DFS to find connected component
            component = []
            stack = [start]
            
            while stack:
                node = stack.pop()
                if not visited[node]:
                    visited[node] = True
                    component.append(node)
                    stack.extend(adjacency[node])
            
            # Add if within size limits (now supports single atoms when min_size=1)
            size = len(component)
            if min_size <= size <= max_size:
                molecules.append(set(component))
        
        return molecules
    
    def _analyze_species(self, indices: Set[int], atoms: Atoms) -> Dict[str, Any]:
        """Analyze a detected species (molecule or single atom)"""
        symbols = atoms.symbols
        positions = atoms.get_positions()
        masses = atoms.get_masses()
        
        # Get composition
        composition = Counter(symbols[i] for i in indices)
        
        # Calculate center of mass
        species_positions = [positions[i] for i in indices]
        species_masses = [masses[i] for i in indices]
        total_mass = sum(species_masses)
        
        if len(species_positions) == 1:
            # Single atom
            com = species_positions[0]
            radius = 0.0
        else:
            # Multi-atom molecule
            com = np.average(species_positions, weights=species_masses, axis=0)
            distances = [np.linalg.norm(pos - com) for pos in species_positions]
            radius = max(distances)
        
        # Generate formula and identify type
        formula = self._generate_formula(composition)
        species_type = self._identify_type(composition)
        
        return {
            'indices': sorted(list(indices)),
            'composition': dict(composition),
            'formula': formula,
            'type': species_type,
            'n_atoms': len(indices),
            'center': com.tolist(),
            'radius': float(radius),
            'mass': float(total_mass)
        }
    
    def _generate_formula(self, composition: Counter) -> str:
        """Generate chemical formula"""
        formula = ""
        
        # Standard order: C, H, then alphabetical
        if 'C' in composition:
            count = composition['C']
            formula += f"C{count if count > 1 else ''}"
        if 'H' in composition:
            count = composition['H']
            formula += f"H{count if count > 1 else ''}"
        
        for element in sorted(composition.keys()):
            if element not in ['C', 'H']:
                count = composition[element]
                formula += f"{element}{count if count > 1 else ''}"
        
        return formula if formula else "Unknown"
    
    def _identify_type(self, composition: Counter) -> str:
        """Identify species type"""
        # Single atoms
        if len(composition) == 1 and sum(composition.values()) == 1:
            element = list(composition.keys())[0]
            return f"Single Atom Adsorbate ({element})"
        
        # Known molecular patterns
        patterns = {
            'CO': {'C': 1, 'O': 1},
            'CO2': {'C': 1, 'O': 2},
            'H2O': {'H': 2, 'O': 1},
            'H3O': {'H': 3, 'O': 1},  # H3O+ recognition
            'CH4': {'C': 1, 'H': 4},
            'NH3': {'N': 1, 'H': 3},
            'O2': {'O': 2},
            'H2': {'H': 2},
        }
        
        for name, pattern in patterns.items():
            if dict(composition) == pattern:
                return name
        
        # General categories
        if 'C' in composition and 'H' in composition:
            return "Hydrocarbon"
        elif 'O' in composition:
            return "Oxide"
        
        return "Unknown"
    
    def _calculate_distance(self, pos1: np.ndarray, pos2: np.ndarray,
                          cell: np.ndarray, pbc: np.ndarray) -> float:
        """Calculate minimum image distance between two positions"""
        delta = pos2 - pos1
        
        if np.any(pbc):
            delta = self._apply_pbc(delta, cell)
        
        return np.linalg.norm(delta)
    
    def _apply_pbc(self, vector: np.ndarray, cell: np.ndarray) -> np.ndarray:
        """Apply periodic boundary conditions to vector"""
        try:
            # Convert to fractional coordinates
            inv_cell = np.linalg.inv(cell.T)
            frac = np.dot(vector, inv_cell)
            # Apply minimum image convention
            frac = frac - np.round(frac)
            # Convert back to Cartesian
            return np.dot(frac, cell.T)
        except np.linalg.LinAlgError:
            return vector
    
    def save_analysis_results(self, results: List[Dict[str, Any]], 
                            output_dir: str, structure_idx: int) -> None:
        """Save molecule analysis results to file"""
        if not results:
            return
            
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{structure_idx}.out")
        
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)


class BaseGPBB(ABC):
    """Base class for GPBB algorithms with common functionality"""
    
    # Default configuration values
    DEFAULT_CONFIG = {
        'log_level': 'INFO',
        'enable_molecule_protection': True,
        'molecule_detection_threshold': 1.5,
        'min_molecule_size': 2,  # Set to 1 to detect single atoms
        'max_molecule_size': 20,
        'molecular_elements': ['C', 'H', 'O', 'N', 'S', 'P'],
        'tolerance': 0.10,
        'confidence_level': 0.80,
        'steps': 2000,
        'num_cores': 1,
        'convergence_check_interval': 50,
        'output_molecule_analysis': False,  # New option
        'molecule_analysis_dir': 'molecule_analysis',  # New option
    }
    
    REQUIRED_CONFIG_KEYS = ['filename', 'elements', 'scale_factors']
    
    def __init__(self, config_file: str):
        """Initialize GPBB algorithm with configuration"""
        self.config = self._load_and_validate_config(config_file)
        self.logger = self._setup_logger()
        self.output_dir = self._setup_output_directory()
        self.molecule_analyzer = MoleculeAnalyzer(self.config)
        
    def _load_and_validate_config(self, config_file: str) -> Dict[str, Any]:
        """Load and validate configuration file"""
        # Load config
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Validate required keys
        missing_keys = [key for key in self.REQUIRED_CONFIG_KEYS if key not in config]
        if missing_keys:
            raise ValueError(f"Missing required configuration keys: {missing_keys}")
        
        # Apply defaults
        for key, value in self.DEFAULT_CONFIG.items():
            config.setdefault(key, value)
        
        # Algorithm-specific validation
        self._validate_algorithm_config(config)
        
        return config
    
    @abstractmethod
    def _validate_algorithm_config(self, config: Dict[str, Any]) -> None:
        """Validate algorithm-specific configuration"""
        pass
    
    def _setup_logger(self) -> logging.Logger:
        """Configure and return logger instance"""
        logger = logging.getLogger('GPBB')
        
        # Configure log level
        log_levels = {
            'DEBUG': logging.DEBUG,
            'INFO': logging.INFO,
            'WARNING': logging.WARNING,
            'ERROR': logging.ERROR
        }
        log_level = log_levels.get(self.config['log_level'].upper(), logging.INFO)
        logger.setLevel(log_level)
        
        # Clear existing handlers and add file handler
        logger.handlers.clear()
        handler = logging.FileHandler('gpbb.log', mode='a')
        handler.setLevel(log_level)
        handler.setFormatter(
            logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        )
        logger.addHandler(handler)
        logger.propagate = False
        
        return logger
    
    def _setup_output_directory(self) -> str:
        """Create and return output directory path"""
        algorithm_name = self.__class__.__name__.lower()
        base_dir = os.path.dirname(self.config['filename'])
        output_dir = os.path.join(base_dir, f'{algorithm_name}_traj')
        
        # Clean and recreate directory
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)
        
        return output_dir
    
    def run(self) -> None:
        """Execute the GPBB algorithm"""
        self._log_run_info()
        
        # Load structures
        images = read(self.config['filename'], ':')
        self.logger.info(f"Loaded {len(images)} structures")
        
        # Process structures in parallel
        self._process_structures_parallel(images)
    
    def _log_run_info(self) -> None:
        """Log algorithm run information"""
        self.logger.info("=" * 60)
        self.logger.info(f"{self.__class__.__name__} - ML Dataset Conversion")
        self.logger.info(f"Molecule Protection: {'Enabled' if self.config['enable_molecule_protection'] else 'Disabled'}")
        self.logger.info(f"Min Molecule Size: {self.config['min_molecule_size']} ({'Single atoms included' if self.config['min_molecule_size'] == 1 else 'Multi-atom only'})")
        self.logger.info(f"Confidence Target: {self.config['confidence_level']}")
        if self.config.get('output_molecule_analysis', False):
            self.logger.info(f"Molecule Analysis Output: {self.config['molecule_analysis_dir']}/")
        self._log_algorithm_info()
        self.logger.info("=" * 60)
    
    @abstractmethod
    def _log_algorithm_info(self) -> None:
        """Log algorithm-specific information"""
        pass
    
    def _process_structures_parallel(self, images: List[Atoms]) -> None:
        """Process structures using multiple processes"""
        num_cores = min(self.config['num_cores'], len(images))
        
        # Split work and create processes
        task_splits = self._split_tasks(len(images), num_cores)
        processes = []
        traj_files = []
        
        for rank, (start, end) in enumerate(task_splits):
            traj_file = os.path.join(self.output_dir, f'rank_{rank}_{start}_{end}.traj')
            traj_files.append(traj_file)
            
            process = Process(
                target=self._process_batch,
                args=(rank, images[start:end], traj_file, start)
            )
            process.start()
            processes.append(process)
        
        # Wait for completion
        for process in processes:
            process.join()
        
        # Combine results
        self._combine_results(traj_files)
    
    def _split_tasks(self, total_tasks: int, num_workers: int) -> List[Tuple[int, int]]:
        """Split tasks evenly among workers"""
        base_size = total_tasks // num_workers
        remainder = total_tasks % num_workers
        
        splits = []
        start = 0
        
        for i in range(num_workers):
            size = base_size + (1 if i < remainder else 0)
            if size > 0:
                splits.append((start, start + size))
                start += size
        
        return splits
    
    def _process_batch(self, rank: int, images: List[Atoms], 
                      output_file: str, start_idx: int) -> None:
        """Process a batch of structures"""
        logger = self._setup_logger()  # Setup logger for child process
        logger.info(f"Process {rank}: Starting {len(images)} structures")
        
        # Initialize molecule analyzer for this process
        molecule_analyzer = MoleculeAnalyzer(self.config)
        
        for i, image in enumerate(images):
            global_idx = start_idx + i
            try:
                adjusted = self._adjust_structure(image, rank, global_idx, molecule_analyzer)
                
                with Trajectory(output_file, 'a') as traj:
                    traj.write(adjusted)
                
                logger.info(f"Process {rank}: Completed structure {global_idx}")
                
            except Exception as e:
                logger.error(f"Process {rank}: Error in structure {global_idx}: {e}")
                # Write original structure on error
                with Trajectory(output_file, 'a') as traj:
                    traj.write(image)
        
        logger.info(f"Process {rank}: Batch complete")
    
    def _adjust_structure(self, image: Atoms, rank: int, idx: int, 
                         molecule_analyzer: MoleculeAnalyzer) -> Atoms:
        """Adjust structure based on number of elements"""
        unique_symbols = image.symbols.species()
        
        if len(unique_symbols) == 1:
            # Single element: apply lattice adjustment
            return self._adjust_single_element(image)
        else:
            # Multi-element: apply algorithm-specific adjustment
            return self._adjust_multi_element(image, rank, idx, molecule_analyzer)
    
    def _adjust_single_element(self, image: Atoms) -> Atoms:
        """Adjust single-element structure by scaling lattice"""
        adjusted = copy.deepcopy(image)
        original_symbol = list(image.symbols.species())[0]
        
        # Get replacement element
        if original_symbol not in self.config['elements']:
            self.logger.warning(f"No mapping for element {original_symbol}")
            return adjusted
        
        new_symbol = self.config['elements'][original_symbol]
        
        # Get scale factor
        bond_key = f"{new_symbol}-{new_symbol}"
        scale_factor = self.config['scale_factors'].get(bond_key, 1.0)
        
        # Apply changes
        adjusted.symbols[:] = new_symbol
        adjusted.set_cell(image.get_cell() * scale_factor, scale_atoms=True)
        
        return adjusted
    
    @abstractmethod
    def _adjust_multi_element(self, image: Atoms, rank: int, idx: int,
                            molecule_analyzer: MoleculeAnalyzer) -> Atoms:
        """Adjust multi-element structure (algorithm-specific)"""
        pass
    
    def _detect_molecules(self, atoms: Atoms) -> List[Set[int]]:
        """Legacy method - redirect to MoleculeAnalyzer"""
        species_list = self.molecule_analyzer.detect_and_analyze(atoms)
        return [set(species['indices']) for species in species_list]
    
    def _calculate_distance(self, pos1: np.ndarray, pos2: np.ndarray,
                          cell: np.ndarray, pbc: np.ndarray) -> float:
        """Calculate minimum image distance between two positions"""
        return self.molecule_analyzer._calculate_distance(pos1, pos2, cell, pbc)
    
    def _apply_pbc(self, vector: np.ndarray, cell: np.ndarray) -> np.ndarray:
        """Apply periodic boundary conditions to vector"""
        return self.molecule_analyzer._apply_pbc(vector, cell)
    
    def _get_bond_key(self, atom1: str, atom2: str) -> str:
        """Generate sorted bond identifier"""
        return "-".join(sorted([atom1, atom2]))
    
    def _get_scale_factor(self, atom1: str, atom2: str) -> float:
        """Get scale factor for bond type"""
        bond_key = self._get_bond_key(atom1, atom2)
        return self.config['scale_factors'].get(bond_key, 1.0)
    
    def _replace_elements(self, atoms: Atoms) -> None:
        """Replace elements in structure according to mapping"""
        symbol_indices = atoms.symbols.indices()
        
        for original, indices in symbol_indices.items():
            if original in self.config['elements']:
                new_symbol = self.config['elements'][original]
                atoms.symbols[indices] = new_symbol
    
    def _combine_results(self, traj_files: List[str]) -> None:
        """Combine trajectory files into final output"""
        all_structures = []
        
        for traj_file in traj_files:
            if os.path.exists(traj_file):
                try:
                    structures = read(traj_file, ':')
                    all_structures.extend(structures)
                except Exception as e:
                    self.logger.error(f"Error reading {traj_file}: {e}")
        
        if all_structures:
            output_file = os.path.join(
                self.output_dir, 
                f'{self.__class__.__name__}_Adjusted.xyz'
            )
            write(output_file, all_structures, format='extxyz')
            
            self.logger.info(f"Algorithm complete: {len(all_structures)} structures")
            self.logger.info(f"Results saved to: {output_file}")