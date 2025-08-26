#!/usr/bin/env python
"""
BL_adjust Algorithm - Improved Version with Better State Management
Enhanced best state tracking and clearer convergence logic
Integrated with unified molecule analysis system
"""

import copy
import os
import logging
import numpy as np
from typing import Dict, List, Tuple, Optional, Set, Any
from dataclasses import dataclass, field
from ase import Atoms
from .base import BaseGPBB, MoleculeAnalyzer


@dataclass
class OptimizationState:
    """Container for complete optimization state"""
    positions: np.ndarray
    rmse: float
    confidence: float
    step: int
    max_error: float = 0.0
    mean_error: float = 0.0
    n_converged: int = 0
    n_total: int = 0
    
    def is_better_than(self, other: 'OptimizationState') -> bool:
        """
        Compare optimization states
        Priority: 1. Confidence, 2. RMSE
        """
        if self.confidence > other.confidence:
            return True
        if self.confidence == other.confidence and self.rmse < other.rmse:
            return True
        return False
    
    def to_dict(self) -> dict:
        """Convert to dictionary for logging"""
        return {
            'step': self.step,
            'confidence': self.confidence,
            'rmse': self.rmse,
            'max_error': self.max_error,
            'mean_error': self.mean_error,
            'n_converged': self.n_converged,
            'n_total': self.n_total
        }


@dataclass
class OptimizationTracker:
    """Track optimization progress and best states"""
    initial_state: OptimizationState
    best_state: OptimizationState
    current_state: OptimizationState
    history: List[OptimizationState] = field(default_factory=list)
    no_improvement_count: int = 0
    
    def update(self, new_state: OptimizationState) -> bool:
        """
        Update tracker with new state
        Returns True if new state is better
        """
        self.current_state = new_state
        self.history.append(new_state)
        
        is_better = new_state.is_better_than(self.best_state)
        if is_better:
            # Deep copy positions to ensure we keep the best
            self.best_state = OptimizationState(
                positions=new_state.positions.copy(),
                rmse=new_state.rmse,
                confidence=new_state.confidence,
                step=new_state.step,
                max_error=new_state.max_error,
                mean_error=new_state.mean_error,
                n_converged=new_state.n_converged,
                n_total=new_state.n_total
            )
            self.no_improvement_count = 0
            return True
        else:
            self.no_improvement_count += 1
            return False
    
    def get_best_positions(self) -> np.ndarray:
        """Get a copy of best positions"""
        return self.best_state.positions.copy()
    
    def should_revert(self, divergence_threshold: float) -> bool:
        """Check if should revert to best state"""
        return self.current_state.rmse > self.initial_state.rmse * divergence_threshold


class BLAdjustGPBB(BaseGPBB):
    """Bond Length Adjustment Algorithm with Improved State Management"""
    
    # Algorithm-specific defaults
    BL_DEFAULTS = {
        'evaluation_distance_cutoff': 6.0,
        'use_adaptive_step': True,
        'initial_step_size': 0.002,
        'min_step_size': 0.0001,
        'max_step_size': 0.01,
        'step_reduction_factor': 0.7,
        'step_increase_factor': 1.05,
        'patience': 50,
        'error_history_size': 30,
        'early_stop_no_improvement': 500,
        'divergence_threshold': 3.0,
        'save_best_every': 100,  # Save best state info every N steps
    }
    
    def _validate_algorithm_config(self, config: Dict[str, Any]) -> None:
        """Apply BL_adjust specific defaults"""
        for key, value in self.BL_DEFAULTS.items():
            config.setdefault(key, value)
    
    def _log_algorithm_info(self) -> None:
        """Log BL_adjust specific information"""
        self.logger.info(f"Algorithm: BL_adjust (Bond Length Optimization)")
        self.logger.info(f"Distance Cutoff: {self.config['evaluation_distance_cutoff']} Ã…")
        self.logger.info(f"Adaptive Step: {'Enabled' if self.config['use_adaptive_step'] else 'Disabled'}")
        self.logger.info(f"Initial Step Size: {self.config['initial_step_size']}")
    
    def _adjust_multi_element(self, image: Atoms, rank: int, idx: int,
                            molecule_analyzer: MoleculeAnalyzer) -> Atoms:
        """Adjust multi-element structure using BL optimization"""
        return self._optimize_structure(image, rank, idx, molecule_analyzer)
    
    def _optimize_structure(self, image: Atoms, rank: int, idx: int,
                          molecule_analyzer: MoleculeAnalyzer) -> Atoms:
        """Main optimization routine with improved state tracking"""
        logger = logging.getLogger('GPBB')
        
        # Prepare working copy
        working_image = copy.deepcopy(image)
        original_symbols = list(image.symbols)
        
        # Replace elements
        self._replace_elements(working_image)
        
        # Detect and analyze molecules using unified analyzer
        molecule_analysis = molecule_analyzer.detect_and_analyze(working_image)
        molecules = [set(species['indices']) for species in molecule_analysis]
        
        # Log molecule information
        if molecules:
            self._log_molecules(molecule_analysis, idx)
        
        # Save molecule analysis if requested
        if self.config.get('output_molecule_analysis', False):
            analysis_dir = self.config.get('molecule_analysis_dir', 'molecule_analysis')
            molecule_analyzer.save_analysis_results(molecule_analysis, analysis_dir, idx)
        
        # Initialize optimization
        tracker = self._initialize_tracker(
            image, working_image, original_symbols, idx
        )
        
        # Log initial state
        logger.info(f"Structure {idx}: Starting optimization - "
                   f"Initial RMSE: {tracker.initial_state.rmse:.4f} Ã…, "
                   f"Confidence: {tracker.initial_state.confidence:.4f}, "
                   f"Evaluating {tracker.initial_state.n_total} bonds")
        
        # Run optimization
        final_tracker = self._run_optimization_loop(
            tracker, image, working_image, 
            original_symbols, molecules, idx
        )
        
        # Apply best positions (not current!)
        best_positions = final_tracker.get_best_positions()
        working_image.set_positions(best_positions)
        
        # Log final state
        logger.info(f"Structure {idx}: Optimization complete - "
                   f"Best step: {final_tracker.best_state.step}, "
                   f"Confidence: {final_tracker.best_state.confidence:.4f}, "
                   f"RMSE: {final_tracker.best_state.rmse:.4f} Ã…, "
                   f"Converged: {final_tracker.best_state.n_converged}/{final_tracker.best_state.n_total}")
        
        # Center if non-periodic
        if not np.any(self._has_pbc(working_image)):
            working_image.center(vacuum=8, axis=(0, 1, 2))
        
        return working_image
    
    def _initialize_tracker(self, original: Atoms, working: Atoms,
                          original_symbols: List[str], idx: int) -> OptimizationTracker:
        """Initialize optimization tracker with initial state"""
        logger = logging.getLogger('GPBB')
        
        positions = original.get_positions()
        cell = original.get_cell()
        
        # Calculate distance matrices
        orig_distances = self._calculate_distances(positions, cell)
        target_distances = self._calculate_targets(
            orig_distances, original_symbols, working.symbols
        )
        
        # Setup evaluation mask
        cutoff = self.config['evaluation_distance_cutoff']
        eval_mask = self._create_evaluation_mask(orig_distances, cutoff)
        n_bonds = np.sum(eval_mask)
        
        # Calculate initial metrics
        errors = np.abs(orig_distances - target_distances)[eval_mask]
        
        initial_state = OptimizationState(
            positions=positions.copy(),
            rmse=np.sqrt(np.mean(errors**2)),
            confidence=np.mean(errors <= self.config['tolerance']),
            step=0,
            max_error=np.max(errors) if len(errors) > 0 else 0,
            mean_error=np.mean(errors) if len(errors) > 0 else 0,
            n_converged=np.sum(errors <= self.config['tolerance']),
            n_total=len(errors)
        )
        
        logger.info(f"Structure {idx}: Evaluating {n_bonds} bonds within {cutoff} Ã…")
        
        # Initialize tracker with initial state as best
        tracker = OptimizationTracker(
            initial_state=initial_state,
            best_state=copy.deepcopy(initial_state),
            current_state=initial_state
        )
        
        return tracker
    
    def _run_optimization_loop(self, tracker: OptimizationTracker,
                              original: Atoms, working: Atoms,
                              original_symbols: List[str],
                              molecules: List[Set[int]],
                              idx: int) -> OptimizationTracker:
        """Execute optimization loop with state tracking"""
        logger = logging.getLogger('GPBB')
        
        # Setup
        cell = original.get_cell()
        positions = tracker.initial_state.positions.copy()
        
        orig_distances = self._calculate_distances(
            original.get_positions(), cell
        )
        target_distances = self._calculate_targets(
            orig_distances, original_symbols, working.symbols
        )
        
        cutoff = self.config['evaluation_distance_cutoff']
        eval_mask = self._create_evaluation_mask(orig_distances, cutoff)
        
        # Optimization parameters
        tolerance = self.config['tolerance']
        confidence_target = self.config['confidence_level']
        max_steps = self.config['steps']
        
        # Adaptive step size
        step_size = self.config['initial_step_size']
        error_history = []
        
        # Main optimization loop
        for step in range(1, max_steps + 1):
            # Calculate current state
            distances = self._calculate_distances(positions, cell)
            errors = np.abs(distances - target_distances)
            eval_errors = errors[eval_mask]
            
            if len(eval_errors) == 0:
                logger.warning(f"Structure {idx}: No bonds to evaluate")
                break
            
            # Create current state
            current_state = OptimizationState(
                positions=positions.copy(),
                rmse=np.sqrt(np.mean(eval_errors**2)),
                confidence=np.mean(eval_errors <= tolerance),
                step=step,
                max_error=np.max(eval_errors),
                mean_error=np.mean(eval_errors),
                n_converged=np.sum(eval_errors <= tolerance),
                n_total=len(eval_errors)
            )
            
            # Update tracker
            is_better = tracker.update(current_state)
            
            # Check divergence
            if tracker.should_revert(self.config['divergence_threshold']):
                logger.warning(f"Structure {idx}: Diverged at step {step}, reverting to best state")
                positions = tracker.get_best_positions()
                step_size = self.config['min_step_size']
                error_history.clear()
                continue
            
            # Log if improvement
            if is_better:
                logger.debug(f"Structure {idx}: New best at step {step} "
                           f"(confidence={current_state.confidence:.4f})")
            
            # Check convergence
            if current_state.confidence >= confidence_target:
                logger.info(f"Structure {idx}: Converged at step {step} - "
                          f"Confidence: {current_state.confidence:.4f}, "
                          f"RMSE: {current_state.rmse:.4f} Ã…")
                break
            
            # Early stopping
            if tracker.no_improvement_count > self.config['early_stop_no_improvement']:
                logger.info(f"Structure {idx}: Early stop at step {step} - "
                          f"No improvement for {tracker.no_improvement_count} steps, "
                          f"Best: step {tracker.best_state.step} "
                          f"(Conf: {tracker.best_state.confidence:.4f}, RMSE: {tracker.best_state.rmse:.4f} Ã…)")
                # Ensure we use best positions
                positions = tracker.get_best_positions()
                break
            
            # High confidence early stop
            if current_state.confidence >= 0.95 and tracker.no_improvement_count > 50:
                logger.info(f"Structure {idx}: Stopping at high confidence "
                          f"{current_state.confidence:.4f}")
                break
            
            # Adaptive step size
            if self.config['use_adaptive_step']:
                step_size = self._adapt_step_size(
                    error_history, eval_errors, step_size
                )
                error_history.append(np.mean(eval_errors))
                if len(error_history) > self.config['error_history_size']:
                    error_history.pop(0)
            
            # Calculate and apply updates
            updates = self._calculate_updates(
                positions, distances, target_distances,
                step_size, tolerance, cell,
                molecules, cutoff
            )
            positions += updates
            
            # Periodic logging
            if step % self.config['convergence_check_interval'] == 0:
                logger.info(f"Structure {idx}: Step {step}, "
                          f"Confidence: {current_state.confidence:.4f}, "
                          f"RMSE: {current_state.rmse:.4f} Ã…, "
                          f"Max Error: {current_state.max_error:.4f} Ã…, "
                          f"Step Size: {step_size:.6f}, "
                          f"Best: step {tracker.best_state.step}")
            
            # Save best state info periodically
            if step % self.config.get('save_best_every', 100) == 0:
                logger.debug(f"Structure {idx}: Best state checkpoint at step {step}")
                logger.debug(f"  {tracker.best_state.to_dict()}")
        
        # Final check - ensure we have best positions
        if tracker.no_improvement_count > 0:
            logger.info(f"Structure {idx}: Reverting to best state from step {tracker.best_state.step}")
        
        return tracker
    
    def _calculate_distances(self, positions: np.ndarray, 
                            cell: np.ndarray) -> np.ndarray:
        """Calculate distance matrix with PBC"""
        n = len(positions)
        distances = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i + 1, n):
                delta = positions[j] - positions[i]
                delta = self._apply_pbc(delta, cell)
                dist = np.linalg.norm(delta)
                distances[i, j] = distances[j, i] = dist
        
        return distances
    
    def _calculate_targets(self, distances: np.ndarray,
                         orig_symbols: List[str],
                         new_symbols) -> np.ndarray:
        """Calculate target distances using scale factors"""
        n = len(orig_symbols)
        targets = np.zeros_like(distances)
        
        for i in range(n):
            for j in range(i + 1, n):
                scale = self._get_scale_factor(orig_symbols[i], orig_symbols[j])
                targets[i, j] = targets[j, i] = distances[i, j] * scale
        
        return targets
    
    def _create_evaluation_mask(self, distances: np.ndarray, 
                               cutoff: float) -> np.ndarray:
        """Create mask for bonds within cutoff"""
        n = distances.shape[0]
        mask = (distances <= cutoff) & (distances > 0)
        # Only upper triangle
        mask = mask & np.triu(np.ones_like(mask, dtype=bool), k=1)
        return mask
    
    def _calculate_updates(self, positions: np.ndarray,
                          distances: np.ndarray,
                          targets: np.ndarray,
                          step_size: float,
                          tolerance: float,
                          cell: np.ndarray,
                          molecules: List[Set[int]],
                          cutoff: float) -> np.ndarray:
        """Calculate position updates"""
        n = len(positions)
        updates = np.zeros_like(positions)
        
        for i in range(n):
            for j in range(i + 1, n):
                dist = distances[i, j]
                
                # Skip if beyond cutoff or too close
                if dist > cutoff or dist < 0.5:
                    continue
                
                target = targets[i, j]
                error = abs(target - dist)
                
                # Skip if within tolerance
                if error < tolerance:
                    continue
                
                # Skip intramolecular bonds
                if self._is_intramolecular(i, j, molecules):
                    continue
                
                # Calculate update
                delta = positions[i] - positions[j]
                delta = self._apply_pbc(delta, cell)
                norm = np.linalg.norm(delta)
                
                if norm < 1e-10:
                    continue
                
                unit = delta / norm
                magnitude = min(step_size, error * 0.1)
                
                # Apply based on direction
                if target > dist:  # Lengthen
                    self._apply_update(updates, i, j, unit * magnitude, 
                                     molecules, lengthen=True)
                else:  # Shorten
                    self._apply_update(updates, i, j, unit * magnitude,
                                     molecules, lengthen=False)
        
        # Limit maximum displacement
        max_disp = np.max(np.linalg.norm(updates, axis=1))
        if max_disp > 0.05:
            updates *= 0.05 / max_disp
        
        return updates
    
    def _apply_update(self, updates: np.ndarray, i: int, j: int,
                     delta: np.ndarray, molecules: List[Set[int]],
                     lengthen: bool) -> None:
        """Apply position update with molecule protection"""
        if not self.config['enable_molecule_protection']:
            # Simple update
            sign = 1 if lengthen else -1
            updates[i] += sign * delta
            updates[j] -= sign * delta
            return
        
        # Find molecules containing atoms
        mol_i = next((m for m in molecules if i in m), None)
        mol_j = next((m for m in molecules if j in m), None)
        
        sign = 1 if lengthen else -1
        
        # Apply based on molecule membership
        if mol_i and not mol_j:
            for atom in mol_i:
                updates[atom] += sign * delta
            updates[j] -= sign * delta
        elif mol_j and not mol_i:
            updates[i] += sign * delta
            for atom in mol_j:
                updates[atom] -= sign * delta
        elif mol_i and mol_j and mol_i != mol_j:
            for atom in mol_i:
                updates[atom] += sign * delta
            for atom in mol_j:
                updates[atom] -= sign * delta
        else:
            updates[i] += sign * delta
            updates[j] -= sign * delta
    
    def _adapt_step_size(self, history: List[float], errors: np.ndarray,
                        current_step: float) -> float:
        """Adaptively adjust step size"""
        config = self.config
        
        if len(history) < config['patience']:
            return current_step
        
        recent = history[-config['patience']:]
        trend = (recent[-1] - recent[0]) / max(recent[0], 1e-10)
        
        if trend > 0.01:  # Getting worse
            new_step = current_step * config['step_reduction_factor']
        elif abs(trend) < 0.001:  # Plateau
            new_step = current_step * config['step_increase_factor']
        else:
            return current_step
        
        return np.clip(new_step, config['min_step_size'], config['max_step_size'])
    
    def _is_intramolecular(self, i: int, j: int, 
                          molecules: List[Set[int]]) -> bool:
        """Check if two atoms are in the same molecule"""
        if not self.config['enable_molecule_protection']:
            return False
        
        return any(i in mol and j in mol for mol in molecules)
    
    def _has_pbc(self, atoms: Atoms) -> np.ndarray:
        """Check which directions have periodic boundaries"""
        cellpar = atoms.cell.cellpar()
        positions = atoms.get_positions()
        cell = atoms.get_cell()
        
        try:
            inv_cell = np.linalg.inv(cell.T)
            frac = np.dot(positions, inv_cell) % 1.0
            ranges = np.ptp(frac, axis=0)
            vacuum_sizes = cellpar[:3] * (1.0 - ranges)
            return vacuum_sizes < 3.0
        except np.linalg.LinAlgError:
            return np.array([True, True, True])
    
    def _log_molecules(self, molecule_analysis: List[Dict[str, Any]], idx: int) -> None:
        """Log detected molecules using unified analysis results"""
        logger = logging.getLogger('GPBB')
        
        # Separate single atoms and multi-atom molecules
        single_atoms = [m for m in molecule_analysis if m['n_atoms'] == 1]
        molecules = [m for m in molecule_analysis if m['n_atoms'] > 1]
        
        # Log molecules
        if molecules:
            logger.debug(f"Structure {idx}: {len(molecules)} molecule(s) detected")
            for i, mol in enumerate(molecules, 1):
                logger.debug(f"  Molecule {i}: {mol['formula']} ({mol['type']})")
        
        # Log single atoms
        if single_atoms:
            logger.debug(f"Structure {idx}: {len(single_atoms)} single atom adsorbate(s) detected")
            for atom in single_atoms:
                logger.debug(f"  Single atom: {atom['formula']} at index {atom['indices'][0]}")


def main():
    """Main entry point"""
    config_file = 'config.yaml'
    
    if not os.path.exists(config_file):
        print(f"Configuration file {config_file} not found")
        return
    
    try:
        algorithm = BLAdjustGPBB(config_file)
        algorithm.run()
    except Exception as e:
        print(f"Error: {e}")
        raise


if __name__ == '__main__':
    main()