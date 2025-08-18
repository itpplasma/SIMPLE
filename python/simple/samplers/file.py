"""
File-based particle initialization for SIMPLE Python API.

Supports loading particle initial conditions from various file formats
with compatibility for existing SIMPLE file formats and scientific data formats.
"""

import numpy as np
from typing import Optional, Dict, Any, Union
from pathlib import Path
import warnings

from ..core.batch import ParticleBatch


class FileSampler:
    """
    File-based particle initialization supporting multiple formats.
    
    Supports:
    - SIMPLE native format (start.dat style)
    - NumPy arrays (.npy, .npz)
    - HDF5 format
    - ASCII text files
    - NetCDF format (for VMEC compatibility)
    """
    
    def __init__(self):
        """Initialize file sampler"""
        pass
    
    def load_simple_format(
        self,
        filename: str,
        max_particles: Optional[int] = None
    ) -> ParticleBatch:
        """
        Load particles from SIMPLE native format (start.dat style).
        
        Expected format: 5 columns (s, theta, phi, v_par, mu) with one particle per row
        
        Args:
            filename: Path to SIMPLE format file
            max_particles: Maximum number of particles to load
            
        Returns:
            ParticleBatch: Loaded particle batch
        """
        file_path = Path(filename)
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {filename}")
        
        try:
            # Load data assuming space/tab separated columns
            data = np.loadtxt(filename)
            
            if data.ndim == 1:
                # Single particle
                data = data.reshape(1, -1)
            
            if data.shape[1] != 5:
                raise ValueError(f"Expected 5 columns, got {data.shape[1]} in {filename}")
            
            n_particles = data.shape[0]
            if max_particles is not None and n_particles > max_particles:
                data = data[:max_particles]
                n_particles = max_particles
                print(f"Loaded first {max_particles} particles from {filename}")
            
            # Create batch and copy data (transpose to SoA format)
            batch = ParticleBatch(n_particles)
            batch.positions[:] = data.T  # Transpose from AoS to SoA
            
            print(f"Loaded {n_particles} particles from {filename}")
            return batch
            
        except Exception as e:
            raise RuntimeError(f"Failed to load SIMPLE format file {filename}: {e}")
    
    def load_numpy(self, filename: str, array_name: Optional[str] = None) -> ParticleBatch:
        """
        Load particles from NumPy file format.
        
        Args:
            filename: Path to .npy or .npz file
            array_name: Array name for .npz files (uses first array if None)
            
        Returns:
            ParticleBatch: Loaded particle batch
        """
        file_path = Path(filename)
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {filename}")
        
        try:
            if filename.endswith('.npy'):
                data = np.load(filename)
            elif filename.endswith('.npz'):
                npz_data = np.load(filename)
                if array_name is not None:
                    if array_name not in npz_data:
                        raise KeyError(f"Array '{array_name}' not found in {filename}")
                    data = npz_data[array_name]
                else:
                    # Use first array
                    array_names = list(npz_data.keys())
                    if not array_names:
                        raise ValueError(f"No arrays found in {filename}")
                    data = npz_data[array_names[0]]
                    print(f"Using array '{array_names[0]}' from {filename}")
            else:
                raise ValueError(f"Unsupported file extension for {filename}")
            
            # Validate and reshape data
            if data.ndim != 2:
                raise ValueError(f"Expected 2D array, got {data.ndim}D")
            
            # Handle both SoA (5, n) and AoS (n, 5) formats
            if data.shape[0] == 5:
                # SoA format
                n_particles = data.shape[1]
                positions = data
            elif data.shape[1] == 5:
                # AoS format - transpose
                n_particles = data.shape[0]
                positions = data.T
            else:
                raise ValueError(f"Expected shape (5, n) or (n, 5), got {data.shape}")
            
            batch = ParticleBatch(n_particles)
            batch.positions[:] = positions
            
            print(f"Loaded {n_particles} particles from {filename}")
            return batch
            
        except Exception as e:
            raise RuntimeError(f"Failed to load NumPy file {filename}: {e}")
    
    def load_hdf5(
        self,
        filename: str,
        dataset_path: str = '/particles',
        group: Optional[str] = None
    ) -> ParticleBatch:
        """
        Load particles from HDF5 format.
        
        Args:
            filename: Path to HDF5 file
            dataset_path: Path to particle dataset within file
            group: HDF5 group name (if different from dataset_path)
            
        Returns:
            ParticleBatch: Loaded particle batch
        """
        try:
            import h5py
        except ImportError:
            raise RuntimeError("h5py required for HDF5 file support")
        
        file_path = Path(filename)
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {filename}")
        
        try:
            with h5py.File(filename, 'r') as f:
                if group is not None:
                    if group not in f:
                        raise KeyError(f"Group '{group}' not found in {filename}")
                    grp = f[group]
                    if dataset_path.startswith('/'):
                        dataset_path = dataset_path[1:]  # Remove leading slash
                    data = grp[dataset_path]
                else:
                    data = f[dataset_path]
                
                # Load data array
                positions = np.array(data)
                
                # Validate format
                if positions.ndim != 2:
                    raise ValueError(f"Expected 2D dataset, got {positions.ndim}D")
                
                if positions.shape[0] == 5:
                    n_particles = positions.shape[1]
                elif positions.shape[1] == 5:
                    n_particles = positions.shape[0]
                    positions = positions.T
                else:
                    raise ValueError(f"Expected shape (5, n) or (n, 5), got {positions.shape}")
                
                batch = ParticleBatch(n_particles)
                batch.positions[:] = positions
                
                print(f"Loaded {n_particles} particles from {filename}:{dataset_path}")
                return batch
                
        except Exception as e:
            raise RuntimeError(f"Failed to load HDF5 file {filename}: {e}")
    
    def load_ascii(
        self,
        filename: str,
        delimiter: Optional[str] = None,
        skiprows: int = 0,
        max_particles: Optional[int] = None,
        column_order: Optional[list] = None
    ) -> ParticleBatch:
        """
        Load particles from ASCII text file.
        
        Args:
            filename: Path to ASCII file
            delimiter: Column delimiter (auto-detect if None)
            skiprows: Number of header rows to skip
            max_particles: Maximum particles to load
            column_order: Order of columns if not [s, theta, phi, v_par, mu]
            
        Returns:
            ParticleBatch: Loaded particle batch
        """
        file_path = Path(filename)
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {filename}")
        
        try:
            # Load with automatic delimiter detection
            data = np.loadtxt(filename, delimiter=delimiter, skiprows=skiprows)
            
            if data.ndim == 1:
                data = data.reshape(1, -1)
            
            if data.shape[1] < 5:
                raise ValueError(f"Need at least 5 columns, got {data.shape[1]}")
            
            # Take only first 5 columns if more exist
            if data.shape[1] > 5:
                data = data[:, :5]
                print(f"Using first 5 columns from {data.shape[1]} available")
            
            # Apply column reordering if specified
            if column_order is not None:
                if len(column_order) != 5:
                    raise ValueError("column_order must have 5 elements")
                data = data[:, column_order]
            
            n_particles = data.shape[0]
            if max_particles is not None and n_particles > max_particles:
                data = data[:max_particles]
                n_particles = max_particles
            
            batch = ParticleBatch(n_particles)
            batch.positions[:] = data.T  # Transpose to SoA
            
            print(f"Loaded {n_particles} particles from {filename}")
            return batch
            
        except Exception as e:
            raise RuntimeError(f"Failed to load ASCII file {filename}: {e}")
    
    def save_simple_format(self, batch: ParticleBatch, filename: str):
        """
        Save particles to SIMPLE native format.
        
        Args:
            batch: ParticleBatch to save
            filename: Output file path
        """
        try:
            # Convert from SoA to AoS format for output
            data = batch.positions.T
            
            # Save with appropriate precision
            np.savetxt(filename, data, fmt='%.12e', 
                      header='s theta phi v_par mu', comments='# ')
            
            print(f"Saved {batch.n_particles} particles to {filename}")
            
        except Exception as e:
            raise RuntimeError(f"Failed to save to {filename}: {e}")
    
    def save_numpy(self, batch: ParticleBatch, filename: str, compressed: bool = True):
        """
        Save particles to NumPy format.
        
        Args:
            batch: ParticleBatch to save
            filename: Output file path (.npy or .npz)
            compressed: Use compression for .npz files
        """
        try:
            if filename.endswith('.npy'):
                np.save(filename, batch.positions)
            elif filename.endswith('.npz'):
                if compressed:
                    np.savez_compressed(filename, positions=batch.positions)
                else:
                    np.savez(filename, positions=batch.positions)
            else:
                # Default to .npz
                if not filename.endswith('.npz'):
                    filename += '.npz'
                if compressed:
                    np.savez_compressed(filename, positions=batch.positions)
                else:
                    np.savez(filename, positions=batch.positions)
            
            print(f"Saved {batch.n_particles} particles to {filename}")
            
        except Exception as e:
            raise RuntimeError(f"Failed to save to {filename}: {e}")
    
    def save_hdf5(
        self,
        batch: ParticleBatch,
        filename: str,
        dataset_path: str = '/particles',
        compression: str = 'gzip'
    ):
        """
        Save particles to HDF5 format.
        
        Args:
            batch: ParticleBatch to save
            filename: Output file path
            dataset_path: HDF5 dataset path
            compression: Compression method
        """
        try:
            import h5py
        except ImportError:
            raise RuntimeError("h5py required for HDF5 file support")
        
        try:
            with h5py.File(filename, 'w') as f:
                dset = f.create_dataset(
                    dataset_path,
                    data=batch.positions,
                    compression=compression,
                    shuffle=True
                )
                
                # Add metadata
                dset.attrs['n_particles'] = batch.n_particles
                dset.attrs['format'] = 'SoA'
                dset.attrs['coordinates'] = ['s', 'theta', 'phi', 'v_par', 'mu']
                
            print(f"Saved {batch.n_particles} particles to {filename}:{dataset_path}")
            
        except Exception as e:
            raise RuntimeError(f"Failed to save HDF5 file {filename}: {e}")
    
    def get_file_info(self, filename: str) -> Dict[str, Any]:
        """
        Get information about a particle file without loading it.
        
        Args:
            filename: Path to file
            
        Returns:
            Dict: File information
        """
        file_path = Path(filename)
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {filename}")
        
        info = {
            'filename': str(file_path),
            'size_bytes': file_path.stat().st_size,
            'format': 'unknown'
        }
        
        # Determine format and try to get basic info
        try:
            if filename.endswith('.npy'):
                info['format'] = 'numpy'
                data = np.load(filename, mmap_mode='r')  # Memory-mapped for large files
                info['shape'] = data.shape
                info['dtype'] = str(data.dtype)
                
            elif filename.endswith('.npz'):
                info['format'] = 'numpy_compressed'
                npz_data = np.load(filename)
                info['arrays'] = list(npz_data.keys())
                if 'positions' in npz_data:
                    info['shape'] = npz_data['positions'].shape
                
            elif filename.endswith('.h5') or filename.endswith('.hdf5'):
                info['format'] = 'hdf5'
                try:
                    import h5py
                    with h5py.File(filename, 'r') as f:
                        info['datasets'] = list(f.keys())
                        if '/particles' in f:
                            info['shape'] = f['/particles'].shape
                except ImportError:
                    info['error'] = 'h5py not available'
                    
            else:
                # Try as ASCII
                info['format'] = 'ascii'
                with open(filename, 'r') as f:
                    lines = f.readlines()
                    info['lines'] = len(lines)
                    if lines:
                        first_line = lines[0].strip()
                        if first_line and not first_line.startswith('#'):
                            info['columns'] = len(first_line.split())
                            
        except Exception as e:
            info['error'] = str(e)
        
        return info
    
    def validate_file_format(self, filename: str, expected_particles: Optional[int] = None) -> bool:
        """
        Validate that file contains properly formatted particle data.
        
        Args:
            filename: Path to file
            expected_particles: Expected number of particles (optional)
            
        Returns:
            bool: True if file format is valid
        """
        try:
            info = self.get_file_info(filename)
            
            if 'error' in info:
                print(f"File validation error: {info['error']}")
                return False
            
            # Check shape if available
            if 'shape' in info:
                shape = info['shape']
                if len(shape) != 2:
                    print(f"Invalid shape: expected 2D array, got {len(shape)}D")
                    return False
                
                if shape[0] != 5 and shape[1] != 5:
                    print(f"Invalid shape: expected one dimension to be 5, got {shape}")
                    return False
                
                if expected_particles is not None:
                    n_particles = shape[1] if shape[0] == 5 else shape[0]
                    if n_particles != expected_particles:
                        print(f"Particle count mismatch: expected {expected_particles}, got {n_particles}")
                        return False
            
            return True
            
        except Exception as e:
            print(f"Validation failed: {e}")
            return False