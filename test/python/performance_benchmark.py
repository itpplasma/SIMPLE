#!/usr/bin/env python3
"""
Performance benchmarking utilities for Python batch API validation

Provides infrastructure to validate <5% overhead requirement and measure
performance characteristics of the Python API vs direct Fortran execution.
"""

import time
import numpy as np
import statistics
from typing import Dict, List, Callable, Any, Optional
from dataclasses import dataclass
from pathlib import Path
import json


@dataclass
class BenchmarkResult:
    """Results from a performance benchmark"""
    name: str
    mean_time: float
    std_time: float
    min_time: float
    max_time: float
    n_samples: int
    overhead_percent: Optional[float] = None
    additional_metrics: Dict[str, Any] = None


class PerformanceBenchmarker:
    """High-precision performance benchmarking for Python API validation"""
    
    def __init__(self, warmup_runs: int = 3, benchmark_runs: int = 10):
        """
        Initialize benchmarker
        
        Args:
            warmup_runs: Number of warmup runs to stabilize performance
            benchmark_runs: Number of measurement runs for statistics
        """
        self.warmup_runs = warmup_runs
        self.benchmark_runs = benchmark_runs
        self.results: List[BenchmarkResult] = []
    
    def benchmark_function(self, func: Callable, name: str, *args, **kwargs) -> BenchmarkResult:
        """
        Benchmark a function with high precision timing
        
        Given: Function to benchmark with arguments
        When: Running multiple iterations with warmup
        Then: Returns statistical timing results
        """
        # Warmup runs
        for _ in range(self.warmup_runs):
            func(*args, **kwargs)
        
        # Measurement runs
        times = []
        for _ in range(self.benchmark_runs):
            start_time = time.perf_counter()
            result = func(*args, **kwargs)
            end_time = time.perf_counter()
            times.append(end_time - start_time)
        
        # Calculate statistics
        mean_time = statistics.mean(times)
        std_time = statistics.stdev(times) if len(times) > 1 else 0.0
        min_time = min(times)
        max_time = max(times)
        
        benchmark_result = BenchmarkResult(
            name=name,
            mean_time=mean_time,
            std_time=std_time,
            min_time=min_time,
            max_time=max_time,
            n_samples=len(times)
        )
        
        self.results.append(benchmark_result)
        return benchmark_result
    
    def compare_implementations(self, 
                              fortran_func: Callable, 
                              python_func: Callable,
                              name: str,
                              *args, **kwargs) -> Dict[str, BenchmarkResult]:
        """
        Compare Fortran vs Python implementation performance
        
        Given: Fortran and Python implementations of same functionality
        When: Benchmarking both with identical inputs
        Then: Returns overhead calculation and detailed comparison
        """
        # Benchmark Fortran implementation
        fortran_result = self.benchmark_function(
            fortran_func, f"{name}_fortran", *args, **kwargs)
        
        # Benchmark Python implementation  
        python_result = self.benchmark_function(
            python_func, f"{name}_python", *args, **kwargs)
        
        # Calculate overhead percentage
        overhead_percent = ((python_result.mean_time - fortran_result.mean_time) 
                          / fortran_result.mean_time * 100)
        
        python_result.overhead_percent = overhead_percent
        
        return {
            'fortran': fortran_result,
            'python': python_result,
            'overhead_percent': overhead_percent
        }
    
    def validate_overhead_requirement(self, overhead_percent: float, 
                                    max_allowed: float = 5.0) -> bool:
        """
        Validate that overhead meets <5% requirement
        
        Given: Measured overhead percentage
        When: Comparing against requirement threshold
        Then: Returns pass/fail validation
        """
        return overhead_percent < max_allowed
    
    def export_results(self, filepath: Path) -> None:
        """Export benchmark results to JSON for analysis"""
        results_data = []
        for result in self.results:
            result_dict = {
                'name': result.name,
                'mean_time': result.mean_time,
                'std_time': result.std_time,
                'min_time': result.min_time,
                'max_time': result.max_time,
                'n_samples': result.n_samples,
                'overhead_percent': result.overhead_percent,
                'additional_metrics': result.additional_metrics or {}
            }
            results_data.append(result_dict)
        
        with open(filepath, 'w') as f:
            json.dump(results_data, f, indent=2)
    
    def print_summary(self) -> None:
        """Print human-readable benchmark summary"""
        print("\n" + "="*80)
        print("PERFORMANCE BENCHMARK SUMMARY")
        print("="*80)
        
        for result in self.results:
            print(f"\n{result.name}:")
            print(f"  Mean time: {result.mean_time:.6f}s ± {result.std_time:.6f}s")
            print(f"  Range: {result.min_time:.6f}s - {result.max_time:.6f}s")
            print(f"  Samples: {result.n_samples}")
            
            if result.overhead_percent is not None:
                status = "✓ PASS" if result.overhead_percent < 5.0 else "✗ FAIL"
                print(f"  Overhead: {result.overhead_percent:.2f}% {status}")


class ParticleBatchBenchmark:
    """Specialized benchmarking for particle batch operations"""
    
    def __init__(self, benchmarker: PerformanceBenchmarker):
        self.benchmarker = benchmarker
    
    def benchmark_particle_initialization(self, n_particles: int) -> BenchmarkResult:
        """
        Benchmark particle batch initialization
        
        Given: Number of particles to initialize
        When: Creating SoA particle arrays
        Then: Returns initialization performance metrics
        """
        def init_particles():
            # SoA initialization (target layout)
            particles = np.zeros((5, n_particles), dtype=np.float64)
            
            # Fill with realistic particle data
            particles[0, :] = np.random.uniform(0.1, 0.9, n_particles)  # s
            particles[1, :] = np.random.uniform(0, 2*np.pi, n_particles)  # theta
            particles[2, :] = np.random.uniform(0, 2*np.pi, n_particles)  # phi
            particles[3, :] = 1.0  # v/v0
            particles[4, :] = np.random.uniform(-1, 1, n_particles)  # v_par/v
            
            return particles
        
        return self.benchmarker.benchmark_function(
            init_particles, f"particle_init_{n_particles}")
    
    def benchmark_memory_access_patterns(self, particles: np.ndarray) -> Dict[str, BenchmarkResult]:
        """
        Benchmark different memory access patterns for SoA validation
        
        Given: SoA particle array
        When: Testing coordinate access patterns
        Then: Returns performance for different access strategies
        """
        n_particles = particles.shape[1]
        
        # Sequential coordinate access (optimal for SoA)
        def access_coordinates_sequential():
            s_coords = particles[0, :]  # All s coordinates
            theta_coords = particles[1, :]  # All theta coordinates
            return s_coords, theta_coords
        
        # Particle-wise access (sub-optimal for SoA)
        def access_particles_individual():
            results = []
            for i in range(min(1000, n_particles)):  # Limit for performance
                particle = particles[:, i]  # Single particle all coordinates
                results.append(particle)
            return results
        
        # Batch computation on coordinates
        def compute_batch_operations():
            # Example: compute cylindrical R coordinates
            s = particles[0, :]
            theta = particles[1, :]
            R = s * np.cos(theta)  # Vectorized operation
            return R
        
        return {
            'sequential': self.benchmarker.benchmark_function(
                access_coordinates_sequential, f"access_sequential_{n_particles}"),
            'individual': self.benchmarker.benchmark_function(
                access_particles_individual, f"access_individual_{n_particles}"),
            'batch_compute': self.benchmarker.benchmark_function(
                compute_batch_operations, f"batch_compute_{n_particles}")
        }
    
    def benchmark_scalability(self, particle_counts: List[int]) -> Dict[int, BenchmarkResult]:
        """
        Benchmark scalability across different particle counts
        
        Given: List of particle counts to test
        When: Running initialization benchmark for each count
        Then: Returns scalability profile
        """
        results = {}
        for count in particle_counts:
            result = self.benchmark_particle_initialization(count)
            results[count] = result
            
            # Calculate throughput (particles per second)
            throughput = count / result.mean_time
            result.additional_metrics = {'throughput_particles_per_sec': throughput}
        
        return results


def create_performance_report(benchmarker: PerformanceBenchmarker, 
                            output_dir: Path) -> None:
    """
    Create comprehensive performance report
    
    Given: Benchmark results
    When: Generating analysis report
    Then: Exports detailed performance analysis
    """
    output_dir.mkdir(exist_ok=True)
    
    # Export raw results
    benchmarker.export_results(output_dir / "benchmark_results.json")
    
    # Create summary report
    report_file = output_dir / "performance_report.txt"
    with open(report_file, 'w') as f:
        f.write("SIMPLE Python API Performance Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("OVERHEAD VALIDATION:\n")
        f.write("-" * 20 + "\n")
        
        for result in benchmarker.results:
            if result.overhead_percent is not None:
                status = "PASS" if result.overhead_percent < 5.0 else "FAIL"
                f.write(f"{result.name}: {result.overhead_percent:.2f}% [{status}]\n")
        
        f.write("\nDETAILED TIMING RESULTS:\n")
        f.write("-" * 25 + "\n")
        
        for result in benchmarker.results:
            f.write(f"\n{result.name}:\n")
            f.write(f"  Mean: {result.mean_time:.6f}s\n")
            f.write(f"  Std:  {result.std_time:.6f}s\n")
            f.write(f"  Min:  {result.min_time:.6f}s\n")
            f.write(f"  Max:  {result.max_time:.6f}s\n")
            
            if result.additional_metrics:
                f.write("  Additional metrics:\n")
                for key, value in result.additional_metrics.items():
                    f.write(f"    {key}: {value}\n")
    
    print(f"Performance report saved to: {report_file}")


if __name__ == "__main__":
    # Example usage and testing
    benchmarker = PerformanceBenchmarker()
    batch_benchmark = ParticleBatchBenchmark(benchmarker)
    
    # Test particle initialization at different scales
    particle_counts = [1000, 10000, 100000]
    print("Testing particle initialization scalability...")
    
    scalability_results = batch_benchmark.benchmark_scalability(particle_counts)
    
    # Test memory access patterns
    test_particles = np.random.random((5, 10000))
    print("Testing memory access patterns...")
    
    access_results = batch_benchmark.benchmark_memory_access_patterns(test_particles)
    
    # Print summary
    benchmarker.print_summary()
    
    # Create report
    output_dir = Path("performance_results")
    create_performance_report(benchmarker, output_dir)