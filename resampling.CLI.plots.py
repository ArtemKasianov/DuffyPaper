import numpy as np
import argparse
import matplotlib.pyplot as plt

def resample_allele_frequencies(initial_freq, num_iterations, sample_size, obs_freq):
    """
    Simulate resampling of alleles given an initial allele C frequency.

    Returns:
        List of float: List of allele C frequencies from each iteration.
    """
    frequencies_C = []
    itresh = 10
    for i in range(num_iterations):
        rng = np.random.default_rng()
        count_C = rng.binomial(n=sample_size, p=initial_freq)
        freq_C = count_C / sample_size
        frequencies_C.append(freq_C)
        
        if (i == itresh):
          plot_histogram(frequencies_C, obs_freq, i)
          itresh = itresh*10
    plot_histogram(frequencies_C, obs_freq, num_iterations)
    return frequencies_C

def plot_histogram(frequencies_C, obs_freq, num_iterations):
    plt.figure(figsize=(8, 6))
    plt.hist(frequencies_C, bins=30, alpha=0.7, color='blue', edgecolor='black')
    plt.axvline(obs_freq, color='red', linestyle='dashed', linewidth=2)
    plt.xlabel("rs2814778 *C frequency", fontsize = 16)
    plt.ylabel("Count", fontsize = 16)
    plt.tick_params(axis='both', which='major', labelsize=11)
    
    # Save the histogram to a file
    filename = f"histogram_{num_iterations}_iterations.tiff"
    plt.savefig(filename)
    filename = f"histogram_{num_iterations}_iterations.svg"
    plt.savefig(filename)
    filename = f"histogram_{num_iterations}_iterations.pdf"
    plt.savefig(filename)
    print(f"Histogram saved as {filename}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Simulate allele frequency resampling.")
    parser.add_argument("--allele_c_frequency", type=float, required=True, help="Frequency of allele C (0 to 1).")
    parser.add_argument("--num_resamplings", type=int, required=True, help="Number of resampling iterations (should be a power of ten).")
    parser.add_argument("--sample_size", type=int, required=True, help="Sample size for each resampling.")
    parser.add_argument("--obs_allele_c_frequency", type=float, required=True, help="Observed frequency of allele C (0 to 1).")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not (0 <= args.allele_c_frequency <= 1):
        raise ValueError("Frequency of allele C must be between 0 and 1.")
    if args.num_resamplings <= 0:
        raise ValueError("Number of resamplings must be a positive integer.")
    if args.sample_size <= 0:
        raise ValueError("Sample size must be a positive integer.")
    if not (0 <= args.obs_allele_c_frequency <= 1):
        raise ValueError("Observed frequency of allele C must be between 0 and 1.")
    
    # Run the resampling simulation
    frequencies_C = resample_allele_frequencies(
        args.allele_c_frequency, args.num_resamplings, args.sample_size, args.obs_allele_c_frequency
    )
    
    # Plot and save the histogram
    
    
if __name__ == '__main__':
    main()