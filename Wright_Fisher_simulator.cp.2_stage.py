import numpy as np
import matplotlib.pyplot as plt
import argparse
import random
import math

def wright_fisher_simulation(N, p0, generations1, generations2, sample_size_val, seed=None):
    """
    Perform two stage Wright-Fisher sampling simulation with continuous gene flow.

    Returns:
        Final frequency of allele A.
    """
    if seed is not None:
        np.random.seed(seed)
        
    # First stage    
    p = p0
    num_A = 0
    for gen in range(1, generations1 + 1):
        num_A = np.random.binomial(N, p)
        p = num_A / N
        outer_num_A = num_A
        if num_A == 0:
            break
    
    
    array = ['A'] * num_A + ['B'] * (N - num_A)
    
    sample = random.sample(array, sample_size_val)
    
    # Calculate frequency of 'A' in the sample
    freq_A_interm = sample.count('A') / sample_size_val
    
    # Second stage    
    p = 1-p
    w = -math.log(0.42)/generations2
    outer_num_A = 0
    countOfNotDuffy = int(N*(1-p))
    array = ['B'] * countOfNotDuffy + ['A'] * (N-countOfNotDuffy)
    for gen in range(1, generations2 + 1):
        num_B = int(N*(1-w))
        num_A = N - num_B
        
        sample_B = random.sample(array, num_B)
        
        sample_A = ['A']*num_A
        
        array = sample_A + sample_B
        
        p = array.count('A') / N
       
        num_A = np.random.binomial(N, p)
        array = ['A'] * num_A + ['B'] * (N - num_A)
        
        
    array = ['B'] * num_A + ['A'] * (N - num_A)
    
    sample = random.sample(array, sample_size_val)
    
    # Calculate frequency of 'A' in the sample
    freq_A = sample.count('A') / sample_size_val
    

    return (freq_A_interm,freq_A)

def main():
    parser = argparse.ArgumentParser(description='Wright-Fisher Sampling Simulation (Multiple Iterations)')
    parser.add_argument('--N', type=int, required=True,
                        help='Population size (number of individuals)')
    parser.add_argument('--p0', type=float, required=True,
                        help='Initial frequency of allele A (between 0 and 1)')
    parser.add_argument('--generations1', type=int, required=True,
                        help='Number of generations to simulate')
    parser.add_argument('--generations2', type=int, required=True,
                        help='Number of generations to simulate')
    parser.add_argument('--iterations', type=int, required=True,
                        help='Number of independent simulations (replicates)')
    parser.add_argument('--sample_size', type=int, required=True,
                        help='Number of final allele frequencies to sample for plotting')
    parser.add_argument('--seed', type=int, default=None,
                        help='Base random seed (optional)')
    parser.add_argument('--freq_to_check', type=float, default=None,
                        help='Frequency to check')
    parser.add_argument('--output', type=str, default=None,
                        help='Optional path to save plot image (e.g. distribution.png)')
    parser.add_argument('--outputDistrVals', type=str, default=None,
                        help='Optional path to distr values file ')
    parser.add_argument('--outputDistrIntermVals', type=str, default=None,
                        help='Optional path to distr values file ')
    
    args = parser.parse_args()

    # Perform multiple iterations
    final_freqs = []
    final_freqs_interm = []
    
    checkFreqs = float(args.freq_to_check)
    
    base_seed = args.seed if args.seed is not None else np.random.randint(0, 1e9)
    count_freqs_greater = 0
    for i in range(args.iterations):
        # Use different seed for each replicate
        sim_seed = base_seed + i
        (interm_p,final_p) = wright_fisher_simulation(args.N, args.p0, args.generations1,args.generations2, int(args.sample_size), sim_seed)
        final_freqs.append(final_p)
        final_freqs_interm.append(interm_p)
        if final_p >= checkFreqs:
            count_freqs_greater = count_freqs_greater + 1
    
    
    
    with open(args.outputDistrVals, "w") as f:
        for curr_final_freq in final_freqs:
            f.write(f"{curr_final_freq}\n")
    with open(args.outputDistrIntermVals, "w") as f:
        for curr_final_freq in final_freqs_interm:
            f.write(f"{curr_final_freq}\n")
    
    
    # Plot the distribution
    
    print("count_freqs_greater - " + str(count_freqs_greater))
    plt.figure(figsize=(10, 6))
    plt.hist(final_freqs, bins=20, edgecolor='black', density=True)
    plt.xlabel('Final frequency of allele A')
    plt.ylabel('Density')
    plt.title(f'Distribution of final allele A frequencies\n'
              f'(N={args.N}, p0={args.p0}, generations1={args.generations1}, generations2={args.generations2}, iterations={args.iterations})')
    plt.grid(True)

    if args.output:
        plt.savefig(args.output)
        print(f'Plot saved to {args.output}')
    else:
        plt.show()

if __name__ == '__main__':
    main()