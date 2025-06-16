#!/usr/bin/env python3
"""
Script to assign random ancestries to samples based on probability distribution.

Usage: python assign_ancestry.py <sample_file>

The sample file should contain one sample name per line.
"""

import sys
import random
import argparse

def main():
    # Define ancestry counts and calculate probabilities
    ancestry_counts = {
        'afr': 84148,
        'amr': 79106,
        'eas': 10099,
        'eur': 234353,
        'mid': 1545,
        'sas': 5579
    }
    
    # Calculate total count
    total_count = sum(ancestry_counts.values())
    
    # Create weighted list for random selection
    # Each ancestry appears in the list proportional to its count
    weighted_ancestries = []
    for ancestry, count in ancestry_counts.items():
        weighted_ancestries.extend([ancestry] * count)
    
    # Alternative approach using random.choices (more memory efficient)
    ancestries = list(ancestry_counts.keys())
    weights = list(ancestry_counts.values())
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Assign random ancestries to samples')
    parser.add_argument('sample_file', help='File containing sample names, one per line')
    parser.add_argument('--seed', type=int, help='Random seed for reproducible results')
    
    args = parser.parse_args()
    
    # Set random seed if provided
    if args.seed:
        random.seed(args.seed)
    
    # Read sample file and assign ancestries
    try:
        with open(args.sample_file, 'r') as f:
            for line in f:
                sample_name = line.strip()
                if sample_name:  # Skip empty lines
                    # Randomly choose ancestry based on weights
                    chosen_ancestry = random.choices(ancestries, weights=weights, k=1)[0]
                    
                    # Output tab-separated line
                    print(f"{sample_name}\t.\t.\t.\t{chosen_ancestry}")
                    
    except FileNotFoundError:
        print(f"Error: Could not find file '{args.sample_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
