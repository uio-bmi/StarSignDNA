#!/usr/bin/env python3
"""
Plotting script for mutational signature analysis results.

This script creates a visualization of the hyperparameter optimization results,
showing the relationship between lambda (regularization) and prediction error
for different numbers of signatures (k).
"""

import matplotlib.pyplot as plt
import pandas as pd


def main():
    """Create hyperparameter optimization plot."""
    # Read the results data
    data = pd.read_csv(snakemake.input[0], delimiter="\t")
    data.columns = ['k', 'Lambda', 'logpmf']
    
    # Normalize logpmf values for better visualization
    max_logpmf = data['logpmf'].max()
    data['logpmf'] = data['logpmf'] / max_logpmf
    
    # Find the best performing parameters (minimum logpmf)
    smallest = data.loc[data['logpmf'].idxmin()]
    highlight_label = f"k = {smallest['k']} and Lambda = {smallest['Lambda']}"
    
    # Get unique values of k for plotting
    unique_k = data['k'].unique()
    
    # Create the plot
    plt.figure(figsize=(15, 10))
    
    # Plot lines for each k value
    for k in unique_k:
        k_data = data[data['k'] == k]
        plt.plot(k_data['Lambda'], k_data['logpmf'], label=f'k = {k}')
    
    # Highlight the best performing point
    plt.scatter(smallest['Lambda'], smallest['logpmf'], color='red', s=100, label='')
    
    # Add annotation for the best point
    plt.annotate(highlight_label, 
                 xy=(smallest['Lambda'], smallest['logpmf']), 
                 xytext=(smallest['Lambda'] + 0.025, smallest['logpmf'] - 0.04), 
                 va='bottom', ha='right', 
                 color='red', fontsize=12,
                 arrowprops=dict(facecolor='red', arrowstyle='->'))
    
    # Set plot properties
    plt.title('Hyperparameter Optimization: logpmf vs Lambda')
    plt.xlabel('Regularization parameter (Lambda)')
    plt.ylabel('Normalized prediction error (logpmf)')
    plt.legend(title='Number of signatures (k)', loc='upper left', bbox_to_anchor=(1.05, 1))
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(snakemake.output[0], dpi=600, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    main()

