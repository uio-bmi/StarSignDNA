#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
snakemake = snakemake 

def main():
    data = pd.read_csv(snakemake.input[0],delimiter="\t") 
    data.columns = ['k', 'Lambda', 'logpmf']
    max_logpmf = data['logpmf'].max()
    data['logpmf'] = data['logpmf'] / max_logpmf
    smallest = data.loc[data['logpmf'].idxmin()]
    highlight_label = f"k = {smallest['k']} and Lambda = {smallest['Lambda']}"
    unique_k = data['k'].unique()

    plt.figure(figsize=(15, 10))
    for k in unique_k:
        k_data = data[data['k'] == k]
        plt.plot(k_data['Lambda'], k_data['logpmf'], label=f'k = {k}')

    plt.scatter(smallest['Lambda'], smallest['logpmf'], color='red', s=100, label='')
    
    # Annotate with the highlight_label at the location of the smallest point
    plt.annotate(highlight_label, 
                 xy=(smallest['Lambda'], smallest['logpmf']), 
                 xytext=(smallest['Lambda'] + 0.025, smallest['logpmf'] - 0.04), 
                 va='bottom', ha='right', 
                 color='red', fontsize=12,
                 arrowprops=dict(facecolor='red', arrowstyle='->'))
    
    plt.title('Plot of logpmf vs Lambda')
    plt.xlabel('Hyperparameter lambda')
    plt.ylabel('Prediction error (logpmf)')
    plt.legend(title='Number of signatures (k)', loc='upper left', bbox_to_anchor=(1.05, 1))
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(snakemake.output[0], dpi=600)

if __name__ == "__main__":
    main()

