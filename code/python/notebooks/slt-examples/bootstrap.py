
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_pairwise_scatter(p1_samples, p2_samples, rho1_samples, rho2_samples, 
                         figsize=(10, 10)):
    """
    Create pairwise scatter plots using seaborn.
    
    Parameters:
    -----------
    p1_samples, p2_samples, rho1_samples, rho2_samples : numpy arrays
        Bootstrap samples for each parameter
    figsize : tuple
        Figure size
    """
    
    # Create DataFrame
    data = pd.DataFrame({
        'ρ₁': rho1_samples,
        'ρ₂': rho2_samples, 
        'p₁': p1_samples,
        'p₂': p2_samples
    })
    
    # Create pairplot
    g = sns.pairplot(data, 
                     diag_kind='hist',      # Histograms on diagonal
                     plot_kws={'alpha': 0.6, 's': 20},  # Scatter plot settings
                     diag_kws={'alpha': 0.7})           # Histogram settings
    
    # Customize the plot
    g.fig.suptitle('Bootstrap Samples: Pairwise Scatter Plots', 
                   fontsize=16, y=1.02)
    
    return g