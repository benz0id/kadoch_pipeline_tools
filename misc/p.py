from utils.utils import *
from components.peaks_pca import generate_pca_plot
from pathlib import Path
import seaborn as sns
import pandas as pd

# Create a 2D DataFrame
df = pd.DataFrame({'x':[1,2,3,4,5], 'y':[4,6,5,8,2]})

# Visualize the DataFrame as a scatterplot using Seaborn
sns.scatterplot(x='x', y='y', data=df)

experiment = ExperimentalDesign(
    {
        'GXATAC001': 'empty',
        'GXATAC002': 'empty',
        'GXATAC003': 'D81WT',
        'GXATAC004': 'D81WT',
        'GXATAC005': 'D81A',
        'GXATAC006': 'D81A',
        'GXATAC007': 'D81N',
        'GXATAC008': 'D81N'
    })

generate_pca_plot(
    "/Users/btudorpr/PycharmProjects/kadoch_pipeline_tools/test_project/counts_matrix.tsv",
    experiment)