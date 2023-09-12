from utils.utils import *
from components.peaks_pca import generate_pca_plot, generate_pca_plot2, generate_pca_plot3
from pathlib import Path
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

# Create a 2D DataFrame
df = pd.DataFrame({'x':[1,2,3,4,5], 'y':[4,6,5,8,2]})

# Visualize the DataFrame as a scatterplot using Seaborn
sns.scatterplot(x='x', y='y', data=df)

"""
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
"""

experiment = ExperimentalDesign({
"KSA041": "NCIH211_DMSO",
"KSA042": "NCIH211_DMSO",
"KSA043": "NCIH211_FHD286",
"KSA044": "NCIH211_FHD286",
"KSA045": "NCIH211_FHD286",
"KSA046": "CORL311_DMSO",
"KSA047": "CORL311_DMSO",
"KSA048": "CORL311_DMSO",
"KSA049": "CORL311_FHD286",
"KSA050": "CORL311_FHD286",
"KSA051": "CORL311_FHD286",
"KSA052": "NCIH526_DMSO",
"KSA053": "NCIH526_DMSO",
"KSA054": "NCIH526_DMSO",
"KSA055": "NCIH526_FHD286",
"KSA056": "NCIH526_FHD286",
"KSA057": "NCIH526_FHD286",
"KSA058": "NCIH1048_DMSO",
"KSA059": "NCIH1048_DMSO",
"KSA060": "NCIH1048_FHD286",
"KSA061": "NCIH1048_FHD286",
"KSA062": "NCIH1048_FHD286",
"KSA063": "NCIH211_DMSO",
"KSA064": "NCIH211_DMSO",
"KSA065": "NCIH211_FHD609",
"KSA066": "NCIH211_FHD609",
"KSA067": "NCIH211_FHD609",
"KSA068": "NCIH1048_DMSO",
"KSA069": "NCIH1048_DMSO",
"KSA070": "NCIH1048_FHD609",
"KSA071": "NCIH1048_FHD609",
"KSA072": "NCIH1048_FHD609"
})

generate_pca_plot(
    "/Users/btudorpr/PycharmProjects/kadoch_pipeline_tools/misc/raw.formatted.tsv",
    experiment)
plt.show()

