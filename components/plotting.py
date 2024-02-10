import matplotlib.pyplot as plt
import pandas as pd
from pandas.core.frame import DataFrame
from typing import Any, Union


def horizontal_stacked_barplot(df: DataFrame, title: str, count_col: str,
                               names_col: str, colour_col: str) -> None:
    """
    Create a horizontal stacked bar plot.

    Parameters:
        df (DataFrame): The pandas DataFrame containing the data.
        count_col (str): The name of the column containing counts.
        names_col (str): The name of the column containing category names.
        colour_col (str): The name of the column containing colors.

    Returns:
        None
    """
    try:
        # Check if columns exist in DataFrame
        if count_col not in df.columns or names_col not in df.columns or colour_col not in df.columns:
            raise ValueError(
                "Invalid column names. Please check your DataFrame columns.")

        # Create horizontal bar plot
        fig, ax = plt.subplots(figsize=(8, 6))
        bars = ax.barh(df[names_col], df[count_col],
                       color=df[colour_col])

        # Add labels and legend
        ax.set_xlabel('Count', fontsize=12)
        ax.set_ylabel('Names', fontsize=12)
        ax.set_title(title, fontsize=14)
        ax.legend(bars, df[colour_col], fontsize=10)

        # Show plot
        plt.show()

    except Exception as e:
        print(f"An error occurred: {str(e)}")