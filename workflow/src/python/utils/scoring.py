import numpy as np
import pandas as pd
from sklearn import metrics

def return_scores(subset_1, subset_2):
    """Function that returns ARI, NMI, homogeneity, and completeness scores.

    Args:
        subset_1 (array): Array of values from subset 1.
        subset_2 (array): Array of values from subset 2.
    """
    ari_val = metrics.adjusted_rand_score(
        subset_1, subset_2
    )
    nmi_val = metrics.adjusted_mutual_info_score(
        subset_1, subset_2
    )
    homogeneity_val = metrics.homogeneity_score(
        subset_1, subset_2
    )
    completeness_val = metrics.completeness_score(
        subset_1, subset_2
    )
    results_dict = {
        "ARI": ari_val,
        "NMI": nmi_val,
        "Homogeneity": homogeneity_val,
        "Completeness": completeness_val
    }
    return results_dict
