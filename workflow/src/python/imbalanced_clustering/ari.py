from .utils import pair_confusion_matrix


def balanced_adjusted_rand_index(labels_true, labels_pred, reweigh=True):
    """Rand index adjusted for chance and balanced across true labels.
    The Rand Index computes a similarity measure between two clusterings
    by considering all pairs of samples and counting pairs that are
    assigned in the same or different clusters in the predicted and
    true clusterings.
    The raw RI score is then "adjusted for chance" into the ARI score
    using the following scheme::
        ARI = (RI - Expected_RI) / (max(RI) - Expected_RI)
    The adjusted Rand index is thus ensured to have a value close to
    0.0 for random labeling independently of the number of clusters and
    samples and exactly 1.0 when the clusterings are identical (up to
    a permutation).
    The original ARI is a symmetric measure:
        adjusted_rand_score(a, b) == adjusted_rand_score(b, a)
    But this does not hold due for the balanced ARI metric.
    The balanced ARI is obtained by reweighing the contingency table
    for all true label marginals, such that they sum to the same nummber,
    while preserving the total number of samples.
    Parameters
    ----------
    labels_true : int array, shape = [n_samples]
        Ground truth class labels to be used as a reference
    labels_pred : array-like of shape (n_samples,)
        Cluster labels to evaluate
    reweigh : bool, default=True
        if `True`, reweighs the contingency table based on the true labels
        such that they all have equal membership. The total number of samples
        is preserved with a round-off error. If 'False', this reverts the
        balanced ARI to the original ARI implementation.
    Returns
    -------
    ARI : float
       Similarity score between -1.0 and 1.0. Random labelings have an ARI
       close to 0.0. 1.0 stands for perfect match.
    References
    ----------
    .. [Hubert1985] L. Hubert and P. Arabie, Comparing Partitions,
      Journal of Classification 1985
      https://link.springer.com/article/10.1007%2FBF01908075
    .. [Steinley2004] D. Steinley, Properties of the Hubert-Arabie
      adjusted Rand index, Psychological Methods 2004
    .. [wk] https://en.wikipedia.org/wiki/Rand_index#Adjusted_Rand_index
    """
    (tn, fp), (fn, tp) = pair_confusion_matrix(
        labels_true, labels_pred, reweigh=reweigh
    )
    # convert to Python integer types, to avoid overflow or underflow
    tn, fp, fn, tp = int(tn), int(fp), int(fn), int(tp)

    # Special cases: empty data or full agreement
    if fn == 0 and fp == 0:
        return 1.0

    return 2.0 * (tp * tn - fn * fp) / ((tp + fn) * (fn + tn) + (tp + fp) * (fp + tn))
