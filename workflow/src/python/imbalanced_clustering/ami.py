import numpy as np

from .utils import (
    check_clusterings,
    contingency_matrix,
    entropy,
    mutual_info_score,
    expected_mutual_information,
    generalized_average,
)


def balanced_adjusted_mutual_info(
    labels_true, labels_pred, *, average_method="arithmetic", reweigh=True
):
    """Mutual Information adjusted for chance and balanced across true labels.
    Adjusted Mutual Information (AMI) is an adjustment of the Mutual
    Information (MI) score to account for chance. It accounts for the fact that
    the MI is generally higher for two clusterings with a larger number of
    clusters, regardless of whether there is actually more information shared.
    For two clusterings :math:`U` and :math:`V`, the AMI is given as::
        AMI(U, V) = [MI(U, V) - E(MI(U, V))] / [avg(H(U), H(V)) - E(MI(U, V))]
    This metric is independent of the absolute values of the labels:
    a permutation of the class or cluster label values won't change the
    score value in any way.
    The original AMI metric is a symmetric measure: switching :math:`U`
    (``label_true``) with :math:`V` (``labels_pred``) will return the same score value,
    but this is not the case for the reweighted and balanced AMI.
    The balanced AMI is obtained by reweighing the contingency table
    for all true label marginals, such that they sum to the same nummber,
    while preserving the total number of samples.
    Be mindful that this function is an order of magnitude slower than other
    metrics, such as the Adjusted Rand Index.
    Parameters
    ----------
    labels_true : int array, shape = [n_samples]
        A clustering of the data into disjoint subsets, called :math:`U` in
        the above formula.
    labels_pred : int array-like of shape (n_samples,)
        A clustering of the data into disjoint subsets, called :math:`V` in
        the above formula.
    average_method : str, default='arithmetic'
        How to compute the normalizer in the denominator. Possible options
        are 'min', 'geometric', 'arithmetic', and 'max'.
        .. versionadded:: 0.20
        .. versionchanged:: 0.22
           The default value of ``average_method`` changed from 'max' to
           'arithmetic'.
    reweigh : bool, default=True
        if `True`, reweighs the contingency table based on the true labels
        such that they all have equal membership. The total number of samples
        is preserved with a round-off error. If 'False', this reverts the
        balanced AMI to the original AMI implementation.
    Returns
    -------
    AMI: float (upperlimited by 1.0)
       The AMI returns a value of 1 when the two partitions are identical
       (ie perfectly matched). Random partitions (independent labellings) have
       an expected AMI around 0 on average hence can be negative. The value is
       in adjusted nats (based on the natural logarithm).
    References
    ----------
    .. [1] `Vinh, Epps, and Bailey, (2010). Information Theoretic Measures for
       Clusterings Comparison: Variants, Properties, Normalization and
       Correction for Chance, JMLR
       <http://jmlr.csail.mit.edu/papers/volume11/vinh10a/vinh10a.pdf>`_
    .. [2] `Wikipedia entry for the Adjusted Mutual Information
       <https://en.wikipedia.org/wiki/Adjusted_Mutual_Information>`_
    """
    labels_true, labels_pred = check_clusterings(labels_true, labels_pred)
    n_samples = labels_true.shape[0]
    classes = np.unique(labels_true)
    clusters = np.unique(labels_pred)
    # Special limit cases: no clustering since the data is not split.
    # This is a perfect match hence return 1.0.
    if (
        classes.shape[0] == clusters.shape[0] == 1
        or classes.shape[0] == clusters.shape[0] == 0
    ):
        return 1.0
    contingency = contingency_matrix(
        labels_true, labels_pred, reweigh=reweigh, sparse=True
    )
    contingency = contingency.astype(np.float64)
    # Calculate the MI for the two clusterings
    mi = mutual_info_score(labels_true, labels_pred, contingency=contingency)
    # Calculate the expected value for the mutual information
    emi = expected_mutual_information(contingency, n_samples)
    # Calculate entropy for each labeling
    h_true, h_pred = entropy(labels_true), entropy(labels_pred)
    normalizer = generalized_average(h_true, h_pred, average_method)
    denominator = normalizer - emi
    # Avoid 0.0 / 0.0 when expectation equals maximum, i.e a perfect match.
    # normalizer should always be >= emi, but because of floating-point
    # representation, sometimes emi is slightly larger. Correct this
    # by preserving the sign.
    if denominator < 0:
        denominator = min(denominator, -np.finfo("float64").eps)
    else:
        denominator = max(denominator, np.finfo("float64").eps)
    ami = (mi - emi) / denominator
    return ami
