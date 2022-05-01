import numpy as np
import scipy.sparse as sp

from .utils import entropy, mutual_info_score, contingency_matrix, check_clusterings


def balanced_homogeneity_completeness_v_measure(
    labels_true, labels_pred, *, beta=1.0, reweigh=True
):
    """Compute balanced the homogeneity, completeness, and V-Measure scores.
    Those metrics are based on normalized conditional entropy measures of
    the clustering labeling to evaluate given the knowledge of a Ground
    Truth class labels of the same samples.
    A clustering result satisfies homogeneity if all of its clusters
    contain only data points which are members of a single class.
    A clustering result satisfies completeness if all the data points
    that are members of a given class are elements of the same cluster.
    Both scores have positive values between 0.0 and 1.0, larger values
    being desirable.
    Those 3 metrics are independent of the absolute values of the labels:
    a permutation of the class or cluster label values won't change the
    score values in any way.
    The imbalanced V-Measure is symmetric: swapping ``labels_true`` and
    ``label_pred`` will give the same score, but this is not the case for the
    balanced V-Measure, as the true labels are reweighted. The symmetric
    property does not hold in general for balanced or imbalanced homogeneity and
    completeness measures. V-Measure is identical to
    :func:`normalized_mutual_info_score` with the arithmetic averaging method.
    The balanced homogeneity, completeness, and v-measure values are obtained by
    reweighing the contingency table for all true label marginals, such that they
    sum to the same nummber, while preserving the total number of samples.
    Parameters
    ----------
    labels_true : int array, shape = [n_samples]
        ground truth class labels to be used as a reference
    labels_pred : array-like of shape (n_samples,)
        cluster labels to evaluate
    beta : float, default=1.0
        Ratio of weight attributed to ``homogeneity`` vs ``completeness``.
        If ``beta`` is greater than 1, ``completeness`` is weighted more
        strongly in the calculation. If ``beta`` is less than 1,
        ``homogeneity`` is weighted more strongly.
    reweigh : bool, default=True
        if `True`, reweighs the contingency table based on the true labels
        such that they all have equal membership. The total number of samples
        is preserved with a round-off error. If 'False', this reverts the
        balanced to the original implementation for all measures.
    Returns
    -------
    balanced_homogeneity : float
       score between 0.0 and 1.0. 1.0 stands for perfectly homogeneous labeling
    balanced_completeness : float
       score between 0.0 and 1.0. 1.0 stands for perfectly complete labeling
    balanced_v_measure : float
        harmonic mean of the first two
    """
    labels_true, labels_pred = check_clusterings(labels_true, labels_pred)

    if len(labels_true) == 0:
        return 1.0, 1.0, 1.0
    
    contingency = contingency_matrix(
        labels_true, labels_pred, reweigh=reweigh, sparse=True
    )
    MI = mutual_info_score(None, None, contingency=contingency)

    # Recalculate labels_true and labels_pred if reweigh is True to
    # factor in the reweighting based on the true class frequencies.
    # These won't preserve order but this is fine since entropy is 
    # invariant to order
    if reweigh is True:
        true_sums = np.squeeze(np.asarray(sp.csc_matrix.sum(contingency, axis = 1)))
        pred_sums = np.squeeze(np.asarray(sp.csc_matrix.sum(contingency, axis = 0)))
        labels_true = np.repeat(
            np.arange(len(true_sums)), true_sums
        )
        labels_pred = np.repeat(
            np.arange(len(pred_sums)), pred_sums
        )

    entropy_C = entropy(labels_true)
    entropy_K = entropy(labels_pred)

    homogeneity = MI / (entropy_C) if entropy_C else 1.0
    completeness = MI / (entropy_K) if entropy_K else 1.0

    if homogeneity + completeness == 0.0:
        v_measure_score = 0.0
    else:
        v_measure_score = (
            (1 + beta)
            * homogeneity
            * completeness
            / (beta * homogeneity + completeness)
        )

    return homogeneity, completeness, v_measure_score


def balanced_homogeneity(labels_true, labels_pred, reweigh=True):
    """Class balanced homogeneity metric of a cluster labeling given a ground truth.
    A clustering result satisfies homogeneity if all of its clusters
    contain only data points which are members of a single class.
    This metric is independent of the absolute values of the labels:
    a permutation of the class or cluster label values won't change the
    score value in any way.
    This metric is not symmetric: switching ``label_true`` with ``label_pred``
    will return the :func:`completeness_score` which will be different in
    general.
    The balanced homogeneity is obtained by reweighing the contingency table for
    all true label marginals, such that they sum to the same nummber, while
    preserving the total number of samples.
    Parameters
    ----------
    labels_true : int array, shape = [n_samples]
        ground truth class labels to be used as a reference
    labels_pred : array-like of shape (n_samples,)
        cluster labels to evaluate
    reweigh : bool, default=True
        if `True`, reweighs the contingency table based on the true labels
        such that they all have equal membership. The total number of samples
        is preserved with a round-off error. If 'False', this reverts the
        balanced homogeneity to the original homogeneity implementation.
    Returns
    -------
    balanced_homogeneity : float
       score between 0.0 and 1.0. 1.0 stands for perfectly homogeneous labeling
    References
    ----------
    .. [1] `Andrew Rosenberg and Julia Hirschberg, 2007. V-Measure: A
       conditional entropy-based external cluster evaluation measure
       <https://aclweb.org/anthology/D/D07/D07-1043.pdf>`_
    """
    return balanced_homogeneity_completeness_v_measure(
        labels_true, labels_pred, reweigh=reweigh
    )[0]


def balanced_completeness(labels_true, labels_pred, reweigh=True):
    """Class balanced Completeness metric of a cluster labeling given a ground truth.
    A clustering result satisfies completeness if all the data points
    that are members of a given class are elements of the same cluster.
    This metric is independent of the absolute values of the labels:
    a permutation of the class or cluster label values won't change the
    score value in any way.
    This metric is not symmetric: switching ``label_true`` with ``label_pred``
    will return the :func:`homogeneity_score` which will be different in
    general.
    The balanced completeness is obtained by reweighing the contingency table for
    all true label marginals, such that they sum to the same nummber, while
    preserving the total number of samples.
    Parameters
    ----------
    labels_true : int array, shape = [n_samples]
        ground truth class labels to be used as a reference
    labels_pred : array-like of shape (n_samples,)
        cluster labels to evaluate
    reweigh : bool, default=True
            if `True`, reweighs the contingency table based on the true labels
            such that they all have equal membership. The total number of samples
            is preserved with a round-off error. If 'False', this reverts the
            balanced completeness to the original completeness implementation.
    Returns
    -------
    balanced_completeness : float
       score between 0.0 and 1.0. 1.0 stands for perfectly complete labeling
    References
    ----------
    .. [1] `Andrew Rosenberg and Julia Hirschberg, 2007. V-Measure: A
       conditional entropy-based external cluster evaluation measure
       <https://aclweb.org/anthology/D/D07/D07-1043.pdf>`_
    """
    return balanced_homogeneity_completeness_v_measure(
        labels_true, labels_pred, reweigh=reweigh
    )[1]


def balanced_v_measure(labels_true, labels_pred, *, beta=1.0, reweigh=True):
    """Class balanced V-measure cluster labeling given a ground truth.
    The V-measure is the harmonic mean between homogeneity and completeness::
        v = (1 + beta) * homogeneity * completeness
             / (beta * homogeneity + completeness)
    This metric is independent of the absolute values of the labels:
    a permutation of the class or cluster label values won't change the
    score value in any way.
    The imbalanced version of this metric is symmetric: switching ``label_true``
    with ``label_pred`` will return the same score value. However, this is not the
    case for the balanaced V-measure, as the true labels are reweighted.
    The balanced V-measure is obtained by reweighing the contingency table for
    all true label marginals, such that they sum to the same nummber, while
    preserving the total number of samples.
    Parameters
    ----------
    labels_true : int array, shape = [n_samples]
        ground truth class labels to be used as a reference
    labels_pred : array-like of shape (n_samples,)
        cluster labels to evaluate
    beta : float, default=1.0
        Ratio of weight attributed to ``homogeneity`` vs ``completeness``.
        If ``beta`` is greater than 1, ``completeness`` is weighted more
        strongly in the calculation. If ``beta`` is less than 1,
        ``homogeneity`` is weighted more strongly.
    reweigh : bool, default=True
            if `True`, reweighs the contingency table based on the true labels
            such that they all have equal membership. The total number of samples
            is preserved with a round-off error. If 'False', this reverts the
            balanced V-measure to the original V-measure implementation.
    Returns
    -------
    balanced_v_measure : float
       score between 0.0 and 1.0. 1.0 stands for perfectly complete labeling
    References
    ----------
    .. [1] `Andrew Rosenberg and Julia Hirschberg, 2007. V-Measure: A
       conditional entropy-based external cluster evaluation measure
       <https://aclweb.org/anthology/D/D07/D07-1043.pdf>`_
    """
    return balanced_homogeneity_completeness_v_measure(
        labels_true, labels_pred, beta=beta, reweigh=reweigh
    )[2]
