import warnings
import os
import threading
import numbers
from collections.abc import Sequence
from contextlib import suppress

import numpy as np
import scipy.sparse as sp
from numpy.core.numeric import ComplexWarning


def _astype_copy_false(X):
    """Returns the copy=False parameter for
    {ndarray, csr_matrix, csc_matrix}.astype when possible,
    otherwise don't specify
    """
    if sp_version >= parse_version("1.1") or not sp.issparse(X):
        return {"copy": False}
    else:
        return {}


def _is_integral_float(y):
    return y.dtype.kind == "f" and np.all(y.astype(int) == y)


def check_consistent_length(*arrays):
    """Check that all arrays have consistent first dimensions.
    Checks whether all objects in arrays have the same shape or length.
    Parameters
    ----------
    *arrays : list or tuple of input objects.
        Objects that will be checked for consistent length.
    """

    lengths = [_num_samples(X) for X in arrays if X is not None]
    uniques = np.unique(lengths)
    if len(uniques) > 1:
        raise ValueError(
            "Found input variables with inconsistent numbers of samples: %r"
            % [int(l) for l in lengths]
        )


def is_multilabel(y):
    """Check if ``y`` is in a multilabel format.
    Parameters
    ----------
    y : ndarray of shape (n_samples,)
        Target values.
    Returns
    -------
    out : bool
        Return ``True``, if ``y`` is in a multilabel format, else ```False``.
    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.utils.multiclass import is_multilabel
    >>> is_multilabel([0, 1, 0, 1])
    False
    >>> is_multilabel([[1], [0, 2], []])
    False
    >>> is_multilabel(np.array([[1, 0], [0, 0]]))
    True
    >>> is_multilabel(np.array([[1], [0], [0]]))
    False
    >>> is_multilabel(np.array([[1, 0, 0]]))
    True
    """
    if hasattr(y, "__array__") or isinstance(y, Sequence):
        # DeprecationWarning will be replaced by ValueError, see NEP 34
        # https://numpy.org/neps/nep-0034-infer-dtype-is-object.html
        with warnings.catch_warnings():
            warnings.simplefilter("error", np.VisibleDeprecationWarning)
            try:
                y = np.asarray(y)
            except np.VisibleDeprecationWarning:
                # dtype=object should be provided explicitly for ragged arrays,
                # see NEP 34
                y = np.array(y, dtype=object)

    if not (hasattr(y, "shape") and y.ndim == 2 and y.shape[1] > 1):
        return False

    if sp.issparse(y):
        if isinstance(y, (sp.dok_matrix, sp.lil_matrix)):
            y = y.tocsr()
        return (
            len(y.data) == 0
            or np.unique(y.data).size == 1
            and (
                y.dtype.kind in "biu"
                or _is_integral_float(np.unique(y.data))  # bool, int, uint
            )
        )
    else:
        labels = np.unique(y)

        return len(labels) < 3 and (
            y.dtype.kind in "biu" or _is_integral_float(labels)  # bool, int, uint
        )


def type_of_target(y):
    """Determine the type of data indicated by the target.
    Note that this type is the most specific type that can be inferred.
    For example:
        * ``binary`` is more specific but compatible with ``multiclass``.
        * ``multiclass`` of integers is more specific but compatible with
          ``continuous``.
        * ``multilabel-indicator`` is more specific but compatible with
          ``multiclass-multioutput``.
    Parameters
    ----------
    y : array-like
    Returns
    -------
    target_type : str
        One of:
        * 'continuous': `y` is an array-like of floats that are not all
          integers, and is 1d or a column vector.
        * 'continuous-multioutput': `y` is a 2d array of floats that are
          not all integers, and both dimensions are of size > 1.
        * 'binary': `y` contains <= 2 discrete values and is 1d or a column
          vector.
        * 'multiclass': `y` contains more than two discrete values, is not a
          sequence of sequences, and is 1d or a column vector.
        * 'multiclass-multioutput': `y` is a 2d array that contains more
          than two discrete values, is not a sequence of sequences, and both
          dimensions are of size > 1.
        * 'multilabel-indicator': `y` is a label indicator matrix, an array
          of two dimensions with at least two columns, and at most 2 unique
          values.
        * 'unknown': `y` is array-like but none of the above, such as a 3d
          array, sequence of sequences, or an array of non-sequence objects.
    Examples
    --------
    >>> from sklearn.utils.multiclass import type_of_target
    >>> import numpy as np
    >>> type_of_target([0.1, 0.6])
    'continuous'
    >>> type_of_target([1, -1, -1, 1])
    'binary'
    >>> type_of_target(['a', 'b', 'a'])
    'binary'
    >>> type_of_target([1.0, 2.0])
    'binary'
    >>> type_of_target([1, 0, 2])
    'multiclass'
    >>> type_of_target([1.0, 0.0, 3.0])
    'multiclass'
    >>> type_of_target(['a', 'b', 'c'])
    'multiclass'
    >>> type_of_target(np.array([[1, 2], [3, 1]]))
    'multiclass-multioutput'
    >>> type_of_target([[1, 2]])
    'multilabel-indicator'
    >>> type_of_target(np.array([[1.5, 2.0], [3.0, 1.6]]))
    'continuous-multioutput'
    >>> type_of_target(np.array([[0, 1], [1, 1]]))
    'multilabel-indicator'
    """
    valid = (
        isinstance(y, Sequence) or sp.issparse(y) or hasattr(y, "__array__")
    ) and not isinstance(y, str)

    if not valid:
        raise ValueError(
            "Expected array-like (array or non-string sequence), got %r" % y
        )

    sparse_pandas = y.__class__.__name__ in ["SparseSeries", "SparseArray"]
    if sparse_pandas:
        raise ValueError("y cannot be class 'SparseSeries' or 'SparseArray'")

    if is_multilabel(y):
        return "multilabel-indicator"

    # DeprecationWarning will be replaced by ValueError, see NEP 34
    # https://numpy.org/neps/nep-0034-infer-dtype-is-object.html
    with warnings.catch_warnings():
        warnings.simplefilter("error", np.VisibleDeprecationWarning)
        try:
            y = np.asarray(y)
        except np.VisibleDeprecationWarning:
            # dtype=object should be provided explicitly for ragged arrays,
            # see NEP 34
            y = np.asarray(y, dtype=object)

    # The old sequence of sequences format
    try:
        if (
            not hasattr(y[0], "__array__")
            and isinstance(y[0], Sequence)
            and not isinstance(y[0], str)
        ):
            raise ValueError(
                "You appear to be using a legacy multi-label data"
                " representation. Sequence of sequences are no"
                " longer supported; use a binary array or sparse"
                " matrix instead - the MultiLabelBinarizer"
                " transformer can convert to this format."
            )
    except IndexError:
        pass

    # Invalid inputs
    if y.ndim > 2 or (y.dtype == object and len(y) and not isinstance(y.flat[0], str)):
        return "unknown"  # [[[1, 2]]] or [obj_1] and not ["label_1"]

    if y.ndim == 2 and y.shape[1] == 0:
        return "unknown"  # [[]]

    if y.ndim == 2 and y.shape[1] > 1:
        suffix = "-multioutput"  # [[1, 2], [1, 2]]
    else:
        suffix = ""  # [1, 2, 3] or [[1], [2], [3]]

    # check float and contains non-integer float values
    if y.dtype.kind == "f" and np.any(y != y.astype(int)):
        # [.1, .2, 3] or [[.1, .2, 3]] or [[1., .2]] and not [1., 2., 3.]
        _assert_all_finite(y)
        return "continuous" + suffix

    if (len(np.unique(y)) > 2) or (y.ndim >= 2 and len(y[0]) > 1):
        return "multiclass" + suffix  # [1, 2, 3] or [[1., 2., 3]] or [[1, 2]]
    else:
        return "binary"  # [1, 2] or [["a"], ["b"]]


def _object_dtype_isnan(X):
    return X != X


_global_config = {
    "assume_finite": bool(os.environ.get("SKLEARN_ASSUME_FINITE", False)),
    "working_memory": int(os.environ.get("SKLEARN_WORKING_MEMORY", 1024)),
    "print_changed_only": True,
    "display": "text",
}
_threadlocal = threading.local()


def _get_threadlocal_config():
    """Get a threadlocal **mutable** configuration. If the configuration
    does not exist, copy the default global configuration."""
    if not hasattr(_threadlocal, "global_config"):
        _threadlocal.global_config = _global_config.copy()
    return _threadlocal.global_config


def _get_config():
    """Retrieve current values for configuration set by :func:`set_config`.
    Returns
    -------
    config : dict
        Keys are parameter names that can be passed to :func:`set_config`.
    See Also
    --------
    config_context : Context manager for global scikit-learn configuration.
    set_config : Set global scikit-learn configuration.
    """
    # Return a copy of the threadlocal configuration so that users will
    # not be able to modify the configuration with the returned dict.
    return _get_threadlocal_config().copy()


def _safe_accumulator_op(op, x, *args, **kwargs):
    """
    This function provides numpy accumulator functions with a float64 dtype
    when used on a floating point input. This prevents accumulator overflow on
    smaller floating point dtypes.
    Parameters
    ----------
    op : function
        A numpy accumulator function such as np.mean or np.sum.
    x : ndarray
        A numpy array to apply the accumulator function.
    *args : positional arguments
        Positional arguments passed to the accumulator function after the
        input x.
    **kwargs : keyword arguments
        Keyword arguments passed to the accumulator function.
    Returns
    -------
    result
        The output of the accumulator function passed to this function.
    """
    if np.issubdtype(x.dtype, np.floating) and x.dtype.itemsize < 8:
        result = op(x, *args, **kwargs, dtype=np.float64)
    else:
        result = op(x, *args, **kwargs)
    return


def _assert_all_finite(X, allow_nan=False, msg_dtype=None):
    """Like assert_all_finite, but only for ndarray."""
    if _get_config()["assume_finite"]:
        return
    X = np.asanyarray(X)
    # First try an O(n) time, O(1) space solution for the common case that
    # everything is finite; fall back to O(n) space np.isfinite to prevent
    # false positives from overflow in sum method. The sum is also calculated
    # safely to reduce dtype induced overflows.
    is_float = X.dtype.kind in "fc"
    if is_float and (np.isfinite(_safe_accumulator_op(np.sum, X))):
        pass
    elif is_float:
        msg_err = "Input contains {} or a value too large for {!r}."
        if (
            allow_nan
            and np.isinf(X).any()
            or not allow_nan
            and not np.isfinite(X).all()
        ):
            type_err = "infinity" if allow_nan else "NaN, infinity"
            raise ValueError(
                msg_err.format(
                    type_err, msg_dtype if msg_dtype is not None else X.dtype
                )
            )
    # for object dtype data, we only check for NaNs (GH-13254)
    elif X.dtype == np.dtype("object") and not allow_nan:
        if _object_dtype_isnan(X).any():
            raise ValueError("Input contains NaN")


def _check_large_sparse(X, accept_large_sparse=False):
    """Raise a ValueError if X has 64bit indices and accept_large_sparse=False"""
    if not accept_large_sparse:
        supported_indices = ["int32"]
        if X.getformat() == "coo":
            index_keys = ["col", "row"]
        elif X.getformat() in ["csr", "csc", "bsr"]:
            index_keys = ["indices", "indptr"]
        else:
            return
        for key in index_keys:
            indices_datatype = getattr(X, key).dtype
            if indices_datatype not in supported_indices:
                raise ValueError(
                    "Only sparse matrices with 32-bit integer"
                    " indices are accepted. Got %s indices." % indices_datatype
                )


def _ensure_sparse_format(
    spmatrix, accept_sparse, dtype, copy, force_all_finite, accept_large_sparse
):
    """Convert a sparse matrix to a given format.
    Checks the sparse format of spmatrix and converts if necessary.
    Parameters
    ----------
    spmatrix : sparse matrix
        Input to validate and convert.
    accept_sparse : str, bool or list/tuple of str
        String[s] representing allowed sparse matrix formats ('csc',
        'csr', 'coo', 'dok', 'bsr', 'lil', 'dia'). If the input is sparse but
        not in the allowed format, it will be converted to the first listed
        format. True allows the input to be any format. False means
        that a sparse matrix input will raise an error.
    dtype : str, type or None
        Data type of result. If None, the dtype of the input is preserved.
    copy : bool
        Whether a forced copy will be triggered. If copy=False, a copy might
        be triggered by a conversion.
    force_all_finite : bool or 'allow-nan'
        Whether to raise an error on np.inf, np.nan, pd.NA in X. The
        possibilities are:
        - True: Force all values of X to be finite.
        - False: accepts np.inf, np.nan, pd.NA in X.
        - 'allow-nan': accepts only np.nan and pd.NA values in X. Values cannot
          be infinite.
        .. versionadded:: 0.20
           ``force_all_finite`` accepts the string ``'allow-nan'``.
        .. versionchanged:: 0.23
           Accepts `pd.NA` and converts it into `np.nan`
    Returns
    -------
    spmatrix_converted : sparse matrix.
        Matrix that is ensured to have an allowed type.
    """
    if dtype is None:
        dtype = spmatrix.dtype

    changed_format = False

    if isinstance(accept_sparse, str):
        accept_sparse = [accept_sparse]

    # Indices dtype validation
    _check_large_sparse(spmatrix, accept_large_sparse)

    if accept_sparse is False:
        raise TypeError(
            "A sparse matrix was passed, but dense "
            "data is required. Use X.toarray() to "
            "convert to a dense numpy array."
        )
    elif isinstance(accept_sparse, (list, tuple)):
        if len(accept_sparse) == 0:
            raise ValueError(
                "When providing 'accept_sparse' "
                "as a tuple or list, it must contain at "
                "least one string value."
            )
        # ensure correct sparse format
        if spmatrix.format not in accept_sparse:
            # create new with correct sparse
            spmatrix = spmatrix.asformat(accept_sparse[0])
            changed_format = True
    elif accept_sparse is not True:
        # any other type
        raise ValueError(
            "Parameter 'accept_sparse' should be a string, "
            "boolean or list of strings. You provided "
            "'accept_sparse={}'.".format(accept_sparse)
        )

    if dtype != spmatrix.dtype:
        # convert dtype
        spmatrix = spmatrix.astype(dtype)
    elif copy and not changed_format:
        # force copy
        spmatrix = spmatrix.copy()

    if force_all_finite:
        if not hasattr(spmatrix, "data"):
            warnings.warn(
                "Can't check %s sparse matrix for nan or inf." % spmatrix.format,
                stacklevel=2,
            )
        else:
            _assert_all_finite(spmatrix.data, allow_nan=force_all_finite == "allow-nan")

    return spmatrix


def _ensure_no_complex_data(array):
    if (
        hasattr(array, "dtype")
        and array.dtype is not None
        and hasattr(array.dtype, "kind")
        and array.dtype.kind == "c"
    ):
        raise ValueError("Complex data not supported\n{}\n".format(array))


def _num_samples(x):
    """Return number of samples in array-like x."""
    message = "Expected sequence or array-like, got %s" % type(x)
    if hasattr(x, "fit") and callable(x.fit):
        # Don't get num_samples from an ensembles length!
        raise TypeError(message)

    if not hasattr(x, "__len__") and not hasattr(x, "shape"):
        if hasattr(x, "__array__"):
            x = np.asarray(x)
        else:
            raise TypeError(message)

    if hasattr(x, "shape") and x.shape is not None:
        if len(x.shape) == 0:
            raise TypeError(
                "Singleton array %r cannot be considered a valid collection." % x
            )
        # Check that shape is returning an integer or default to len
        # Dask dataframes may not return numeric shape[0] value
        if isinstance(x.shape[0], numbers.Integral):
            return x.shape[0]

    try:
        return len(x)
    except TypeError as type_error:
        raise TypeError(message) from type_error


def check_array(
    array,
    accept_sparse=False,
    *,
    accept_large_sparse=True,
    dtype="numeric",
    order=None,
    copy=False,
    force_all_finite=True,
    ensure_2d=True,
    allow_nd=False,
    ensure_min_samples=1,
    ensure_min_features=1,
    estimator=None,
):

    """Input validation on an array, list, sparse matrix or similar.
    By default, the input is checked to be a non-empty 2D array containing
    only finite values. If the dtype of the array is object, attempt
    converting to float, raising on failure.
    Parameters
    ----------
    array : object
        Input object to check / convert.
    accept_sparse : str, bool or list/tuple of str, default=False
        String[s] representing allowed sparse matrix formats, such as 'csc',
        'csr', etc. If the input is sparse but not in the allowed format,
        it will be converted to the first listed format. True allows the input
        to be any format. False means that a sparse matrix input will
        raise an error.
    accept_large_sparse : bool, default=True
        If a CSR, CSC, COO or BSR sparse matrix is supplied and accepted by
        accept_sparse, accept_large_sparse=False will cause it to be accepted
        only if its indices are stored with a 32-bit dtype.
        .. versionadded:: 0.20
    dtype : 'numeric', type, list of type or None, default='numeric'
        Data type of result. If None, the dtype of the input is preserved.
        If "numeric", dtype is preserved unless array.dtype is object.
        If dtype is a list of types, conversion on the first type is only
        performed if the dtype of the input is not in the list.
    order : {'F', 'C'} or None, default=None
        Whether an array will be forced to be fortran or c-style.
        When order is None (default), then if copy=False, nothing is ensured
        about the memory layout of the output array; otherwise (copy=True)
        the memory layout of the returned array is kept as close as possible
        to the original array.
    copy : bool, default=False
        Whether a forced copy will be triggered. If copy=False, a copy might
        be triggered by a conversion.
    force_all_finite : bool or 'allow-nan', default=True
        Whether to raise an error on np.inf, np.nan, pd.NA in array. The
        possibilities are:
        - True: Force all values of array to be finite.
        - False: accepts np.inf, np.nan, pd.NA in array.
        - 'allow-nan': accepts only np.nan and pd.NA values in array. Values
          cannot be infinite.
        .. versionadded:: 0.20
           ``force_all_finite`` accepts the string ``'allow-nan'``.
        .. versionchanged:: 0.23
           Accepts `pd.NA` and converts it into `np.nan`
    ensure_2d : bool, default=True
        Whether to raise a value error if array is not 2D.
    allow_nd : bool, default=False
        Whether to allow array.ndim > 2.
    ensure_min_samples : int, default=1
        Make sure that the array has a minimum number of samples in its first
        axis (rows for a 2D array). Setting to 0 disables this check.
    ensure_min_features : int, default=1
        Make sure that the 2D array has some minimum number of features
        (columns). The default value of 1 rejects empty datasets.
        This check is only enforced when the input data has effectively 2
        dimensions or is originally 1D and ``ensure_2d`` is True. Setting to 0
        disables this check.
    estimator : str or estimator instance, default=None
        If passed, include the name of the estimator in warning messages.
    Returns
    -------
    array_converted : object
        The converted and validated array.
    """
    if isinstance(array, np.matrix):
        warnings.warn(
            "np.matrix usage is deprecated in 1.0 and will raise a TypeError "
            "in 1.2. Please convert to a numpy array with np.asarray. For "
            "more information see: "
            "https://numpy.org/doc/stable/reference/generated/numpy.matrix.html",  # noqa
            FutureWarning,
        )

    # store reference to original array to check if copy is needed when
    # function returns
    array_orig = array

    # store whether originally we wanted numeric dtype
    dtype_numeric = isinstance(dtype, str) and dtype == "numeric"

    dtype_orig = getattr(array, "dtype", None)
    if not hasattr(dtype_orig, "kind"):
        # not a data type (e.g. a column named dtype in a pandas DataFrame)
        dtype_orig = None

    # check if the object contains several dtypes (typically a pandas
    # DataFrame), and store them. If not, store None.
    dtypes_orig = None
    has_pd_integer_array = False
    if hasattr(array, "dtypes") and hasattr(array.dtypes, "__array__"):
        # throw warning if columns are sparse. If all columns are sparse, then
        # array.sparse exists and sparsity will be preserved (later).
        with suppress(ImportError):
            from pandas.api.types import is_sparse

            if not hasattr(array, "sparse") and array.dtypes.apply(is_sparse).any():
                warnings.warn(
                    "pandas.DataFrame with sparse columns found."
                    "It will be converted to a dense numpy array."
                )

        dtypes_orig = list(array.dtypes)
        # pandas boolean dtype __array__ interface coerces bools to objects
        for i, dtype_iter in enumerate(dtypes_orig):
            if dtype_iter.kind == "b":
                dtypes_orig[i] = np.dtype(object)
            elif dtype_iter.name.startswith(("Int", "UInt")):
                # name looks like an Integer Extension Array, now check for
                # the dtype
                with suppress(ImportError):
                    from pandas import (
                        Int8Dtype,
                        Int16Dtype,
                        Int32Dtype,
                        Int64Dtype,
                        UInt8Dtype,
                        UInt16Dtype,
                        UInt32Dtype,
                        UInt64Dtype,
                    )

                    if isinstance(
                        dtype_iter,
                        (
                            Int8Dtype,
                            Int16Dtype,
                            Int32Dtype,
                            Int64Dtype,
                            UInt8Dtype,
                            UInt16Dtype,
                            UInt32Dtype,
                            UInt64Dtype,
                        ),
                    ):
                        has_pd_integer_array = True

        if all(isinstance(dtype, np.dtype) for dtype in dtypes_orig):
            dtype_orig = np.result_type(*dtypes_orig)

    if dtype_numeric:
        if dtype_orig is not None and dtype_orig.kind == "O":
            # if input is object, convert to float.
            dtype = np.float64
        else:
            dtype = None

    if isinstance(dtype, (list, tuple)):
        if dtype_orig is not None and dtype_orig in dtype:
            # no dtype conversion required
            dtype = None
        else:
            # dtype conversion required. Let's select the first element of the
            # list of accepted types.
            dtype = dtype[0]

    if has_pd_integer_array:
        # If there are any pandas integer extension arrays,
        array = array.astype(dtype)

    if force_all_finite not in (True, False, "allow-nan"):
        raise ValueError(
            'force_all_finite should be a bool or "allow-nan". Got {!r} instead'.format(
                force_all_finite
            )
        )

    if estimator is not None:
        if isinstance(estimator, str):
            estimator_name = estimator
        else:
            estimator_name = estimator.__class__.__name__
    else:
        estimator_name = "Estimator"
    context = " by %s" % estimator_name if estimator is not None else ""

    # When all dataframe columns are sparse, convert to a sparse array
    if hasattr(array, "sparse") and array.ndim > 1:
        # DataFrame.sparse only supports `to_coo`
        array = array.sparse.to_coo()
        if array.dtype == np.dtype("object"):
            unique_dtypes = set([dt.subtype.name for dt in array_orig.dtypes])
            if len(unique_dtypes) > 1:
                raise ValueError(
                    "Pandas DataFrame with mixed sparse extension arrays "
                    "generated a sparse matrix with object dtype which "
                    "can not be converted to a scipy sparse matrix."
                    "Sparse extension arrays should all have the same "
                    "numeric type."
                )

    if sp.issparse(array):
        _ensure_no_complex_data(array)
        array = _ensure_sparse_format(
            array,
            accept_sparse=accept_sparse,
            dtype=dtype,
            copy=copy,
            force_all_finite=force_all_finite,
            accept_large_sparse=accept_large_sparse,
        )
    else:
        # If np.array(..) gives ComplexWarning, then we convert the warning
        # to an error. This is needed because specifying a non complex
        # dtype to the function converts complex to real dtype,
        # thereby passing the test made in the lines following the scope
        # of warnings context manager.
        with warnings.catch_warnings():
            try:
                warnings.simplefilter("error", ComplexWarning)
                if dtype is not None and np.dtype(dtype).kind in "iu":
                    # Conversion float -> int should not contain NaN or
                    # inf (numpy#14412). We cannot use casting='safe' because
                    # then conversion float -> int would be disallowed.
                    array = np.asarray(array, order=order)
                    if array.dtype.kind == "f":
                        _assert_all_finite(array, allow_nan=False, msg_dtype=dtype)
                    array = array.astype(dtype, casting="unsafe", copy=False)
                else:
                    array = np.asarray(array, order=order, dtype=dtype)
            except ComplexWarning as complex_warning:
                raise ValueError(
                    "Complex data not supported\n{}\n".format(array)
                ) from complex_warning

        # It is possible that the np.array(..) gave no warning. This happens
        # when no dtype conversion happened, for example dtype = None. The
        # result is that np.array(..) produces an array of complex dtype
        # and we need to catch and raise exception for such cases.
        _ensure_no_complex_data(array)

        if ensure_2d:
            # If input is scalar raise error
            if array.ndim == 0:
                raise ValueError(
                    "Expected 2D array, got scalar array instead:\narray={}.\n"
                    "Reshape your data either using array.reshape(-1, 1) if "
                    "your data has a single feature or array.reshape(1, -1) "
                    "if it contains a single sample.".format(array)
                )
            # If input is 1D raise error
            if array.ndim == 1:
                raise ValueError(
                    "Expected 2D array, got 1D array instead:\narray={}.\n"
                    "Reshape your data either using array.reshape(-1, 1) if "
                    "your data has a single feature or array.reshape(1, -1) "
                    "if it contains a single sample.".format(array)
                )

        # make sure we actually converted to numeric:
        if dtype_numeric and array.dtype.kind in "OUSV":
            warnings.warn(
                "Arrays of bytes/strings is being converted to decimal "
                "numbers if dtype='numeric'. This behavior is deprecated in "
                "0.24 and will be removed in 1.1 (renaming of 0.26). Please "
                "convert your data to numeric values explicitly instead.",
                FutureWarning,
                stacklevel=2,
            )
            try:
                array = array.astype(np.float64)
            except ValueError as e:
                raise ValueError(
                    "Unable to convert array of bytes/strings "
                    "into decimal numbers with dtype='numeric'"
                ) from e
        if not allow_nd and array.ndim >= 3:
            raise ValueError(
                "Found array with dim %d. %s expected <= 2."
                % (array.ndim, estimator_name)
            )

        if force_all_finite:
            _assert_all_finite(array, allow_nan=force_all_finite == "allow-nan")

    if ensure_min_samples > 0:
        n_samples = _num_samples(array)
        if n_samples < ensure_min_samples:
            raise ValueError(
                "Found array with %d sample(s) (shape=%s) while a"
                " minimum of %d is required%s."
                % (n_samples, array.shape, ensure_min_samples, context)
            )

    if ensure_min_features > 0 and array.ndim == 2:
        n_features = array.shape[1]
        if n_features < ensure_min_features:
            raise ValueError(
                "Found array with %d feature(s) (shape=%s) while"
                " a minimum of %d is required%s."
                % (n_features, array.shape, ensure_min_features, context)
            )

    if copy and np.may_share_memory(array, array_orig):
        array = np.array(array, dtype=dtype, order=order)

    return array


def check_clusterings(labels_true, labels_pred):
    """Check that the labels arrays are 1D and of same dimension.
    Parameters
    ----------
    labels_true : array-like of shape (n_samples,)
        The true labels.
    labels_pred : array-like of shape (n_samples,)
        The predicted labels.
    """
    labels_true = check_array(
        labels_true,
        ensure_2d=False,
        ensure_min_samples=0,
        dtype=None,
    )

    labels_pred = check_array(
        labels_pred,
        ensure_2d=False,
        ensure_min_samples=0,
        dtype=None,
    )

    type_label = type_of_target(labels_true)
    type_pred = type_of_target(labels_pred)

    if "continuous" in (type_pred, type_label):
        msg = (
            "Clustering metrics expects discrete values but received"
            f" {type_label} values for label, and {type_pred} values "
            "for target"
        )
        warnings.warn(msg, UserWarning)

    # input checks
    if labels_true.ndim != 1:
        raise ValueError("labels_true must be 1D: shape is %r" % (labels_true.shape,))
    if labels_pred.ndim != 1:
        raise ValueError("labels_pred must be 1D: shape is %r" % (labels_pred.shape,))
    check_consistent_length(labels_true, labels_pred)

    return labels_true, labels_pred
