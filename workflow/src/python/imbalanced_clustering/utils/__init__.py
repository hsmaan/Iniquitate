import pyximport
import numpy

pyximport.install(setup_args={"include_dirs": numpy.get_include()}, reload_support=True)
from ._emi_cython import expected_mutual_information
from .contingency import pair_confusion_matrix, contingency_matrix
from .checks import check_clusterings
from .mi import mutual_info_score, entropy
from .avg import generalized_average
