CP_MoCap algorithms

This is demo of CP decomposition-based algorithms for completion problem of motion capture data.
_Pattern Analysis and Applications. https://doi.org/10.1007/s10044-024-01342-4_

Algorithms:
Implementation of three CP-based tensor completion algorithms for reconstructing missing markers and frames in motion capture data:

CP: The CANDECOMP/PARAFAC decomposition
Smooth CP: CP decomposition with smoothness constraints on the vector factors
Sparse CP: CP decomposition with sparsity regularization on the vector factors


Requirements:
This code requires the following toolboxes:
1. [Tensor Toolbox](http://www.tensortoolbox.org/)
2. [Mocap Toolbox](https://github.com/mocaptoolbox).
