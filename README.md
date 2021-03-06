## [Nearly radial Neumann eigenfunctions on symmetric domains](http://arxiv.org/abs/1508.07019)
Collecion of scripts developed by Bartek Siudeja and Ben Young for the paper by Nilima Nigam, Bartek Siudeja and Ben Young (linked above). A separate set of Matlab scripts was developed by Nilima.

Most of the scripts are written in Python with [FEniCS](http://fenicsproject.org) ver. 1.4 as the major prerequisite (**these do not work with ver. 1.6+**). 

Some scripts also require scipy, sympy, mpmath and/or scikits.sparse.cholmod.

Finally, one script uses Sage interpreter in Python mode (no preparsing).

### Results

All sparse matrices generated by various scripts are included in appropriately named subfolders of the results folder.

Other text files contain computed eigenvalues/diagonal elements together with other performance metrics and error measurements.

The PDF file in the root folder shows all subdomains needed in our proof, together with their eigenvalues.

Finally, the mp4 movie shows the mesh, Dirichlet DOFs (yellow) and subdomains (blue) considered in the process of discovering the optimal sequence of 20 subdomains.
