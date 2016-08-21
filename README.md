This toolbox reproduces the numerical results of the paper:

Justin Solomon, Fernando de Goes, Gabriel Peyré, Marco Cuturi, Adrian Butscher, Andy Nguyen, Tao Du, Leonidas Guibas, _Convolutional Wasserstein Distances: Efficient Optimal Transportation on Geometric Domains_, Proc. [SIGGRAPH 2015](http://s2015.siggraph.org/).

![Wasserstein barycenters of volumetric histograms](imgs/triangleinterp.jpg)


Content
-------

The main directories are:
* data/: images and meshes datasets.
* code/: code directory, with the following sub-directories:
    - cpp/: C++ implementation of the algorithm.
    - figures/: Matlab scripts to reproduce the figure of the article.
    - tests/: Matlab scripts to reproduce some further examples not shown in the article.
    - convolutional_wasserstein/: Matlab main functions implementing the algorithms.
    - toolbox/: Matlab helper functions.
    - blur_functions/ and mesh_functions/: Matlab function to compute heat kernels.
    - colors_functions/: exernal library (c) Pascal Getreuer.
    - image_blur/: external library imgaussian (c) Dirk-Jan Kroon

Copyright
-------

Copyright (c) 2015, Justin Solomon, Fernando de Goes, Gabriel Peyré, Marco Cuturi, Adrian Butscher, Andy Nguyen, Tao Du, Leonidas Guibas  
