# Computing the square root of a low-rank perturbation of the scaled identity matrix



### About

The MATLAB code in this repository was used to produce the figures and tables in sections [sect. 1, 1] and [sect. 5, 1]. The results in the manuscript were produced by running these codes in MATLAB 2021b on a machine equipped with an AMD Ryzen 7 PRO 5850U and 32 GiB of RAM.



### Dependencies

The code leverages the capabilities of the Advanpix Multiprecision Computing Toolbox for MATLAB [3] to compute residual and forward error in higher-than-working precision.



### Running the tests

The `run_tests` script sets up the test environment and runs the test suite. Our code is organized in three folders:
* `data`, which contains the matrices provided as part of the Lingvo Framework [2] for TensorFlow in `.mat` format;
* `include`, which contains utility functions for computing the residual and forward error, for formatting numbers as LaTeX strings, and for timing the execution of MATLAB functions;
* `methods`, which contains our implementation of the algorithms for computing the square root of low-rank perturbations of the scaled identity; and
* `tests`, which contains the functions that run the algorithms, compute the execution time, residual, and forward error, and store the results to file.



### References

[1] M. Fasi, N. J. Higham, and X. Liu. [*Computing the square root of a low-rank perturbation of the scaled identity matrix*](http://eprints.maths.manchester.ac.uk/XXXX/). MIMS EPrint 2021.X, Manchester Institute for Mathematical Sciences, The University of Manchester, UK, December 2021.

[2] [Lingvo](https://github.com/tensorflow/lingvo/)
[3] [* Multiprecision Computing Toolbox for MATLAB](https://www.advanpix.com/). Advanpix LLC. Yokohama, Japan.

If you use the code in this repository, please cite the preprint [1].


### License

This software is distributed under the terms of the BSD 2-Clause "Simplified" license (see [LICENSE.md](./LICENSE.md)).
