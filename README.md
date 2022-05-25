# Probabilistic methods for nuclear data 2022

This presentation was given at the
*Joint ICTP-IAEA Advanced School/Workshop on Computational Nuclear Science and Engineering*
in Trieste, Italy. It follows the argument of Edwin Thompson Jaynes to motivate Bayesian
statistics, introduces the Bayesian update formula and specializes it to the case of
the multivariate normal distribution. The resulting equations are commonly known 
as the Generalized Least Squares (GLS) method. It is discussed how the required 
prior distribution can be constructed from a set of model predictions or by using a
covariance function. The latter possibility is at the heart of Gaussian process regression,
a non-parametric regression technique, another useful instrument in the toolbox of a
data scientist.

The GLS method is a general method and can be broadly applied beyond the domain
of nuclear data evaluation. To highlight this point, it is also shown how the
GLS method can be applied to the restoration of pictures with faces and the
tracking of a schematic mars rover. Even though state-of-the-art approaches would
certainly rely on other methods, the examples demonstrate the applicability
of statistical methods over a variety of application scenarios.

This repository also contains the scripts to reproduce the Bayesian inference
examples of the presentation including the plots and videos.
The resulting videos are stored in a separate folder `videos`. 

