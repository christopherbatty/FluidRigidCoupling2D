# FluidRigidCoupling2D
A simple implementation of fluid-rigid-body coupling in 2D, using the method of Batty et al. 2007

The original paper, "A Fast Variational Framework for Accurate Solid-Fluid Coupling", is here:
https://www.cs.ubc.ca/labs/imager/tr/2007/Batty_VariationalFluids/

This version uses face area weights (length in 2D) rather than volume fraction weights (area in 2D), as the former have been shown to be slightly more accurate by Ng et al, due to their connection with finite volume discretizations:
http://math.ewha.ac.kr/~chohong/publications/article_12_fluid_solid_coupling.pdf
