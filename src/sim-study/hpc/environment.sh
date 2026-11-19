# Shared by preparation, array tasks and collection. Edit the module version
# after `module avail math/R` on bwUniCluster; use the SAME one throughout.
module load "${R_MODULE:-math/R/4.5.1}"

# One CPU per run: stop BLAS, OpenMP and Stan from silently adding threads.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 STAN_NUM_THREADS=1 TBB_NUM_THREADS=1
export MAKEFLAGS=-j1
# Rscript below reads the project's .Rprofile and activates renv. Packages
# must already be installed on Linux; array tasks never install packages.
