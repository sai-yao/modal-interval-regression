# Modal Interval Regression (MIR) in MATLAB

This repository provides a MATLAB implementation of a **Modal Interval Regression (MIR)** method.
The workflow is script-based: **run a single main script** and adjust parameters inside it.

## Contents
- **Main script**: runs the full MIR pipeline (the method entry point).
- **Three core functions**:
  1. **Conditional distribution estimation**.
  2. **Conditional MI estimation** (extract the conditional modal interval from the estimated conditional distribution).
  3. **ADMM prox for quantile loss** (proximal operator used in the ADMM solver for the quantile-loss formulation).


## Files
- `main_mir.m` — main script (run this)
- `cdfe.m` — conditional distribution estimator
- `cmie.m` — conditional MI estimator
- `per_prox.m` — prox operator for quantile loss used in ADMM

## Requirements
- MATLAB
- Require **Statistics and Machine Learning Toolbox**

## How to Use

You only need to modify the parameters in `main_mir.m` and run the script.

Before running, provide the input data vectors:
- `X` — traininiig set covariate vector
- `d_data` — training set response (observed) vector
- `test_X` — test set covariate vector
- `test_Y` — test set response (observed) vector

Typical parameters you may want to adjust include:

- **Coverage level (`alpha`)**: the MI coverage level (e.g., `alpha = 0.5` for the 50% MI).
- **Number of polynomial pieces (`b`)**: the number of polynomial segments over the domain.
- **Domain endpoints (`start_value`, `end_value`)**: the left and right endpoints of the covariate domain.
- **Polynomial degree (`d`)**: the degree of each polynomial piece.
- **Spline smoothness (`rho`)**: the smoothness/continuity level imposed on the spline.
- **Smoothing parameter (`lambda`)**: the regularization strength (larger values typically yield smoother curves).
- **ADMM iterations (`iter`)**: the number of ADMM iterations (maximum iteration count) used by the solver.
- **CWC penalty parameter (`eta`)**: the penalty strength parameter in mCWC.

## Notes on Spline Knots (Important)

In this implementation, spline knots are selected **uniformly**.

⚠️ **To avoid a singular matrix**, when choosing knots, make sure that:
- **Between any two neighboring knots, there are at least two data points.**


## Output

Running `main_mir.m` produces a figure showing the estimated conditional MI bounds (lower and upper curves), saves the fitted coefficient vector `c_star` to `mir_coefficient.mat`, and reports interval quality metrics **mCWC** value.

## Related publication

This code is associated with the following paper:

Yao, S., Araki, Y., and Iwata, O. (2026). Nonlinear modal interval regression for bivariate data analysis. *Journal of Applied Statistics*, 1–28. https://doi.org/10.1080/02664763.2026.2667949


