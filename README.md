# Modal Interval Regression (MIR) in MATLAB

This repository provides a MATLAB implementation of a **Modal Interval Regression (MIR)** method.
The workflow is script-based: **run a single main script** and adjust parameters inside it.

## Contents
- **Main script**: runs the full MIR pipeline (the method entry point).
- **Three core functions**:
  1. **Conditional distribution estimation** (given \(X\), estimate \(\hat F_{Y|X}\).
  2. **Conditional MI estimation** (extract the conditional modal interval from the estimated conditional distribution).
  3. **ADMM prox for quantile loss** (proximal operator used in the ADMM solver for the quantile-loss formulation).


## Files
- `main_mir.m` — main script (run this)
- `cdfe.m` — conditional distribution estimator
- `cmie.m` — conditional MI estimator from \(\hat F_{Y|X}\) / \(\hat f_{Y|X}\)
- `per_prox.m` — prox operator for quantile loss used in ADMM

## Requirements
- MATLAB
- Require **Statistics and Machine Learning Toolbox**

## How to Use

You only need to modify the parameters in `main_mir.m` and run the script.

Typical parameters you may want to adjust include:

- **Coverage level** \(\alpha\) for the conditional MI (e.g., \(\alpha=0.5\) for the 50% MI)
- **Spline settings**, such as the spline degree/order, the number of knots, and any regularization/smoothing parameters
- **ADMM settings**, such as the penalty parameter, maximum number of iterations, and convergence tolerance
- **Smoothing parameters**


## Notes on Spline Knots (Important)

In this implementation, spline knots are selected **uniformly**.
⚠️ **To avoid a singular matrix**, when choosing knots, make sure that:
- **Between any two neighboring knots, there are at least two data points.**


## Output

Running `main_mir.m` typically returns and/or plots:

- Estimated **lower** and **upper** bound curves of the conditional MI
- 


