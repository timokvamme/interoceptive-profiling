# Environment

These are the versions the deposited outputs were produced with. The whole
pipeline was run in a clean clone of this repository and every artifact was
compared against the deposited copy, so this list is what was actually used,
not what was intended.

## R

R 4.4.0.

| package | version |
|---|---|
| tidyverse | 2.0.0 |
| tidyLPA | 1.1.0 |
| mclust | 6.1.1 |
| lavaan | 0.6.21 |
| psych | 2.4.3 |
| igraph | 2.0.3 |

## Python, steps other than 15

Python 3.12.1.

| package | version |
|---|---|
| numpy | 2.4.6 |
| scipy | 1.17.0 |
| scikit-learn | 1.8.0 |
| pandas | 3.0.0 |
| matplotlib | 3.10.0 |

## Python, step 15 only

Step 15 needs `scikit-learn-extra`, which provides `KMedoids`. That package was
last released in 2022 and is compiled against NumPy 1. Importing it under NumPy
2 fails with `ImportError: numpy.core.multiarray failed to import`. So step 15
needs its own environment.

Python 3.10.18.

| package | version |
|---|---|
| numpy | 1.26.4 |
| scipy | 1.15.3 |
| scikit-learn | 1.7.2 |
| scikit-learn-extra | 0.3.0 |
| pandas | 2.3.3 |
| matplotlib | 3.10.9 |
| gower | latest |

Build it once:

```
python -m venv .venv15
.venv15/Scripts/python -m pip install "numpy<2" scikit-learn-extra gower pandas matplotlib scipy scikit-learn
```

Then point the runner at it for that step. `PYTHON_EXE` overrides the
interpreter used for Python steps, and `RSCRIPT` does the same for R:

```
set PYTHON_EXE=.venv15\Scripts\python.exe
python scripts/run_analysis.py 15
```

Running the whole pipeline without this override fails at step 15 and continues
past it. Nothing downstream reads step 15's output, so the rest of the results
are unaffected, but Supplementary Figure 16 will not be regenerated.

## Reproducibility

Two runs of the complete pipeline, in separate checkouts, produce identical
output with these exceptions:

- Three `sem_sensitivity_*` files and four `.xlsx` workbooks carry a creation
  timestamp. Their content is identical; only the recorded date differs.
- Step 15 agrees to about 1 part in 10^12 across NumPy builds. Test statistics,
  effect sizes, means and standard deviations are identical; only the far
  decimal places of very small p-values move. Within one environment it is
  exact.

Script 05 sets a seed, because the bootstrapped likelihood ratio test resamples.
Without it the `BLRT p` column changed on every run. The seed fixes the random
stream only; the models, log likelihood, AIC, BIC, entropy and class sizes were
identical across runs before it was added, and are unchanged by it.
