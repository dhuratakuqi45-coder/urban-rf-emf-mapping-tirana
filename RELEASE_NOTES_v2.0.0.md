# Version 2.0.0

This release adds the reproducible secondary spatial-analysis workflow used for the revised WSEAS manuscript.

## Added

- Global Moran's I permutation tests;
- sensitivity analysis for k = 4, 6, 8 and 10;
- Local Moran cluster classification with Benjamini-Hochberg FDR control;
- permutation-based distance-decay statistics;
- nearest-neighbour versus IDW leave-one-out validation;
- exact paired sign tests and relative RMSE reduction;
- Voronoi and Local Moran figures for KT and MKA;
- documentation and a dedicated dependency file.

## Reproducibility

All permutation procedures use fixed seeds and 999 permutations. The workflow reads the same raw campaign workbooks archived for the original study; no new field measurements are introduced.

