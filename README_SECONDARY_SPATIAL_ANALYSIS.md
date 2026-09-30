# Secondary spatial analysis workflow

This package reproduces the additional spatial analyses used in the WSEAS manuscript based on the 171-session Tirana RF-EMF campaign previously archived with the project.

## Analyses included

- reconstruction of six band-integrated electric-field values from 531 spectral bins using root-sum-square aggregation;
- session-wise multi-frequency exposure quotient calculation;
- Voronoi partitions clipped to the sampled convex hull;
- Global Moran's I with 999 permutations;
- neighbourhood sensitivity for k = 4, 6, 8 and 10 nearest neighbours;
- Local Moran's I with Benjamini-Hochberg false-discovery-rate correction;
- permutation-based distance-decay analysis using Spearman rank correlation;
- leave-one-out comparison of nearest-neighbour and inverse-distance-weighting prediction;
- exact two-sided sign tests for paired absolute prediction errors.

## Input files

Place the script in a directory containing the original campaign workbooks, or in a parent directory from which they can be found recursively.

Expected coordinate workbooks:

- `Flete_Matje*.xlsx`

Expected spectrum workbooks:

- `Fusha_E_*.xlsx`

The script matches measurement sessions to coordinates by site, campaign date and within-day acquisition order, following the structure of the archived campaign files.

## Installation

Python 3.11 or later is recommended.

```bash
python -m pip install -r requirements-secondary-spatial.txt
```

## Execution

```bash
python secondary_spatial_analysis.py
```

The script uses fixed random seeds and 999 permutations for reproducibility.

## Outputs

Results are written to an `analysis` subdirectory:

- `session_band_values.csv`: reconstructed session-level band values;
- `summary.json`: descriptive statistics, spatial tests and cross-validation results;
- `voronoi_KT.png` and `voronoi_MKA.png`: Voronoi measurement partitions;
- `local_moran_KT.png` and `local_moran_MKA.png`: FDR-controlled Local Moran cluster maps.

## Interpretation

Voronoi cells are descriptive nearest-site partitions and should not be interpreted as a continuous prediction surface. The distance-decay statistic is a permutation-based association diagnostic, not a propagation model. Global and Local Moran results depend on the spatial-weights definition; the included k-sensitivity analysis should therefore be reported alongside the primary k = 6 result.

## Data citation

The measurement dataset underlying this workflow is archived at:

https://doi.org/10.5281/zenodo.21625221

The original campaign and its primary analysis are described in:

Kuqi D, Myzeqari B, Zeneli E, Kuqali M, Askushaj D, Plaku M. Spatial and Temporal Characterization of RF-EMF Exposure from Mobile Base Stations in Urban Tirana: A Multi-Band Field Measurement Study. Applied Sciences. 2026;16:8676. https://doi.org/10.3390/app16178676.

