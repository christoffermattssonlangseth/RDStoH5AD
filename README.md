# RDStoH5AD

`RDStoH5AD` is a Python CLI with an R backend for converting `.rds` inputs to `.h5ad`.

This design is deliberate:

- R should own object deserialization for `Seurat`, `SingleCellExperiment`, and `SpatialExperiment`
- Python should own CLI ergonomics and process orchestration
- `.h5ad` writing should happen from a normalized `SingleCellExperiment` via `zellkonverter`

That is the shortest path to a converter that is both fast enough and correct.

## Why this design

For this problem, the expensive and fragile part is not CLI startup. It is:

- loading large R objects from disk
- understanding S4/Seurat assay structure
- preserving sparse matrices without accidental densification
- mapping metadata and embeddings into AnnData-compatible structures

Rewriting those parts in Rust would increase complexity quickly while giving limited benefit, because you still need R semantics to read and interpret many real-world `.rds` objects.

## Optimization priorities

The implementation is optimized around a few practical rules:

- keep sparse matrices sparse end-to-end
- avoid JSON/CSV intermediates for matrix payloads
- separate `inspect` from `convert`
- keep conversion focused on assay/layer/export work only
- let `zellkonverter` handle final `.h5ad` layout

## Requirements

Python:

- Python 3.10+

R packages:

- `jsonlite`
- `SingleCellExperiment`
- `SummarizedExperiment`
- `S4Vectors`
- `zellkonverter`

Optional, depending on input:

- `SeuratObject` for reading many Seurat objects
- `SpatialExperiment` for explicit spatial coordinate extraction from `SpatialExperiment`

Install the writer dependency in R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("zellkonverter", "SingleCellExperiment"))
```

## Install

From a local checkout:

```bash
python3 -m pip install -e .
```

From GitHub:

```bash
python3 -m pip install "git+https://github.com/christoffermattssonlangseth/RDStoH5AD.git"
```

## Usage

Inspect an input:

```bash
rds2h5ad inspect input.rds
```

Convert with defaults:

```bash
rds2h5ad convert input.rds output.h5ad
```

Validate a written H5AD:

```bash
rds2h5ad validate output.h5ad
```

Run the smoke tests:

```bash
python3 -m unittest discover -s tests -v
```

Choose an assay and AnnData `X` source explicitly:

```bash
rds2h5ad convert input.rds output.h5ad --assay SCT --x-layer data
```

Restrict exported layers and embeddings:

```bash
rds2h5ad convert input.rds output.h5ad \
  --assay RNA \
  --layers counts,data \
  --reduced-dims pca,umap
```

## Current support

- `Seurat`
- `SingleCellExperiment`
- `SpatialExperiment`
- plain lists with `obs` plus `expression` or `expr`, and optional `umap` / `coordinates`

## Notes

- `inspect` does not require `zellkonverter`.
- `convert` fails fast with a clear message if required R packages are missing.
- `validate` reads the generated `.h5ad` back through `zellkonverter` and prints a compact JSON summary.
- Seurat layers with incompatible dimensions, such as some `scale.data` cases, are skipped rather than densified or forced into the export.

## Repository Checklist

Before publishing, you should still make an explicit choice about:

- repository license
- GitHub repository URL in this README
