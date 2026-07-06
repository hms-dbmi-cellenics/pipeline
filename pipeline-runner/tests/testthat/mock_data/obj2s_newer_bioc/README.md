# obj2s fixtures serialized with a newer Bioconductor

`sce_seqinfo.rds` is a `SingleCellExperiment` (built from `SeuratObject::pbmc_small`)
whose `rowRanges` carry a `Seqinfo` slot. It was saved on a machine running
Bioconductor >= 3.21, where the `Seqinfo` S4 class lives in its own standalone
`Seqinfo` package (it used to live in `GenomeInfoDb`). The serialized class
therefore records `package = "Seqinfo"`.

The pipeline image runs Bioconductor 3.19 (GenomeInfoDb 1.40, no standalone
`Seqinfo` package), so reading this object and touching `rowRanges` (via
`updateObject()`) used to fail with `unable to find required package 'Seqinfo'`.
`reconstruct_sce()` now drops `rowRanges` before any accessor runs, so the object
loads regardless.

Regenerate on a machine with a recent Bioconductor:

```r
library(SingleCellExperiment); library(GenomicRanges); library(Seurat)
data("pbmc_small", package = "SeuratObject")
pbmc_small$seurat_clusters <- pbmc_small$RNA_snn_res.1
sce <- Seurat::as.SingleCellExperiment(pbmc_small)
gr <- GRanges(rep("chr1", nrow(sce)), IRanges(seq_len(nrow(sce)), width = 1))
names(gr) <- rownames(sce)
rowRanges(sce) <- split(gr, factor(rownames(sce), levels = rownames(sce)))
saveRDS(sce, "sce_seqinfo.rds")
```
