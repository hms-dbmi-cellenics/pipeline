# Unisample integration works

    Code
      str(integrated_scdata)
    Output
      Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        ..@ assays      :List of 1
        .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
        .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. ..@ i       : int [1:4456] 1 5 8 11 22 30 33 34 36 38 ...
        .. .. .. .. .. ..@ p       : int [1:81] 0 47 99 149 205 258 306 342 387 423 ...
        .. .. .. .. .. ..@ Dim     : int [1:2] 230 80
        .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. ..$ : chr [1:80] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. ..@ x       : num [1:4456] 1 1 3 1 1 4 1 5 1 1 ...
        .. .. .. .. .. ..@ factors : list()
        .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. ..@ i       : int [1:4456] 1 5 8 11 22 30 33 34 36 38 ...
        .. .. .. .. .. ..@ p       : int [1:81] 0 47 99 149 205 258 306 342 387 423 ...
        .. .. .. .. .. ..@ Dim     : int [1:2] 230 80
        .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. ..$ : chr [1:80] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. ..@ x       : num [1:4456] 4.97 4.97 6.06 4.97 4.97 ...
        .. .. .. .. .. ..@ factors : list()
        .. .. .. ..@ scale.data   : num [1:230, 1:80] -0.409 1.64 -0.428 -1.375 -0.329 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. ..$ : chr [1:80] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. ..@ assay.orig   : NULL
        .. .. .. ..@ var.features : chr [1:230] "PPBP" "IGLL5" "VDAC3" "CD1C" ...
        .. .. .. ..@ meta.features:'data.frame':	230 obs. of  5 variables:
        .. .. .. .. ..$ vst.mean                 : num [1:230] 0.388 0.6 0.7 13.425 0.3 ...
        .. .. .. .. ..$ vst.variance             : num [1:230] 1.025 1.281 4.365 725.463 0.871 ...
        .. .. .. .. ..$ vst.variance.expected    : num [1:230] 1.141 2.708 4.018 745.526 0.642 ...
        .. .. .. .. ..$ vst.variance.standardized: num [1:230] 0.898 0.473 1.086 0.973 1.356 ...
        .. .. .. .. ..$ vst.variable             : logi [1:230] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. ..@ misc         : Named list()
        .. .. .. ..@ key          : Named chr "rna_"
        .. .. .. .. ..- attr(*, "names")= chr ""
        ..@ meta.data   :'data.frame':	80 obs. of  5 variables:
        .. ..$ orig.ident  : chr [1:80] "SeuratProject" "SeuratProject" "SeuratProject" "SeuratProject" ...
        .. ..$ nCount_RNA  : num [1:80] 70 85 87 127 173 70 64 72 52 100 ...
        .. ..$ nFeature_RNA: int [1:80] 47 52 50 56 53 48 36 45 36 41 ...
        .. ..$ samples     : chr [1:80] "123abc" "123abc" "123abc" "123abc" ...
        .. ..$ cells_id    : int [1:80] 0 1 2 3 4 5 6 7 8 9 ...
        ..@ active.assay: chr "RNA"
        ..@ active.ident: Factor w/ 1 level "SeuratProject": 1 1 1 1 1 1 1 1 1 1 ...
        .. ..- attr(*, "names")= chr [1:80] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        ..@ graphs      : list()
        ..@ neighbors   : list()
        ..@ reductions  :List of 1
        .. ..$ pca:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
        .. .. .. ..@ cell.embeddings           : num [1:80, 1:39] 3.12 3.56 2.4 3.43 2.78 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:80] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. ..$ : chr [1:39] "PC_1" "PC_2" "PC_3" "PC_4" ...
        .. .. .. ..@ feature.loadings          : num [1:230, 1:39] 0.05711 0.00738 0.03005 -0.04766 0.05598 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:230] "PPBP" "IGLL5" "VDAC3" "CD1C" ...
        .. .. .. .. .. ..$ : chr [1:39] "PC_1" "PC_2" "PC_3" "PC_4" ...
        .. .. .. ..@ feature.loadings.projected: num[0 , 0 ] 
        .. .. .. ..@ assay.used                : chr "RNA"
        .. .. .. ..@ global                    : logi FALSE
        .. .. .. ..@ stdev                     : num [1:39] 5.75 5.21 4.32 3.62 2.77 ...
        .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
        .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ] 
        .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ] 
        .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ] 
        .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ] 
        .. .. .. ..@ misc                      :List of 1
        .. .. .. .. ..$ total.variance: num 230
        .. .. .. ..@ key                       : chr "PC_"
        ..@ images      : list()
        ..@ project.name: chr "SeuratProject"
        ..@ misc        :List of 6
        .. ..$ gene_annotations:'data.frame':	230 obs. of  3 variables:
        .. .. ..$ input        : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. ..$ name         : chr [1:230] "SYMBOL - MS4A1" "SYMBOL - CD79B" "SYMBOL - CD79A" "SYMBOL - HLA-DRA" ...
        .. .. ..$ original_name: chr [1:230] "SYMBOL - MS4A1" "SYMBOL - CD79B" "SYMBOL - CD79A" "SYMBOL - HLA-DRA" ...
        .. ..$ color_pool      : chr [1:383] "#77aadd" "#ee8866" "#eedd88" "#ffaabb" ...
        .. ..$ ingestionDate   : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. ..$ gene_dispersion :'data.frame':	230 obs. of  5 variables:
        .. .. ..$ mean                 : num [1:230] 0.388 0.6 0.7 13.425 0.3 ...
        .. .. ..$ variance             : num [1:230] 1.025 1.281 4.365 725.463 0.871 ...
        .. .. ..$ variance.standardized: num [1:230] 0.898 0.473 1.086 0.973 1.356 ...
        .. .. ..$ SYMBOL               : chr [1:230] "SYMBOL - MS4A1" "SYMBOL - CD79B" "SYMBOL - CD79A" "SYMBOL - HLA-DRA" ...
        .. .. ..$ ENSEMBL              : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. ..$ active.reduction: chr "pca"
        .. ..$ numPCs          : num 2
        ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        .. ..$ : int [1:3] 5 0 2
        ..@ commands    :List of 4
        .. ..$ NormalizeData.RNA       :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "NormalizeData.RNA"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : chr "RNA"
        .. .. .. ..@ call.string: chr [1:2] "NormalizeData(scdata, normalization.method = normalization, " "    verbose = FALSE)"
        .. .. .. ..@ params     :List of 5
        .. .. .. .. ..$ assay               : chr "RNA"
        .. .. .. .. ..$ normalization.method: chr "LogNormalize"
        .. .. .. .. ..$ scale.factor        : num 10000
        .. .. .. .. ..$ margin              : num 1
        .. .. .. .. ..$ verbose             : logi FALSE
        .. ..$ FindVariableFeatures.RNA:Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "FindVariableFeatures.RNA"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : chr "RNA"
        .. .. .. ..@ call.string: chr [1:2] "FindVariableFeatures(Seurat::NormalizeData(scdata, normalization.method = normalization, " "    verbose = FALSE), assay = \"RNA\", nfeatures = nfeatures, verbose = FALSE)"
        .. .. .. ..@ params     :List of 12
        .. .. .. .. ..$ assay              : chr "RNA"
        .. .. .. .. ..$ selection.method   : chr "vst"
        .. .. .. .. ..$ loess.span         : num 0.3
        .. .. .. .. ..$ clip.max           : chr "auto"
        .. .. .. .. ..$ mean.function      :function (mat, display_progress)  
        .. .. .. .. .. ..- attr(*, "srcref")= 'srcref' int [1:8] 48 16 50 1 16 1 50 52
        .. .. .. .. .. .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x63a6d062c000> 
        .. .. .. .. ..$ dispersion.function:function (mat, display_progress)  
        .. .. .. .. .. ..- attr(*, "srcref")= 'srcref' int [1:8] 60 15 62 1 15 1 62 64
        .. .. .. .. .. .. ..- attr(*, "srcfile")=Classes 'srcfilealias', 'srcfile' <environment: 0x63a6d062c000> 
        .. .. .. .. ..$ num.bin            : num 20
        .. .. .. .. ..$ binning.method     : chr "equal_width"
        .. .. .. .. ..$ nfeatures          : num 1000
        .. .. .. .. ..$ mean.cutoff        : num [1:2] 0.1 8
        .. .. .. .. ..$ dispersion.cutoff  : num [1:2] 1 Inf
        .. .. .. .. ..$ verbose            : logi FALSE
        .. ..$ ScaleData.RNA           :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "ScaleData.RNA"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : chr "RNA"
        .. .. .. ..@ call.string: chr [1:3] "ScaleData(Seurat::FindVariableFeatures(Seurat::NormalizeData(scdata, " "    normalization.method = normalization, verbose = FALSE), assay = \"RNA\", " "    nfeatures = nfeatures, verbose = FALSE), verbose = FALSE)"
        .. .. .. ..@ params     :List of 10
        .. .. .. .. ..$ features          : chr [1:230] "PPBP" "IGLL5" "VDAC3" "CD1C" ...
        .. .. .. .. ..$ assay             : chr "RNA"
        .. .. .. .. ..$ model.use         : chr "linear"
        .. .. .. .. ..$ use.umi           : logi FALSE
        .. .. .. .. ..$ do.scale          : logi TRUE
        .. .. .. .. ..$ do.center         : logi TRUE
        .. .. .. .. ..$ scale.max         : num 10
        .. .. .. .. ..$ block.size        : num 1000
        .. .. .. .. ..$ min.cells.to.block: num 80
        .. .. .. .. ..$ verbose           : logi FALSE
        .. ..$ RunPCA.RNA              :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "RunPCA.RNA"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : chr "RNA"
        .. .. .. ..@ call.string: chr [1:2] "RunPCA(scdata, npcs = npcs_for_pca, features = Seurat::VariableFeatures(object = scdata), " "    verbose = FALSE)"
        .. .. .. ..@ params     :List of 11
        .. .. .. .. ..$ assay          : chr "RNA"
        .. .. .. .. ..$ features       : chr [1:230] "PPBP" "IGLL5" "VDAC3" "CD1C" ...
        .. .. .. .. ..$ npcs           : num 39
        .. .. .. .. ..$ rev.pca        : logi FALSE
        .. .. .. .. ..$ weight.by.var  : logi TRUE
        .. .. .. .. ..$ verbose        : logi FALSE
        .. .. .. .. ..$ ndims.print    : int [1:5] 1 2 3 4 5
        .. .. .. .. ..$ nfeatures.print: num 30
        .. .. .. .. ..$ reduction.name : chr "pca"
        .. .. .. .. ..$ reduction.key  : chr "PC_"
        .. .. .. .. ..$ seed.use       : num 42
        ..@ tools       : list()

