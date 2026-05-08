# FastMNN works with or without bpcells

    Code
      str(integrated_scdata)
    Output
      Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        ..@ assays      :List of 1
        .. ..$ RNA:Formal class 'Assay5' [package "SeuratObject"] with 8 slots
        .. .. .. ..@ layers    :List of 4
        .. .. .. .. ..$ counts.1:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. ..@ i       : int [1:2060] 1 5 8 11 22 30 33 34 36 38 ...
        .. .. .. .. .. .. ..@ p       : int [1:41] 0 47 99 149 205 258 306 342 387 423 ...
        .. .. .. .. .. .. ..@ Dim     : int [1:2] 230 40
        .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..@ x       : num [1:2060] 1 1 3 1 1 4 1 5 1 1 ...
        .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. ..$ counts.2:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. ..@ i       : int [1:2396] 22 31 42 43 46 47 48 52 55 60 ...
        .. .. .. .. .. .. ..@ p       : int [1:41] 0 60 121 165 197 244 284 332 387 439 ...
        .. .. .. .. .. .. ..@ Dim     : int [1:2] 230 40
        .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..@ x       : num [1:2396] 7 1 2 2 1 1 1 2 2 1 ...
        .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. ..$ data.1  :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. ..@ i       : int [1:2060] 1 5 8 11 22 30 33 34 36 38 ...
        .. .. .. .. .. .. ..@ p       : int [1:41] 0 47 99 149 205 258 306 342 387 423 ...
        .. .. .. .. .. .. ..@ Dim     : int [1:2] 230 40
        .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..@ x       : num [1:2060] 4.97 4.97 6.06 4.97 4.97 ...
        .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. ..$ data.2  :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. ..@ i       : int [1:2396] 22 31 42 43 46 47 48 52 55 60 ...
        .. .. .. .. .. .. ..@ p       : int [1:41] 0 60 121 165 197 244 284 332 387 439 ...
        .. .. .. .. .. .. ..@ Dim     : int [1:2] 230 40
        .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..@ x       : num [1:2396] 6.11 4.18 4.86 4.86 4.18 ...
        .. .. .. .. .. .. ..@ factors : list()
        .. .. .. ..@ cells     :Formal class 'LogMap' [package "SeuratObject"] with 1 slot
        .. .. .. .. .. ..@ .Data: logi [1:80, 1:4] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. .. ..$ : chr [1:80] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. ..$ : chr [1:4] "counts.1" "counts.2" "data.1" "data.2"
        .. .. .. .. .. ..$ dim     : int [1:2] 80 4
        .. .. .. .. .. ..$ dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:80] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. ..$ : chr [1:4] "counts.1" "counts.2" "data.1" "data.2"
        .. .. .. ..@ features  :Formal class 'LogMap' [package "SeuratObject"] with 1 slot
        .. .. .. .. .. ..@ .Data: logi [1:230, 1:4] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. ..$ : chr [1:4] "counts.1" "counts.2" "data.1" "data.2"
        .. .. .. .. .. ..$ dim     : int [1:2] 230 4
        .. .. .. .. .. ..$ dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. ..$ : chr [1:4] "counts.1" "counts.2" "data.1" "data.2"
        .. .. .. ..@ default   : int 2
        .. .. .. ..@ assay.orig: chr(0) 
        .. .. .. ..@ meta.data :'data.frame':	230 obs. of  14 variables:
        .. .. .. .. ..$ vf_vst_counts.1_mean                 : num [1:230] 0.75 0.825 1.175 6.225 0.5 ...
        .. .. .. .. ..$ vf_vst_counts.1_variance             : num [1:230] 1.78 1.69 6.76 72.13 1.28 ...
        .. .. .. .. ..$ vf_vst_counts.1_variance.expected    : num [1:230] 2.04 2.951 6.299 126.599 0.984 ...
        .. .. .. .. ..$ vf_vst_counts.1_variance.standardized: num [1:230] 0.874 0.572 1.074 0.57 1.303 ...
        .. .. .. .. ..$ vf_vst_counts.1_variable             : logi [1:230] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. .. ..$ vf_vst_counts.1_rank                 : int [1:230] 122 195 71 197 44 123 101 100 185 102 ...
        .. .. .. .. ..$ vf_vst_counts.2_mean                 : num [1:230] 0.025 0.375 0.225 20.625 0.1 ...
        .. .. .. .. ..$ vf_vst_counts.2_variance             : num [1:230] 0.025 0.804 1.615 1291.061 0.4 ...
        .. .. .. .. ..$ vf_vst_counts.2_variance.expected    : num [1:230] 2.47e-02 6.88e-01 3.14e-01 1.40e+03 1.28e-01 ...
        .. .. .. .. ..$ vf_vst_counts.2_variance.standardized: num [1:230] 1.011 1.169 1.232 0.922 1.104 ...
        .. .. .. .. ..$ vf_vst_counts.2_variable             : logi [1:230] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. .. ..$ vf_vst_counts.2_rank                 : int [1:230] 87 52 38 106 67 102 110 193 195 217 ...
        .. .. .. .. ..$ var.features                         : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. ..$ var.features.rank                    : int [1:230] 88 137 20 182 24 107 92 176 220 190 ...
        .. .. .. ..@ misc      : list()
        .. .. .. ..@ key       : chr "rna_"
        ..@ meta.data   :'data.frame':	80 obs. of  9 variables:
        .. ..$ orig.ident    : chr [1:80] "SeuratProject" "SeuratProject" "SeuratProject" "SeuratProject" ...
        .. ..$ nCount_RNA    : num [1:80] 70 85 87 127 173 70 64 72 52 100 ...
        .. ..$ nFeature_RNA  : int [1:80] 47 52 50 56 53 48 36 45 36 41 ...
        .. ..$ cells_id      : int [1:80] 0 1 2 3 4 5 6 7 8 9 ...
        .. ..$ samples       : chr [1:80] "123abc" "123abc" "123abc" "123abc" ...
        .. ..$ doublet_scores: num [1:80] 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 ...
        .. ..$ doublet_class : chr [1:80] "singlet" "singlet" "singlet" "singlet" ...
        .. ..$ percent.mt    : num [1:80] 5.44 5.77 7.56 6.07 6.13 ...
        .. ..$ emptyDrops_FDR: num [1:80] 0.009 0.009 0.009 0.009 0.009 0.009 0.009 0.009 0.009 0.009 ...
        ..@ active.assay: chr "RNA"
        ..@ active.ident: Factor w/ 1 level "SeuratProject": 1 1 1 1 1 1 1 1 1 1 ...
        .. ..- attr(*, "names")= chr [1:80] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        ..@ graphs      : list()
        ..@ neighbors   : list()
        ..@ reductions  :List of 1
        .. ..$ mnn:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
        .. .. .. ..@ cell.embeddings           : num [1:80, 1:5] -0.055 -0.119 -0.144 -0.246 -0.221 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:80] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. ..$ : chr [1:5] "mnn_1" "mnn_2" "mnn_3" "mnn_4" ...
        .. .. .. ..@ feature.loadings          : num [1:230, 1:5] 0.01134 0.03001 0.02309 0.14153 0.00647 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. ..$ : chr [1:5] "mnn_1" "mnn_2" "mnn_3" "mnn_4" ...
        .. .. .. ..@ feature.loadings.projected: num[0 , 0 ] 
        .. .. .. ..@ assay.used                : chr "RNA"
        .. .. .. ..@ global                    : logi FALSE
        .. .. .. ..@ stdev                     : num(0) 
        .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
        .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ] 
        .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ] 
        .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ] 
        .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ] 
        .. .. .. ..@ misc                      : list()
        .. .. .. ..@ key                       : chr "mnn_"
        ..@ images      : list()
        ..@ project.name: chr "SeuratProject"
        ..@ misc        :List of 6
        .. ..$ gene_annotations:'data.frame':	230 obs. of  3 variables:
        .. .. ..$ input        : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. ..$ name         : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. ..$ original_name: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. ..$ color_pool      : chr [1:383] "#77aadd" "#ee8866" "#eedd88" "#ffaabb" ...
        .. ..$ ingestionDate   : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. ..$ gene_dispersion :'data.frame':	230 obs. of  6 variables:
        .. .. ..$ mean                 : num [1:230] 0.75 0.825 1.175 6.225 0.5 ...
        .. .. ..$ variance             : num [1:230] 1.78 1.69 6.76 72.13 1.28 ...
        .. .. ..$ variance.expected    : num [1:230] 2.04 2.951 6.299 126.599 0.984 ...
        .. .. ..$ variance.standardized: num [1:230] 0.874 0.572 1.074 0.57 1.303 ...
        .. .. ..$ SYMBOL               : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. ..$ ENSEMBL              : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. ..$ numPCs          : num 5
        .. ..$ active.reduction: chr "mnn"
        ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        .. ..$ : int [1:3] 5 3 0
        ..@ commands    :List of 2
        .. ..$ NormalizeData.RNA       :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "NormalizeData.RNA"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : chr "RNA"
        .. .. .. ..@ call.string: chr [1:2] "NormalizeData(scdata, normalization.method = \"LogNormalize\", " "    verbose = FALSE)"
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
        .. .. .. ..@ call.string: chr [1:2] "FindVariableFeatures(Seurat::NormalizeData(scdata, normalization.method = \"LogNormalize\", " "    verbose = FALSE), nfeatures = nfeatures, verbose = FALSE)"
        .. .. .. ..@ params     :List of 10
        .. .. .. .. ..$ assay            : chr "RNA"
        .. .. .. .. ..$ selection.method : chr "vst"
        .. .. .. .. ..$ loess.span       : num 0.3
        .. .. .. .. ..$ clip.max         : chr "auto"
        .. .. .. .. ..$ num.bin          : num 20
        .. .. .. .. ..$ binning.method   : chr "equal_width"
        .. .. .. .. ..$ nfeatures        : num 1000
        .. .. .. .. ..$ mean.cutoff      : num [1:2] 0.1 8
        .. .. .. .. ..$ dispersion.cutoff: num [1:2] 1 Inf
        .. .. .. .. ..$ verbose          : logi FALSE
        ..@ tools       :List of 1
        .. ..$ SeuratWrappers::RunFastMNN:List of 2
        .. .. ..$ merge.info:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. .. ..@ rownames       : NULL
        .. .. .. .. ..@ nrows          : int 1
        .. .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. .. ..@ elementMetadata: NULL
        .. .. .. .. ..@ metadata       : list()
        .. .. .. .. ..@ listData       :List of 6
        .. .. .. .. .. ..$ left      :Formal class 'SimpleList' [package "S4Vectors"] with 4 slots
        .. .. .. .. .. .. .. ..@ listData       :List of 1
        .. .. .. .. .. .. .. .. ..$ : chr "123abc"
        .. .. .. .. .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. .. .. .. .. ..@ elementMetadata: NULL
        .. .. .. .. .. .. .. ..@ metadata       : list()
        .. .. .. .. .. ..$ right     :Formal class 'SimpleList' [package "S4Vectors"] with 4 slots
        .. .. .. .. .. .. .. ..@ listData       :List of 1
        .. .. .. .. .. .. .. .. ..$ : chr "123def"
        .. .. .. .. .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. .. .. .. .. ..@ elementMetadata: NULL
        .. .. .. .. .. .. .. ..@ metadata       : list()
        .. .. .. .. .. ..$ pairs     :Formal class 'SimpleDFrameList' [package "IRanges"] with 4 slots
        .. .. .. .. .. .. .. ..@ elementType    : chr "DFrame"
        .. .. .. .. .. .. .. ..@ elementMetadata: NULL
        .. .. .. .. .. .. .. ..@ metadata       : list()
        .. .. .. .. .. .. .. ..@ listData       :List of 1
        .. .. .. .. .. .. .. .. ..$ :Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. .. .. .. .. .. .. .. ..@ rownames       : NULL
        .. .. .. .. .. .. .. .. .. .. ..@ nrows          : int 558
        .. .. .. .. .. .. .. .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. .. .. .. .. .. .. .. ..@ elementMetadata: NULL
        .. .. .. .. .. .. .. .. .. .. ..@ metadata       : list()
        .. .. .. .. .. .. .. .. .. .. ..@ listData       :List of 2
        .. .. .. .. .. .. .. .. .. .. .. ..$ left : int [1:558] 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. .. .. .. .. .. .. .. ..$ right: int [1:558] 50 45 48 43 42 41 54 60 46 47 ...
        .. .. .. .. .. ..$ batch.size: num 0.661
        .. .. .. .. .. ..$ skipped   : logi FALSE
        .. .. .. .. .. ..$ lost.var  : num [1, 1:2] 0.1163 0.0904
        .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. ..$ : chr [1:2] "123abc" "123def"
        .. .. ..$ pca.info  :List of 3
        .. .. .. ..$ centers      : Named num [1:230] 0.0239 0.0393 0.026 0.1118 0.0173 ...
        .. .. .. .. ..- attr(*, "names")= chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. ..$ var.explained: num [1:5] 0.1119 0.0939 0.0616 0.0394 0.019
        .. .. .. ..$ var.total    : num 0.652

