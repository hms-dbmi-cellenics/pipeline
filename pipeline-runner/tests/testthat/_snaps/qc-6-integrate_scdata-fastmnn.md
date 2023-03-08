# FastMNN works

    Code
      str(integrated_scdata)
    Output
      Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        ..@ assays      :List of 2
        .. ..$ RNA              :Formal class 'Assay' [package "SeuratObject"] with 8 slots
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
        .. .. .. ..@ scale.data   : num[0 , 0 ] 
        .. .. .. ..@ key          : Named chr "rna_"
        .. .. .. .. ..- attr(*, "names")= chr ""
        .. .. .. ..@ assay.orig   : NULL
        .. .. .. ..@ var.features : logi(0) 
        .. .. .. ..@ meta.features:'data.frame':	230 obs. of  0 variables
        .. .. .. ..@ misc         : list()
        .. ..$ mnn.reconstructed:Formal class 'Assay' [package "SeuratObject"] with 8 slots
        .. .. .. ..@ counts       : num[0 , 0 ] 
        .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. ..@ i       : int [1:18400] 0 1 2 3 4 5 6 7 8 9 ...
        .. .. .. .. .. ..@ p       : int [1:81] 0 230 460 690 920 1150 1380 1610 1840 2070 ...
        .. .. .. .. .. ..@ Dim     : int [1:2] 230 80
        .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. ..$ : chr [1:80] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. ..@ x       : num [1:18400] -0.0248 -0.0135 -0.0249 -0.0476 -0.0187 ...
        .. .. .. .. .. ..@ factors : list()
        .. .. .. ..@ scale.data   : num[0 , 0 ] 
        .. .. .. ..@ key          : chr "mnnreconstructed_"
        .. .. .. ..@ assay.orig   : NULL
        .. .. .. ..@ var.features : chr [1:230] "PPBP" "IGLL5" "VDAC3" "CD1C" ...
        .. .. .. ..@ meta.features:'data.frame':	230 obs. of  0 variables
        .. .. .. ..@ misc         : list()
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
        .. ..$ mnn:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
        .. .. .. ..@ cell.embeddings           : num [1:80, 1:5] 0.055 0.119 0.144 0.246 0.221 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:80] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. ..$ : chr [1:5] "mnn_1" "mnn_2" "mnn_3" "mnn_4" ...
        .. .. .. ..@ feature.loadings          : num [1:230, 1:5] -0.01134 -0.03001 -0.02309 -0.14153 -0.00647 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. ..$ : chr [1:5] "mnn_1" "mnn_2" "mnn_3" "mnn_4" ...
        .. .. .. ..@ feature.loadings.projected: num[0 , 0 ] 
        .. .. .. ..@ assay.used                : chr "RNA"
        .. .. .. ..@ global                    : logi FALSE
        .. .. .. ..@ stdev                     : num(0) 
        .. .. .. ..@ key                       : chr "mnn_"
        .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
        .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ] 
        .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ] 
        .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ] 
        .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ] 
        .. .. .. ..@ misc                      : list()
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
        .. ..$ numPCs          : num 5
        .. ..$ active.reduction: chr "mnn"
        ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        .. ..$ : int [1:3] 4 1 0
        ..@ commands    :List of 1
        .. ..$ SeuratWrappers..RunFastMNN.RNA:Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "SeuratWrappers::RunFastMNN.RNA"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : chr "RNA"
        .. .. .. ..@ call.string: chr [1:2] "SeuratWrappers::RunFastMNN(scdata_split, features = nfeatures, " "    d = npcs, get.variance = TRUE)"
        .. .. .. ..@ params     :List of 6
        .. .. .. .. ..$ assay              : chr "RNA"
        .. .. .. .. ..$ features           : chr [1:230] "PPBP" "IGLL5" "VDAC3" "CD1C" ...
        .. .. .. .. ..$ reduction.name     : chr "mnn"
        .. .. .. .. ..$ reduction.key      : chr "mnn_"
        .. .. .. .. ..$ reconstructed.assay: chr "mnn.reconstructed"
        .. .. .. .. ..$ verbose            : logi TRUE
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

