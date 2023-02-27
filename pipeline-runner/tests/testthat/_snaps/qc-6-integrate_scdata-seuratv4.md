# misc slot is complete after Seurat V4 integration

    Code
      str(integrated_scdata)
    Output
      Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        ..@ assays      :List of 2
        .. ..$ RNA       :Formal class 'Assay' [package "SeuratObject"] with 8 slots
        .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. ..@ i       : int [1:13368] 1 5 8 11 22 30 33 34 36 38 ...
        .. .. .. .. .. ..@ p       : int [1:241] 0 47 99 149 205 258 306 342 387 423 ...
        .. .. .. .. .. ..@ Dim     : int [1:2] 230 240
        .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. ..$ : chr [1:240] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. ..@ x       : num [1:13368] 1 1 3 1 1 4 1 5 1 1 ...
        .. .. .. .. .. ..@ factors : list()
        .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. ..@ i       : int [1:13368] 1 5 8 11 22 30 33 34 36 38 ...
        .. .. .. .. .. ..@ p       : int [1:241] 0 47 99 149 205 258 306 342 387 423 ...
        .. .. .. .. .. ..@ Dim     : int [1:2] 230 240
        .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. ..$ : chr [1:240] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. ..@ x       : num [1:13368] 4.97 4.97 6.06 4.97 4.97 ...
        .. .. .. .. .. ..@ factors : list()
        .. .. .. ..@ scale.data   : num[0 , 0 ] 
        .. .. .. ..@ key          : Named chr "rna_"
        .. .. .. .. ..- attr(*, "names")= chr ""
        .. .. .. ..@ assay.orig   : NULL
        .. .. .. ..@ var.features : chr [1:230] "IGLL5" "VDAC3" "TSC22D1" "PTCRA" ...
        .. .. .. ..@ meta.features:'data.frame':	230 obs. of  5 variables:
        .. .. .. .. ..$ vst.mean                 : num [1:230] 0.388 0.6 0.7 13.425 0.3 ...
        .. .. .. .. ..$ vst.variance             : num [1:230] 1.017 1.27 4.328 719.392 0.864 ...
        .. .. .. .. ..$ vst.variance.expected    : num [1:230] 1.132 2.685 3.984 739.287 0.637 ...
        .. .. .. .. ..$ vst.variance.standardized: num [1:230] 0.898 0.473 1.086 0.973 1.356 ...
        .. .. .. .. ..$ vst.variable             : logi [1:230] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. ..@ misc         : list()
        .. ..$ integrated:Formal class 'Assay' [package "SeuratObject"] with 8 slots
        .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. ..@ i       : int(0) 
        .. .. .. .. .. ..@ p       : int 0
        .. .. .. .. .. ..@ Dim     : int [1:2] 0 0
        .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. ..@ x       : num(0) 
        .. .. .. .. .. ..@ factors : list()
        .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. ..@ i       : int [1:28387] 14 27 32 45 57 59 63 67 93 104 ...
        .. .. .. .. .. ..@ p       : int [1:241] 0 47 99 149 205 258 306 342 387 423 ...
        .. .. .. .. .. ..@ Dim     : int [1:2] 230 240
        .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. .. .. .. ..$ : chr [1:240] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. ..@ x       : num [1:28387] 5.66 4.97 4.97 5.66 4.97 ...
        .. .. .. .. .. ..@ factors : list()
        .. .. .. ..@ scale.data   : num [1:230, 1:240] -0.542 -0.194 -0.37 -0.431 -0.349 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. .. .. ..$ : chr [1:240] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. ..@ key          : chr "integrated_"
        .. .. .. ..@ assay.orig   : NULL
        .. .. .. ..@ var.features : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. ..@ meta.features:'data.frame':	230 obs. of  0 variables
        .. .. .. ..@ misc         : NULL
        ..@ meta.data   :'data.frame':	240 obs. of  5 variables:
        .. ..$ orig.ident  : chr [1:240] "SeuratProject" "SeuratProject" "SeuratProject" "SeuratProject" ...
        .. ..$ nCount_RNA  : num [1:240] 70 85 87 127 173 70 64 72 52 100 ...
        .. ..$ nFeature_RNA: int [1:240] 47 52 50 56 53 48 36 45 36 41 ...
        .. ..$ samples     : chr [1:240] "123abc" "123abc" "123abc" "123abc" ...
        .. ..$ cells_id    : int [1:240] 0 1 2 3 4 5 6 7 8 9 ...
        ..@ active.assay: chr "integrated"
        ..@ active.ident: Factor w/ 1 level "SeuratProject": 1 1 1 1 1 1 1 1 1 1 ...
        .. ..- attr(*, "names")= chr [1:240] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        ..@ graphs      : list()
        ..@ neighbors   : list()
        ..@ reductions  :List of 1
        .. ..$ pca:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
        .. .. .. ..@ cell.embeddings           : num [1:240, 1:30] -3.13 -3.64 -2.5 -3.65 -2.9 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:240] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. ..$ : chr [1:30] "PC_1" "PC_2" "PC_3" "PC_4" ...
        .. .. .. ..@ feature.loadings          : num [1:230, 1:30] -0.03297 -0.00653 -0.05748 -0.05327 -0.05767 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. .. .. ..$ : chr [1:30] "PC_1" "PC_2" "PC_3" "PC_4" ...
        .. .. .. ..@ feature.loadings.projected: num[0 , 0 ] 
        .. .. .. ..@ assay.used                : chr "integrated"
        .. .. .. ..@ global                    : logi FALSE
        .. .. .. ..@ stdev                     : num [1:30] 5.77 5.22 4.35 3.65 2.75 ...
        .. .. .. ..@ key                       : chr "PC_"
        .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
        .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ] 
        .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ] 
        .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ] 
        .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ] 
        .. .. .. ..@ misc                      :List of 1
        .. .. .. .. ..$ total.variance: num 230
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
        .. .. ..$ variance             : num [1:230] 1.017 1.27 4.328 719.392 0.864 ...
        .. .. ..$ variance.standardized: num [1:230] 0.898 0.473 1.086 0.973 1.356 ...
        .. .. ..$ SYMBOL               : chr [1:230] "SYMBOL - MS4A1" "SYMBOL - CD79B" "SYMBOL - CD79A" "SYMBOL - HLA-DRA" ...
        .. .. ..$ ENSEMBL              : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. ..$ active.reduction: chr "pca"
        .. ..$ numPCs          : num 30
        ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        .. ..$ : int [1:3] 4 1 0
        ..@ commands    :List of 5
        .. ..$ FindIntegrationAnchors  :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "Seurat::FindIntegrationAnchors"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : Named chr [1:2] "RNA" "RNA"
        .. .. .. .. ..- attr(*, "names")= chr [1:2] "123abc" "123def"
        .. .. .. ..@ call.string: chr [1:3] "Seurat::FindIntegrationAnchors(object.list = scdata_list, dims = 1:npcs, " "    k.filter = k.filter, normalization.method = normalization, " "    verbose = TRUE, reduction = reduction)"
        .. .. .. ..@ params     :List of 15
        .. .. .. .. ..$ assay               : Named chr [1:2] "RNA" "RNA"
        .. .. .. .. .. ..- attr(*, "names")= chr [1:2] "123abc" "123def"
        .. .. .. .. ..$ anchor.features     : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. .. ..$ scale               : logi TRUE
        .. .. .. .. ..$ normalization.method: chr "LogNormalize"
        .. .. .. .. ..$ reduction           : chr "pca"
        .. .. .. .. ..$ l2.norm             : logi TRUE
        .. .. .. .. ..$ dims                : int [1:30] 1 2 3 4 5 6 7 8 9 10 ...
        .. .. .. .. ..$ k.anchor            : num 5
        .. .. .. .. ..$ k.filter            : logi NA
        .. .. .. .. ..$ k.score             : num 30
        .. .. .. .. ..$ max.features        : num 200
        .. .. .. .. ..$ nn.method           : chr "annoy"
        .. .. .. .. ..$ n.trees             : num 50
        .. .. .. .. ..$ eps                 : num 0
        .. .. .. .. ..$ verbose             : logi TRUE
        .. ..$ withCallingHandlers     :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "withCallingHandlers"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : NULL
        .. .. .. ..@ call.string: chr [1:2] "withCallingHandlers(expr, warning = function(w) if (inherits(w, " "    classes)) tryInvokeRestart(\"muffleWarning\"))"
        .. .. .. ..@ params     :List of 9
        .. .. .. .. ..$ new.assay.name      : chr "integrated"
        .. .. .. .. ..$ normalization.method: chr "LogNormalize"
        .. .. .. .. ..$ features            : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. .. ..$ dims                : int [1:30] 1 2 3 4 5 6 7 8 9 10 ...
        .. .. .. .. ..$ k.weight            : num 59
        .. .. .. .. ..$ sd.weight           : num 1
        .. .. .. .. ..$ preserve.order      : logi FALSE
        .. .. .. .. ..$ eps                 : num 0
        .. .. .. .. ..$ verbose             : logi TRUE
        .. ..$ FindVariableFeatures.RNA:Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "FindVariableFeatures.RNA"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : chr "RNA"
        .. .. .. ..@ call.string: chr [1:2] "FindVariableFeatures(scdata, assay = \"RNA\", nfeatures = nfeatures, " "    verbose = FALSE)"
        .. .. .. ..@ params     :List of 12
        .. .. .. .. ..$ assay              : chr "RNA"
        .. .. .. .. ..$ selection.method   : chr "vst"
        .. .. .. .. ..$ loess.span         : num 0.3
        .. .. .. .. ..$ clip.max           : chr "auto"
        .. .. .. .. ..$ mean.function      :function (mat, display_progress)  
        .. .. .. .. ..$ dispersion.function:function (mat, display_progress)  
        .. .. .. .. ..$ num.bin            : num 20
        .. .. .. .. ..$ binning.method     : chr "equal_width"
        .. .. .. .. ..$ nfeatures          : num 1000
        .. .. .. .. ..$ mean.cutoff        : num [1:2] 0.1 8
        .. .. .. .. ..$ dispersion.cutoff  : num [1:2] 1 Inf
        .. .. .. .. ..$ verbose            : logi FALSE
        .. ..$ ScaleData.integrated    :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "ScaleData.integrated"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : chr "integrated"
        .. .. .. ..@ call.string: chr "ScaleData(scdata, verbose = FALSE)"
        .. .. .. ..@ params     :List of 10
        .. .. .. .. ..$ features          : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. .. ..$ assay             : chr "integrated"
        .. .. .. .. ..$ model.use         : chr "linear"
        .. .. .. .. ..$ use.umi           : logi FALSE
        .. .. .. .. ..$ do.scale          : logi TRUE
        .. .. .. .. ..$ do.center         : logi TRUE
        .. .. .. .. ..$ scale.max         : num 10
        .. .. .. .. ..$ block.size        : num 1000
        .. .. .. .. ..$ min.cells.to.block: num 240
        .. .. .. .. ..$ verbose           : logi FALSE
        .. ..$ RunPCA.integrated       :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "RunPCA.integrated"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : chr "integrated"
        .. .. .. ..@ call.string: chr [1:2] "RunPCA(scdata, npcs = npcs, features = Seurat::VariableFeatures(object = scdata), " "    verbose = FALSE)"
        .. .. .. ..@ params     :List of 11
        .. .. .. .. ..$ assay          : chr "integrated"
        .. .. .. .. ..$ features       : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. .. ..$ npcs           : num 30
        .. .. .. .. ..$ rev.pca        : logi FALSE
        .. .. .. .. ..$ weight.by.var  : logi TRUE
        .. .. .. .. ..$ verbose        : logi FALSE
        .. .. .. .. ..$ ndims.print    : int [1:5] 1 2 3 4 5
        .. .. .. .. ..$ nfeatures.print: num 30
        .. .. .. .. ..$ reduction.name : chr "pca"
        .. .. .. .. ..$ reduction.key  : chr "PC_"
        .. .. .. .. ..$ seed.use       : num 42
        ..@ tools       :List of 1
        .. ..$ Integration:Formal class 'IntegrationData' [package "Seurat"] with 7 slots
        .. .. .. ..@ neighbors         : NULL
        .. .. .. ..@ weights           : NULL
        .. .. .. ..@ integration.matrix: NULL
        .. .. .. ..@ anchors           :'data.frame':	626 obs. of  5 variables:
        .. .. .. .. ..$ cell1   : num [1:626] 1 2 3 3 3 4 5 5 5 6 ...
        .. .. .. .. ..$ cell2   : num [1:626] 41 42 43 49 45 44 45 42 43 46 ...
        .. .. .. .. ..$ score   : num [1:626] 0.579 0.632 0.842 0 0.316 ...
        .. .. .. .. ..$ dataset1: int [1:626] 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..$ dataset2: int [1:626] 2 2 2 2 2 2 2 2 2 2 ...
        .. .. .. ..@ offsets           : NULL
        .. .. .. ..@ objects.ncell     : NULL
        .. .. .. ..@ sample.tree       : int [1, 1:2] -1 -2

---

    Code
      str(integrated_scdata@misc)
    Output
      List of 6
       $ gene_annotations:'data.frame':	230 obs. of  3 variables:
        ..$ input        : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        ..$ name         : chr [1:230] "SYMBOL - MS4A1" "SYMBOL - CD79B" "SYMBOL - CD79A" "SYMBOL - HLA-DRA" ...
        ..$ original_name: chr [1:230] "SYMBOL - MS4A1" "SYMBOL - CD79B" "SYMBOL - CD79A" "SYMBOL - HLA-DRA" ...
       $ color_pool      : chr [1:383] "#77aadd" "#ee8866" "#eedd88" "#ffaabb" ...
       $ ingestionDate   : POSIXct[1:1], format: "1991-12-19 05:23:00"
       $ gene_dispersion :'data.frame':	230 obs. of  5 variables:
        ..$ mean                 : num [1:230] 0.388 0.6 0.7 13.425 0.3 ...
        ..$ variance             : num [1:230] 1.017 1.27 4.328 719.392 0.864 ...
        ..$ variance.standardized: num [1:230] 0.898 0.473 1.086 0.973 1.356 ...
        ..$ SYMBOL               : chr [1:230] "SYMBOL - MS4A1" "SYMBOL - CD79B" "SYMBOL - CD79A" "SYMBOL - HLA-DRA" ...
        ..$ ENSEMBL              : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
       $ active.reduction: chr "pca"
       $ numPCs          : num 30

