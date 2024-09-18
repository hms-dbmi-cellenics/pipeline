# misc slot is complete after Seurat V4 integration

    Code
      str(integrated_scdata)
    Output
      Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        ..@ assays      :List of 2
        .. ..$ RNA       :Formal class 'Assay5' [package "SeuratObject"] with 8 slots
        .. .. .. ..@ layers    :List of 3
        .. .. .. .. ..$ scale.data: num [1:230, 1:240] -0.476 1.47 -0.487 -1.34 -0.374 ...
        .. .. .. .. ..$ data      :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. ..@ i       : int [1:13368] 1 5 8 11 22 30 33 34 36 38 ...
        .. .. .. .. .. .. ..@ p       : int [1:241] 0 47 99 149 205 258 306 342 387 423 ...
        .. .. .. .. .. .. ..@ Dim     : int [1:2] 230 240
        .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..@ x       : num [1:13368] 4.97 4.97 6.06 4.97 4.97 ...
        .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. ..$ counts    :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. ..@ i       : int [1:13368] 1 5 8 11 22 30 33 34 36 38 ...
        .. .. .. .. .. .. ..@ p       : int [1:241] 0 47 99 149 205 258 306 342 387 423 ...
        .. .. .. .. .. .. ..@ Dim     : int [1:2] 230 240
        .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..@ x       : num [1:13368] 1 1 3 1 1 4 1 5 1 1 ...
        .. .. .. .. .. .. ..@ factors : list()
        .. .. .. ..@ cells     :Formal class 'LogMap' [package "SeuratObject"] with 1 slot
        .. .. .. .. .. ..@ .Data: logi [1:240, 1:3] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. .. ..$ : chr [1:240] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. ..$ : chr [1:3] "counts" "data" "scale.data"
        .. .. .. .. .. ..$ dim     : int [1:2] 240 3
        .. .. .. .. .. ..$ dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:240] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. .. ..$ : chr [1:3] "counts" "data" "scale.data"
        .. .. .. ..@ features  :Formal class 'LogMap' [package "SeuratObject"] with 1 slot
        .. .. .. .. .. ..@ .Data: logi [1:230, 1:3] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. ..$ : chr [1:3] "counts" "data" "scale.data"
        .. .. .. .. .. ..$ dim     : int [1:2] 230 3
        .. .. .. .. .. ..$ dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. ..$ : chr [1:3] "counts" "data" "scale.data"
        .. .. .. ..@ default   : int 1
        .. .. .. ..@ assay.orig: chr(0) 
        .. .. .. ..@ meta.data :'data.frame':	230 obs. of  14 variables:
        .. .. .. .. ..$ vf_vst_counts.1_mean                 : num [1:230] 0.508 0.675 0.858 11.025 0.367 ...
        .. .. .. .. ..$ vf_vst_counts.1_variance             : num [1:230] 1.29 1.41 5.16 516.86 1.01 ...
        .. .. .. .. ..$ vf_vst_counts.1_variance.expected    : num [1:230] 1.68 3.91 5.87 506.74 1.11 ...
        .. .. .. .. ..$ vf_vst_counts.1_variance.standardized: num [1:230] 0.772 0.362 0.88 1.02 0.908 ...
        .. .. .. .. ..$ vf_vst_counts.1_variable             : logi [1:230] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. .. ..$ vf_vst_counts.1_rank                 : int [1:230] 139 221 115 93 112 96 152 186 172 100 ...
        .. .. .. .. ..$ vf_vst_counts.2_mean                 : num [1:230] 0.267 0.525 0.542 15.825 0.233 ...
        .. .. .. .. ..$ vf_vst_counts.2_variance             : num [1:230] 0.718 1.125 3.477 916.347 0.718 ...
        .. .. .. .. ..$ vf_vst_counts.2_variance.expected    : num [1:230] 0.526 2.432 2.632 971.807 0.432 ...
        .. .. .. .. ..$ vf_vst_counts.2_variance.standardized: num [1:230] 1.364 0.463 1.321 0.943 1.662 ...
        .. .. .. .. ..$ vf_vst_counts.2_variable             : logi [1:230] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. .. ..$ vf_vst_counts.2_rank                 : int [1:230] 51 213 58 109 32 107 75 219 209 94 ...
        .. .. .. .. ..$ var.features                         : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. ..$ var.features.rank                    : int [1:230] 93 220 86 104 77 105 115 210 194 97 ...
        .. .. .. ..@ misc      : list()
        .. .. .. ..@ key       : chr "rna_"
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
        .. .. .. .. .. ..@ i       : int [1:28305] 14 27 32 45 57 59 63 67 93 104 ...
        .. .. .. .. .. ..@ p       : int [1:241] 0 47 99 149 205 258 306 342 387 423 ...
        .. .. .. .. .. ..@ Dim     : int [1:2] 230 240
        .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. .. .. .. ..$ : chr [1:240] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. ..@ x       : num [1:28305] 5.66 4.97 4.97 5.66 4.97 ...
        .. .. .. .. .. ..@ factors : list()
        .. .. .. ..@ scale.data   : num [1:230, 1:240] -0.548 -0.192 -0.378 -0.429 -0.355 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. .. .. ..$ : chr [1:240] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. ..@ assay.orig   : NULL
        .. .. .. ..@ var.features : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. ..@ meta.features:'data.frame':	230 obs. of  0 variables
        .. .. .. ..@ misc         : NULL
        .. .. .. ..@ key          : chr "integrated_"
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
        .. .. .. ..@ cell.embeddings           : num [1:240, 1:50] 3.27 3.69 2.67 3.72 3.07 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:240] "ATGCCAGAACGACT" "CATGGCCTGTGCAT" "GAACCTGATGAACC" "TGACTGGATTCTCA" ...
        .. .. .. .. .. ..$ : chr [1:50] "PC_1" "PC_2" "PC_3" "PC_4" ...
        .. .. .. ..@ feature.loadings          : num [1:230, 1:50] 0.02823 0.00596 0.04809 0.04561 0.05145 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. .. .. ..$ : chr [1:50] "PC_1" "PC_2" "PC_3" "PC_4" ...
        .. .. .. ..@ feature.loadings.projected: num[0 , 0 ] 
        .. .. .. ..@ assay.used                : chr "integrated"
        .. .. .. ..@ global                    : logi FALSE
        .. .. .. ..@ stdev                     : num [1:50] 5.82 5.25 4.34 3.64 2.75 ...
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
        .. ..$ gene_dispersion :'data.frame':	230 obs. of  6 variables:
        .. .. ..$ mean                 : num [1:230] 0.508 0.675 0.858 11.025 0.367 ...
        .. .. ..$ variance             : num [1:230] 1.29 1.41 5.16 516.86 1.01 ...
        .. .. ..$ variance.standardized: num [1:230] 1.68 3.91 5.87 506.74 1.11 ...
        .. .. ..$ NA                   : num [1:230] 0.772 0.362 0.88 1.02 0.908 ...
        .. .. ..$ SYMBOL               : chr [1:230] "SYMBOL - MS4A1" "SYMBOL - CD79B" "SYMBOL - CD79A" "SYMBOL - HLA-DRA" ...
        .. .. ..$ ENSEMBL              : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. ..$ active.reduction: chr "pca"
        .. ..$ numPCs          : num 50
        ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        .. ..$ : int [1:3] 5 0 2
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
        .. .. .. .. ..$ dims                : int [1:50] 1 2 3 4 5 6 7 8 9 10 ...
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
        .. .. .. .. ..$ dims                : int [1:50] 1 2 3 4 5 6 7 8 9 10 ...
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
        .. .. .. ..@ call.string: chr [1:2] "RunPCA(scdata, npcs = npcs_for_pca, features = Seurat::VariableFeatures(object = scdata), " "    verbose = FALSE)"
        .. .. .. ..@ params     :List of 11
        .. .. .. .. ..$ assay          : chr "integrated"
        .. .. .. .. ..$ features       : chr [1:230] "VDAC3" "IGLL5" "PF4" "PPBP" ...
        .. .. .. .. ..$ npcs           : num 50
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
        .. .. .. ..@ anchors           :'data.frame':	640 obs. of  5 variables:
        .. .. .. .. ..$ cell1   : num [1:640] 1 2 2 3 3 4 5 6 7 7 ...
        .. .. .. .. ..$ cell2   : num [1:640] 41 42 45 43 49 44 45 46 47 46 ...
        .. .. .. .. ..$ score   : num [1:640] 0.4 0.6 0.2 0.6 0.35 0.5 0.5 0.75 0.6 0.4 ...
        .. .. .. .. ..$ dataset1: int [1:640] 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..$ dataset2: int [1:640] 2 2 2 2 2 2 2 2 2 2 ...
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
       $ gene_dispersion :'data.frame':	230 obs. of  6 variables:
        ..$ mean                 : num [1:230] 0.508 0.675 0.858 11.025 0.367 ...
        ..$ variance             : num [1:230] 1.29 1.41 5.16 516.86 1.01 ...
        ..$ variance.standardized: num [1:230] 1.68 3.91 5.87 506.74 1.11 ...
        ..$ NA                   : num [1:230] 0.772 0.362 0.88 1.02 0.908 ...
        ..$ SYMBOL               : chr [1:230] "SYMBOL - MS4A1" "SYMBOL - CD79B" "SYMBOL - CD79A" "SYMBOL - HLA-DRA" ...
        ..$ ENSEMBL              : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
       $ active.reduction: chr "pca"
       $ numPCs          : num 50

