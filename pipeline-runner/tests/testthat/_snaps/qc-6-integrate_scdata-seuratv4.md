# misc slot is complete after Seurat V4 integration

    Code
      str(integrated_scdata@misc)
    Output
      List of 6
       $ gene_annotations:'data.frame':	230 obs. of  3 variables:
        ..$ input        : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        ..$ name         : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        ..$ original_name: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
       $ color_pool      : chr [1:383] "#77aadd" "#ee8866" "#eedd88" "#ffaabb" ...
       $ ingestionDate   : POSIXct[1:1], format: "1991-12-19 05:23:00"
       $ gene_dispersion :'data.frame':	230 obs. of  6 variables:
        ..$ mean                 : num [1:230] 0.508 0.675 0.858 11.025 0.367 ...
        ..$ variance             : num [1:230] 1.29 1.41 5.16 516.86 1.01 ...
        ..$ variance.expected    : num [1:230] 1.68 3.91 5.87 506.74 1.15 ...
        ..$ variance.standardized: num [1:230] 0.772 0.362 0.88 1.02 0.874 ...
        ..$ SYMBOL               : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        ..$ ENSEMBL              : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
       $ numPCs          : num 50
       $ active.reduction: chr "integrated.rpca"

---

    Code
      str(integrated_scdata)
    Output
      Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        ..@ assays      :List of 1
        .. ..$ RNA:Formal class 'Assay5' [package "SeuratObject"] with 8 slots
        .. .. .. ..@ layers    :List of 3
        .. .. .. .. ..$ data      :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. ..@ matrix   :Formal class 'ColBindMatrices' [package "BPCells"] with 5 slots
        .. .. .. .. .. .. .. .. ..@ matrix_list:List of 2
        .. .. .. .. .. .. .. .. .. ..$ :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'TransformLog1p' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'TransformScaleShift' [package "BPCells"] with 8 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ active_transforms: logi [1:3, 1:2] FALSE TRUE TRUE FALSE FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:3] "row" "col" "global"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:2] "scale" "shift"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix           :Formal class 'ConvertMatrixType' [package "BPCells"] with 5 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'MatrixSubset' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'MatrixDir' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dir        : chr "/tmp/RtmpNY4Wop/matrix_dir"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ compressed : logi TRUE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ buffer_size: int 8192
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ type       : chr "uint32_t"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim        : int [1:2] 230 240
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose  : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames   :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_selection: int [1:230] 1 2 3 4 5 6 7 8 9 10 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_selection: int [1:120] 1 2 3 4 5 6 7 8 9 10 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ zero_dims    : logi [1:2] FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ cells   : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ type     : chr "double"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_params       : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_params       : num [1:2, 1:120] 0.0143 0 0.0118 0 0.0115 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ global_params    : num [1:2] 10000 0
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim              : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose        : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames         :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_params   : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_params   : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ global_params: num(0) 
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. ..$ cells   : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. ..$ :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'TransformLog1p' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'TransformScaleShift' [package "BPCells"] with 8 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ active_transforms: logi [1:3, 1:2] FALSE TRUE TRUE FALSE FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:3] "row" "col" "global"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:2] "scale" "shift"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix           :Formal class 'ConvertMatrixType' [package "BPCells"] with 5 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'MatrixSubset' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'MatrixDir' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dir        : chr "/tmp/RtmpNY4Wop/matrix_dir"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ compressed : logi TRUE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ buffer_size: int 8192
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ type       : chr "uint32_t"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim        : int [1:2] 230 240
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose  : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames   :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_selection: int [1:230] 1 2 3 4 5 6 7 8 9 10 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_selection: int [1:120] 121 122 123 124 125 126 127 128 129 130 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ zero_dims    : logi [1:2] FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ cells   : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ type     : chr "double"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_params       : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_params       : num [1:2, 1:120] 0.00641 0 0.00719 0 0.00926 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ global_params    : num [1:2] 10000 0
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim              : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose        : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames         :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_params   : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_params   : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ global_params: num(0) 
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. ..$ cells   : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. ..@ threads    : int 0
        .. .. .. .. .. .. .. .. ..@ dim        : int [1:2] 230 240
        .. .. .. .. .. .. .. .. ..@ transpose  : logi FALSE
        .. .. .. .. .. .. .. .. ..@ dimnames   :List of 2
        .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..@ dim      : int [1:2] 230 240
        .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. ..$ cells   : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. ..$ counts    :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. ..@ matrix   :Formal class 'ColBindMatrices' [package "BPCells"] with 5 slots
        .. .. .. .. .. .. .. .. ..@ matrix_list:List of 2
        .. .. .. .. .. .. .. .. .. ..$ :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'MatrixSubset' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'MatrixDir' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dir        : chr "/tmp/RtmpNY4Wop/matrix_dir"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ compressed : logi TRUE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ buffer_size: int 8192
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ type       : chr "uint32_t"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim        : int [1:2] 230 240
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose  : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames   :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_selection: int [1:230] 1 2 3 4 5 6 7 8 9 10 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_selection: int [1:120] 1 2 3 4 5 6 7 8 9 10 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ zero_dims    : logi [1:2] FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. ..$ cells   : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. ..$ :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'MatrixSubset' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'MatrixDir' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dir        : chr "/tmp/RtmpNY4Wop/matrix_dir"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ compressed : logi TRUE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ buffer_size: int 8192
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ type       : chr "uint32_t"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim        : int [1:2] 230 240
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose  : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames   :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_selection: int [1:230] 1 2 3 4 5 6 7 8 9 10 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_selection: int [1:120] 121 122 123 124 125 126 127 128 129 130 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ zero_dims    : logi [1:2] FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. ..$ cells   : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. ..@ threads    : int 0
        .. .. .. .. .. .. .. .. ..@ dim        : int [1:2] 230 240
        .. .. .. .. .. .. .. .. ..@ transpose  : logi FALSE
        .. .. .. .. .. .. .. .. ..@ dimnames   :List of 2
        .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..@ dim      : int [1:2] 230 240
        .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. ..$ cells   : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. ..$ scale.data:Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. ..@ matrix   :Formal class 'MatrixSubset' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'TransformScaleShift' [package "BPCells"] with 8 slots
        .. .. .. .. .. .. .. .. .. .. ..@ active_transforms: logi [1:3, 1:2] TRUE FALSE FALSE TRUE FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:3] "row" "col" "global"
        .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:2] "scale" "shift"
        .. .. .. .. .. .. .. .. .. .. ..@ matrix           :Formal class 'TransformMinByRow' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'MatrixSubset' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'ColBindMatrices' [package "BPCells"] with 5 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix_list:List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'TransformLog1p' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'TransformScaleShift' [package "BPCells"] with 8 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ active_transforms: logi [1:3, 1:2] FALSE TRUE TRUE FALSE FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:3] "row" "col" "global"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:2] "scale" "shift"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix           :Formal class 'ConvertMatrixType' [package "BPCells"] with 5 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'MatrixSubset' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'MatrixDir' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dir        : chr "/tmp/RtmpNY4Wop/matrix_dir"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ compressed : logi TRUE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ buffer_size: int 8192
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ type       : chr "uint32_t"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim        : int [1:2] 230 240
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose  : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames   :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_selection: int [1:230] 1 2 3 4 5 6 7 8 9 10 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_selection: int [1:120] 1 2 3 4 5 6 7 8 9 10 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ zero_dims    : logi [1:2] FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ cells   : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ type     : chr "double"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_params       : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_params       : num [1:2, 1:120] 0.0143 0 0.0118 0 0.0115 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ global_params    : num [1:2] 10000 0
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim              : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose        : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames         :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_params   : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_params   : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ global_params: num(0) 
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ cells   : chr [1:120] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'TransformLog1p' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'TransformScaleShift' [package "BPCells"] with 8 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ active_transforms: logi [1:3, 1:2] FALSE TRUE TRUE FALSE FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:3] "row" "col" "global"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:2] "scale" "shift"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix           :Formal class 'ConvertMatrixType' [package "BPCells"] with 5 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'RenameDims' [package "BPCells"] with 4 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix   :Formal class 'MatrixSubset' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ matrix       :Formal class 'MatrixDir' [package "BPCells"] with 7 slots
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dir        : chr "/tmp/RtmpNY4Wop/matrix_dir"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ compressed : logi TRUE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ buffer_size: int 8192
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ type       : chr "uint32_t"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim        : int [1:2] 230 240
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose  : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames   :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_selection: int [1:230] 1 2 3 4 5 6 7 8 9 10 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_selection: int [1:120] 121 122 123 124 125 126 127 128 129 130 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ zero_dims    : logi [1:2] FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ cells   : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ type     : chr "double"
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_params       : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_params       : num [1:2, 1:120] 0.00641 0 0.00719 0 0.00926 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ global_params    : num [1:2] 10000 0
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim              : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose        : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames         :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_params   : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_params   : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ global_params: num(0) 
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim      : int [1:2] 230 120
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ cells   : chr [1:120] "123defTAGGGACTGAACTC.1" "123defGCTCCATGAGAAGT.1" "123defTACAATGATGCTAG.1" "123defCTTCATGACCGAAT.1" ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ threads    : int 0
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim        : int [1:2] 230 240
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose  : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames   :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_selection: int [1:230] 203 177 28 174 190 197 133 194 200 202 ...
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_selection: int(0) 
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ zero_dims    : logi [1:2] FALSE FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 240
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. ..@ row_params   : num [1, 1:230] 20.3 21.9 13 26.5 18.1 ...
        .. .. .. .. .. .. .. .. .. .. .. .. ..@ col_params   : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. .. .. ..@ global_params: num(0) 
        .. .. .. .. .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 240
        .. .. .. .. .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. ..@ row_params       : num [1:2, 1:230] 0.519 -0.527 0.473 -0.372 0.784 ...
        .. .. .. .. .. .. .. .. .. .. ..@ col_params       : num[0 , 0 ] 
        .. .. .. .. .. .. .. .. .. .. ..@ global_params    : num(0) 
        .. .. .. .. .. .. .. .. .. .. ..@ dim              : int [1:2] 230 240
        .. .. .. .. .. .. .. .. .. .. ..@ transpose        : logi FALSE
        .. .. .. .. .. .. .. .. .. .. ..@ dimnames         :List of 2
        .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. ..@ row_selection: int [1:230] 93 220 85 103 78 105 115 210 194 97 ...
        .. .. .. .. .. .. .. .. ..@ col_selection: int(0) 
        .. .. .. .. .. .. .. .. ..@ zero_dims    : logi [1:2] FALSE FALSE
        .. .. .. .. .. .. .. .. ..@ dim          : int [1:2] 230 240
        .. .. .. .. .. .. .. .. ..@ transpose    : logi FALSE
        .. .. .. .. .. .. .. .. ..@ dimnames     :List of 2
        .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..@ dim      : int [1:2] 230 240
        .. .. .. .. .. .. ..@ transpose: logi FALSE
        .. .. .. .. .. .. ..@ dimnames :List of 2
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. .. .. ..$ features: chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. .. .. ..$ cells   : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. ..@ cells     :Formal class 'LogMap' [package "SeuratObject"] with 1 slot
        .. .. .. .. .. ..@ .Data: logi [1:240, 1:3] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. .. ..$ : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. .. .. ..$ : chr [1:3] "counts" "data" "scale.data"
        .. .. .. .. .. ..$ dim     : int [1:2] 240 3
        .. .. .. .. .. ..$ dimnames:List of 2
        .. .. .. .. .. .. ..$ : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
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
        .. .. .. .. ..$ vf_vst_counts.123abc_mean                 : num [1:230] 0.508 0.675 0.858 11.025 0.367 ...
        .. .. .. .. ..$ vf_vst_counts.123abc_variance             : num [1:230] 1.29 1.41 5.16 516.86 1.01 ...
        .. .. .. .. ..$ vf_vst_counts.123abc_variance.expected    : num [1:230] 1.68 3.91 5.87 506.74 1.15 ...
        .. .. .. .. ..$ vf_vst_counts.123abc_variance.standardized: num [1:230] 0.772 0.362 0.88 1.02 0.874 ...
        .. .. .. .. ..$ vf_vst_counts.123abc_variable             : logi [1:230] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. .. ..$ vf_vst_counts.123abc_rank                 : int [1:230] 136 221 114 93 117 96 151 186 172 101 ...
        .. .. .. .. ..$ vf_vst_counts.123def_mean                 : num [1:230] 0.267 0.525 0.542 15.825 0.233 ...
        .. .. .. .. ..$ vf_vst_counts.123def_variance             : num [1:230] 0.718 1.125 3.477 916.347 0.718 ...
        .. .. .. .. ..$ vf_vst_counts.123def_variance.expected    : num [1:230] 0.531 2.432 2.632 971.807 0.433 ...
        .. .. .. .. ..$ vf_vst_counts.123def_variance.standardized: num [1:230] 1.353 0.463 1.321 0.943 1.659 ...
        .. .. .. .. ..$ vf_vst_counts.123def_variable             : logi [1:230] TRUE TRUE TRUE TRUE TRUE TRUE ...
        .. .. .. .. ..$ vf_vst_counts.123def_rank                 : int [1:230] 54 213 58 108 32 107 74 220 209 94 ...
        .. .. .. .. ..$ var.features                              : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. .. .. ..$ var.features.rank                         : int [1:230] 93 220 85 103 78 105 115 210 194 97 ...
        .. .. .. ..@ misc      : list()
        .. .. .. ..@ key       : chr "rna_"
        ..@ meta.data   :'data.frame':	240 obs. of  9 variables:
        .. ..$ orig.ident    : chr [1:240] "SeuratProject" "SeuratProject" "SeuratProject" "SeuratProject" ...
        .. ..$ nCount_RNA    : num [1:240] 70 85 87 127 173 70 64 72 52 100 ...
        .. ..$ nFeature_RNA  : num [1:240] 47 52 50 56 53 48 36 45 36 41 ...
        .. ..$ cells_id      : int [1:240] 0 1 2 3 4 5 6 7 8 9 ...
        .. ..$ samples       : chr [1:240] "123abc" "123abc" "123abc" "123abc" ...
        .. ..$ doublet_scores: num [1:240] 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 ...
        .. ..$ doublet_class : chr [1:240] "singlet" "singlet" "singlet" "singlet" ...
        .. ..$ percent.mt    : num [1:240] 5.44 5.77 7.56 6.07 6.13 ...
        .. ..$ emptyDrops_FDR: num [1:240] 0.009 0.009 0.009 0.009 0.009 0.009 0.009 0.009 0.009 0.009 ...
        ..@ active.assay: chr "RNA"
        ..@ active.ident: Factor w/ 1 level "SeuratProject": 1 1 1 1 1 1 1 1 1 1 ...
        .. ..- attr(*, "names")= chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        ..@ graphs      : list()
        ..@ neighbors   : list()
        ..@ reductions  :List of 2
        .. ..$ pca            :Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
        .. .. .. ..@ cell.embeddings           : num [1:240, 1:50] 3.13 3.58 2.41 3.45 2.79 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. ..$ : chr [1:50] "PC_1" "PC_2" "PC_3" "PC_4" ...
        .. .. .. ..@ feature.loadings          : num [1:230, 1:50] 0.03005 0.06095 0.00738 0.05711 0.06203 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:230] "VDAC3" "PF4" "IGLL5" "PPBP" ...
        .. .. .. .. .. ..$ : chr [1:50] "PC_1" "PC_2" "PC_3" "PC_4" ...
        .. .. .. ..@ feature.loadings.projected: num[0 , 0 ] 
        .. .. .. ..@ assay.used                : chr "RNA"
        .. .. .. ..@ global                    : logi FALSE
        .. .. .. ..@ stdev                     : num [1:50] 5.75 5.21 4.32 3.62 2.77 ...
        .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
        .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ] 
        .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ] 
        .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ] 
        .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ] 
        .. .. .. ..@ misc                      :List of 1
        .. .. .. .. ..$ total.variance: num 230
        .. .. .. ..@ key                       : chr "PC_"
        .. ..$ integrated.rpca:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
        .. .. .. ..@ cell.embeddings           : num [1:240, 1:50] 3.13 3.58 2.41 3.45 2.79 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:240] "123abcATGCCAGAACGACT" "123abcCATGGCCTGTGCAT" "123abcGAACCTGATGAACC" "123abcTGACTGGATTCTCA" ...
        .. .. .. .. .. ..$ : chr [1:50] "integratedrpca_1" "integratedrpca_2" "integratedrpca_3" "integratedrpca_4" ...
        .. .. .. ..@ feature.loadings          : num [1:230, 1:50] 0.03005 0.06095 0.00738 0.05711 0.06203 ...
        .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. ..$ : chr [1:230] "VDAC3" "PF4" "IGLL5" "PPBP" ...
        .. .. .. .. .. ..$ : chr [1:50] "integratedrpca_1" "integratedrpca_2" "integratedrpca_3" "integratedrpca_4" ...
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
        .. .. .. ..@ key                       : chr "integratedrpca_"
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
        .. .. ..$ mean                 : num [1:230] 0.508 0.675 0.858 11.025 0.367 ...
        .. .. ..$ variance             : num [1:230] 1.29 1.41 5.16 516.86 1.01 ...
        .. .. ..$ variance.expected    : num [1:230] 1.68 3.91 5.87 506.74 1.15 ...
        .. .. ..$ variance.standardized: num [1:230] 0.772 0.362 0.88 1.02 0.874 ...
        .. .. ..$ SYMBOL               : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. .. ..$ ENSEMBL              : chr [1:230] "MS4A1" "CD79B" "CD79A" "HLA-DRA" ...
        .. ..$ numPCs          : num 50
        .. ..$ active.reduction: chr "integrated.rpca"
        ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        .. ..$ : int [1:3] 5 3 0
        ..@ commands    :List of 4
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
        .. ..$ ScaleData.RNA           :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "ScaleData.RNA"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : chr "RNA"
        .. .. .. ..@ call.string: chr "ScaleData(scdata, verbose = FALSE)"
        .. .. .. ..@ params     :List of 9
        .. .. .. .. ..$ assay             : chr "RNA"
        .. .. .. .. ..$ model.use         : chr "linear"
        .. .. .. .. ..$ use.umi           : logi FALSE
        .. .. .. .. ..$ do.scale          : logi TRUE
        .. .. .. .. ..$ do.center         : logi TRUE
        .. .. .. .. ..$ scale.max         : num 10
        .. .. .. .. ..$ block.size        : num 1000
        .. .. .. .. ..$ min.cells.to.block: num 3000
        .. .. .. .. ..$ verbose           : logi FALSE
        .. ..$ RunPCA.RNA              :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
        .. .. .. ..@ name       : chr "RunPCA.RNA"
        .. .. .. ..@ time.stamp : POSIXct[1:1], format: "1991-12-19 05:23:00"
        .. .. .. ..@ assay.used : chr "RNA"
        .. .. .. ..@ call.string: chr [1:2] "RunPCA(scdata, npcs = npcs, reduction.name = reduction_name, " "    verbose = FALSE)"
        .. .. .. ..@ params     :List of 10
        .. .. .. .. ..$ assay          : chr "RNA"
        .. .. .. .. ..$ npcs           : num 50
        .. .. .. .. ..$ rev.pca        : logi FALSE
        .. .. .. .. ..$ weight.by.var  : logi TRUE
        .. .. .. .. ..$ verbose        : logi FALSE
        .. .. .. .. ..$ ndims.print    : int [1:5] 1 2 3 4 5
        .. .. .. .. ..$ nfeatures.print: num 30
        .. .. .. .. ..$ reduction.name : chr "pca"
        .. .. .. .. ..$ reduction.key  : chr "PC_"
        .. .. .. .. ..$ seed.use       : num 42
        ..@ tools       : list()

