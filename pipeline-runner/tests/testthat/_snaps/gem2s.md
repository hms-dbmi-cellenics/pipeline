# gem2s is reproducible

    Code
      task_name
    Output
      [1] "downloadGem"
    Code
      rlang::hash(res)
    Output
      [1] "8c13d45d1c7104ff7e5f9d04485622fc"
    Code
      str(res)
    Output
      List of 2
       $ data  : list()
       $ output:List of 1
        ..$ config:List of 5
        .. ..$ name         : chr "mock_experiment"
        .. ..$ samples      :List of 2
        .. .. ..$ : chr "mock_sample_1_id"
        .. .. ..$ : chr "mock_sample_2_id"
        .. ..$ organism     : Named list()
        .. ..$ input        :List of 1
        .. .. ..$ type: chr "10x"
        .. ..$ sampleOptions: NULL

---

    Code
      task_name
    Output
      [1] "preproc"
    Code
      rlang::hash(res)
    Output
      [1] "234cece1f989e63edff7f55ab22bcff4"
    Code
      str(res)
    Output
      List of 2
       $ data  : list()
       $ output:List of 3
        ..$ config     :List of 5
        .. ..$ name         : chr "mock_experiment"
        .. ..$ samples      :List of 2
        .. .. ..$ : chr "mock_sample_1_id"
        .. .. ..$ : chr "mock_sample_2_id"
        .. ..$ organism     : Named list()
        .. ..$ input        :List of 1
        .. .. ..$ type: chr "10x"
        .. ..$ sampleOptions: NULL
        ..$ counts_list:List of 2
        .. ..$ mock_sample_1_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:36762] 66 454 1662 1694 1964 70 194 290 307 435 ...
        .. .. .. ..@ p       : int [1:501] 0 5 41 179 201 204 270 421 503 699 ...
        .. .. .. ..@ Dim     : int [1:2] 2000 500
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" ...
        .. .. .. ..@ x       : num [1:36762] 1 1 2 10 1 1 2 1 1 1 ...
        .. .. .. ..@ factors : list()
        .. ..$ mock_sample_2_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:33141] 23 177 207 341 451 470 663 834 838 952 ...
        .. .. .. ..@ p       : int [1:502] 0 24 30 110 169 264 315 349 457 496 ...
        .. .. .. ..@ Dim     : int [1:2] 2000 501
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" ...
        .. .. .. ..@ x       : num [1:33141] 1 1 3 1 2 5 1 1 1 1 ...
        .. .. .. ..@ factors : list()
        ..$ annot      :'data.frame':	2000 obs. of  3 variables:
        .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...

---

    Code
      task_name
    Output
      [1] "emptyDrops"
    Code
      rlang::hash(res)
    Output
      [1] "b0d3cba10eae73eb3d113cfca3b29384"
    Code
      str(res)
    Output
      List of 2
       $ data  : list()
       $ output:List of 4
        ..$ config     :List of 5
        .. ..$ name         : chr "mock_experiment"
        .. ..$ samples      :List of 2
        .. .. ..$ : chr "mock_sample_1_id"
        .. .. ..$ : chr "mock_sample_2_id"
        .. ..$ organism     : Named list()
        .. ..$ input        :List of 1
        .. .. ..$ type: chr "10x"
        .. ..$ sampleOptions: NULL
        ..$ counts_list:List of 2
        .. ..$ mock_sample_1_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:36762] 66 454 1662 1694 1964 70 194 290 307 435 ...
        .. .. .. ..@ p       : int [1:501] 0 5 41 179 201 204 270 421 503 699 ...
        .. .. .. ..@ Dim     : int [1:2] 2000 500
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" ...
        .. .. .. ..@ x       : num [1:36762] 1 1 2 10 1 1 2 1 1 1 ...
        .. .. .. ..@ factors : list()
        .. ..$ mock_sample_2_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:33141] 23 177 207 341 451 470 663 834 838 952 ...
        .. .. .. ..@ p       : int [1:502] 0 24 30 110 169 264 315 349 457 496 ...
        .. .. .. ..@ Dim     : int [1:2] 2000 501
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" ...
        .. .. .. ..@ x       : num [1:33141] 1 1 3 1 2 5 1 1 1 1 ...
        .. .. .. ..@ factors : list()
        ..$ annot      :'data.frame':	2000 obs. of  3 variables:
        .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        ..$ edrops     :List of 2
        .. ..$ mock_sample_1_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" ...
        .. .. .. ..@ nrows          : int 500
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:963, 1] 6.16e-05 1.75e-04 6.16e-05 6.16e-05 1.75e-04 ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:963] "CERNA1" "GPR174" "DCLK2" "AC001226.2" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 68.2
        .. .. .. .. ..$ retain : num 515
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:500] 15 208 422 49 6 520 578 329 637 442 ...
        .. .. .. .. ..$ LogProb: num [1:500] NA -173 -587 NA NA ...
        .. .. .. .. ..$ PValue : num [1:500] NA 0.9235 0.0001 NA NA ...
        .. .. .. .. ..$ Limited: logi [1:500] NA FALSE TRUE NA NA FALSE ...
        .. .. .. .. ..$ FDR    : num [1:500] NA 0.988069 0.000141 NA NA ...
        .. ..$ mock_sample_2_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" ...
        .. .. .. ..@ nrows          : int 501
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:893, 1] 8.09e-05 8.09e-05 8.09e-05 8.09e-05 2.60e-04 ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:893] "CERNA1" "GPR174" "AC001226.2" "AC080129.2" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 85.4
        .. .. .. .. ..$ retain : num 195
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:501] 46 21 289 249 263 193 93 472 101 112 ...
        .. .. .. .. ..$ LogProb: num [1:501] NA NA -367 -268 -404 ...
        .. .. .. .. ..$ PValue : num [1:501] NA NA 1e-04 9e-04 1e-04 ...
        .. .. .. .. ..$ Limited: logi [1:501] NA NA TRUE FALSE TRUE TRUE ...
        .. .. .. .. ..$ FDR    : num [1:501] NA NA 0 0 0 ...

---

    Code
      task_name
    Output
      [1] "doubletScores"
    Code
      rlang::hash(res)
    Output
      [1] "597d5c11ceb195b121e7b768aafccb4c"
    Code
      str(res)
    Output
      List of 2
       $ data  : list()
       $ output:List of 5
        ..$ config        :List of 5
        .. ..$ name         : chr "mock_experiment"
        .. ..$ samples      :List of 2
        .. .. ..$ : chr "mock_sample_1_id"
        .. .. ..$ : chr "mock_sample_2_id"
        .. ..$ organism     : Named list()
        .. ..$ input        :List of 1
        .. .. ..$ type: chr "10x"
        .. ..$ sampleOptions: NULL
        ..$ counts_list   :List of 2
        .. ..$ mock_sample_1_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:36762] 66 454 1662 1694 1964 70 194 290 307 435 ...
        .. .. .. ..@ p       : int [1:501] 0 5 41 179 201 204 270 421 503 699 ...
        .. .. .. ..@ Dim     : int [1:2] 2000 500
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" ...
        .. .. .. ..@ x       : num [1:36762] 1 1 2 10 1 1 2 1 1 1 ...
        .. .. .. ..@ factors : list()
        .. ..$ mock_sample_2_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:33141] 23 177 207 341 451 470 663 834 838 952 ...
        .. .. .. ..@ p       : int [1:502] 0 24 30 110 169 264 315 349 457 496 ...
        .. .. .. ..@ Dim     : int [1:2] 2000 501
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" ...
        .. .. .. ..@ x       : num [1:33141] 1 1 3 1 2 5 1 1 1 1 ...
        .. .. .. ..@ factors : list()
        ..$ annot         :'data.frame':	2000 obs. of  3 variables:
        .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        ..$ edrops        :List of 2
        .. ..$ mock_sample_1_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" ...
        .. .. .. ..@ nrows          : int 500
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:963, 1] 6.16e-05 1.75e-04 6.16e-05 6.16e-05 1.75e-04 ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:963] "CERNA1" "GPR174" "DCLK2" "AC001226.2" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 68.2
        .. .. .. .. ..$ retain : num 515
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:500] 15 208 422 49 6 520 578 329 637 442 ...
        .. .. .. .. ..$ LogProb: num [1:500] NA -173 -587 NA NA ...
        .. .. .. .. ..$ PValue : num [1:500] NA 0.9235 0.0001 NA NA ...
        .. .. .. .. ..$ Limited: logi [1:500] NA FALSE TRUE NA NA FALSE ...
        .. .. .. .. ..$ FDR    : num [1:500] NA 0.988069 0.000141 NA NA ...
        .. ..$ mock_sample_2_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" ...
        .. .. .. ..@ nrows          : int 501
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:893, 1] 8.09e-05 8.09e-05 8.09e-05 8.09e-05 2.60e-04 ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:893] "CERNA1" "GPR174" "AC001226.2" "AC080129.2" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 85.4
        .. .. .. .. ..$ retain : num 195
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:501] 46 21 289 249 263 193 93 472 101 112 ...
        .. .. .. .. ..$ LogProb: num [1:501] NA NA -367 -268 -404 ...
        .. .. .. .. ..$ PValue : num [1:501] NA NA 1e-04 9e-04 1e-04 ...
        .. .. .. .. ..$ Limited: logi [1:501] NA NA TRUE FALSE TRUE TRUE ...
        .. .. .. .. ..$ FDR    : num [1:501] NA NA 0 0 0 ...
        ..$ doublet_scores:List of 2
        .. ..$ mock_sample_1_id:'data.frame':	216 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:216] "GGTCACGAGTGAGCCA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 1 1 1 1 1 1 1 1 ...
        .. .. ..$ doublet_scores: num [1:216] 0.337 0.379 0.28 0.741 0.189 ...
        .. ..$ mock_sample_2_id:'data.frame':	242 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:242] "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "CAGCGTGGTACGATCT-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 2 2 1 2 1 2 2 1 ...
        .. .. ..$ doublet_scores: num [1:242] 0.13 0.295 0.923 0.984 0.406 ...

---

    Code
      task_name
    Output
      [1] "createSeurat"
    Code
      rlang::hash(res)
    Output
      [1] "2c31159cdfa0c963cf61adb1008b0314"
    Code
      str(res)
    Output
      List of 2
       $ data  : list()
       $ output:List of 7
        ..$ config            :List of 5
        .. ..$ name         : chr "mock_experiment"
        .. ..$ samples      :List of 2
        .. .. ..$ : chr "mock_sample_1_id"
        .. .. ..$ : chr "mock_sample_2_id"
        .. ..$ organism     : Named list()
        .. ..$ input        :List of 1
        .. .. ..$ type: chr "10x"
        .. ..$ sampleOptions: NULL
        ..$ counts_list       :List of 2
        .. ..$ mock_sample_1_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:36762] 66 454 1662 1694 1964 70 194 290 307 435 ...
        .. .. .. ..@ p       : int [1:501] 0 5 41 179 201 204 270 421 503 699 ...
        .. .. .. ..@ Dim     : int [1:2] 2000 500
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" ...
        .. .. .. ..@ x       : num [1:36762] 1 1 2 10 1 1 2 1 1 1 ...
        .. .. .. ..@ factors : list()
        .. ..$ mock_sample_2_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:33141] 23 177 207 341 451 470 663 834 838 952 ...
        .. .. .. ..@ p       : int [1:502] 0 24 30 110 169 264 315 349 457 496 ...
        .. .. .. ..@ Dim     : int [1:2] 2000 501
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" ...
        .. .. .. ..@ x       : num [1:33141] 1 1 3 1 2 5 1 1 1 1 ...
        .. .. .. ..@ factors : list()
        ..$ annot             :'data.frame':	2000 obs. of  3 variables:
        .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        ..$ edrops            :List of 2
        .. ..$ mock_sample_1_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" ...
        .. .. .. ..@ nrows          : int 500
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:963, 1] 6.16e-05 1.75e-04 6.16e-05 6.16e-05 1.75e-04 ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:963] "CERNA1" "GPR174" "DCLK2" "AC001226.2" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 68.2
        .. .. .. .. ..$ retain : num 515
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:500] 15 208 422 49 6 520 578 329 637 442 ...
        .. .. .. .. ..$ LogProb: num [1:500] NA -173 -587 NA NA ...
        .. .. .. .. ..$ PValue : num [1:500] NA 0.9235 0.0001 NA NA ...
        .. .. .. .. ..$ Limited: logi [1:500] NA FALSE TRUE NA NA FALSE ...
        .. .. .. .. ..$ FDR    : num [1:500] NA 0.988069 0.000141 NA NA ...
        .. ..$ mock_sample_2_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" ...
        .. .. .. ..@ nrows          : int 501
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:893, 1] 8.09e-05 8.09e-05 8.09e-05 8.09e-05 2.60e-04 ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:893] "CERNA1" "GPR174" "AC001226.2" "AC080129.2" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 85.4
        .. .. .. .. ..$ retain : num 195
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:501] 46 21 289 249 263 193 93 472 101 112 ...
        .. .. .. .. ..$ LogProb: num [1:501] NA NA -367 -268 -404 ...
        .. .. .. .. ..$ PValue : num [1:501] NA NA 1e-04 9e-04 1e-04 ...
        .. .. .. .. ..$ Limited: logi [1:501] NA NA TRUE FALSE TRUE TRUE ...
        .. .. .. .. ..$ FDR    : num [1:501] NA NA 0 0 0 ...
        ..$ doublet_scores    :List of 2
        .. ..$ mock_sample_1_id:'data.frame':	216 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:216] "GGTCACGAGTGAGCCA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 1 1 1 1 1 1 1 1 ...
        .. .. ..$ doublet_scores: num [1:216] 0.337 0.379 0.28 0.741 0.189 ...
        .. ..$ mock_sample_2_id:'data.frame':	242 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:242] "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "CAGCGTGGTACGATCT-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 2 2 1 2 1 2 2 1 ...
        .. .. ..$ doublet_scores: num [1:242] 0.13 0.295 0.923 0.984 0.406 ...
        ..$ scdata_list       :List of 2
        .. ..$ mock_sample_1_id:Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        .. .. .. ..@ assays      :List of 1
        .. .. .. .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
        .. .. .. .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:36426] 70 194 290 307 435 467 506 512 573 618 ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:437] 0 36 174 196 262 413 495 691 821 877 ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 436
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:36426] 1 2 1 1 1 1 3 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:36426] 70 194 290 307 435 467 506 512 573 618 ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:437] 0 36 174 196 262 413 495 691 821 877 ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 436
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:36426] 1 2 1 1 1 1 3 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ scale.data   : num[0 , 0 ] 
        .. .. .. .. .. .. ..@ key          : chr "rna_"
        .. .. .. .. .. .. ..@ assay.orig   : NULL
        .. .. .. .. .. .. ..@ var.features : logi(0) 
        .. .. .. .. .. .. ..@ meta.features:'data.frame':	2000 obs. of  0 variables
        .. .. .. .. .. .. ..@ misc         : list()
        .. .. .. ..@ meta.data   :'data.frame':	436 obs. of  13 variables:
        .. .. .. .. ..$ barcode           : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" ...
        .. .. .. .. ..$ orig.ident        : Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..$ nCount_RNA        : num [1:436] 208 422 49 520 578 329 637 442 229 521 ...
        .. .. .. .. ..$ nFeature_RNA      : int [1:436] 36 138 22 66 151 82 196 130 56 57 ...
        .. .. .. .. ..$ samples           : chr [1:436] "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" ...
        .. .. .. .. ..$ percent.mt        : num [1:436] 76.4 25.4 53.1 64.2 11.2 ...
        .. .. .. .. ..$ doublet_scores    : num [1:436] NA 0.337 NA 0.379 0.28 ...
        .. .. .. .. ..$ doublet_class     : Factor w/ 2 levels "singlet","doublet": NA 1 NA 1 1 1 1 1 NA 1 ...
        .. .. .. .. ..$ emptyDrops_Total  : int [1:436] 208 422 49 520 578 329 637 442 229 521 ...
        .. .. .. .. ..$ emptyDrops_LogProb: num [1:436] -173 -587 NA -333 -667 ...
        .. .. .. .. ..$ emptyDrops_PValue : num [1:436] 0.9235 0.0001 NA 0.1191 0.0001 ...
        .. .. .. .. ..$ emptyDrops_Limited: logi [1:436] FALSE TRUE NA FALSE TRUE TRUE ...
        .. .. .. .. ..$ emptyDrops_FDR    : num [1:436] 0.988069 0.000141 NA 0 0 ...
        .. .. .. ..@ active.assay: chr "RNA"
        .. .. .. ..@ active.ident: Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..- attr(*, "names")= chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" ...
        .. .. .. ..@ graphs      : list()
        .. .. .. ..@ neighbors   : list()
        .. .. .. ..@ reductions  : list()
        .. .. .. ..@ images      : list()
        .. .. .. ..@ project.name: chr "mock_experiment"
        .. .. .. ..@ misc        : list()
        .. .. .. ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        .. .. .. .. ..$ : int [1:3] 4 1 0
        .. .. .. ..@ commands    : list()
        .. .. .. ..@ tools       :List of 1
        .. .. .. .. ..$ flag_filtered: logi FALSE
        .. ..$ mock_sample_2_id:Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        .. .. .. ..@ assays      :List of 1
        .. .. .. .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
        .. .. .. .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:32918] 23 177 207 341 451 470 663 834 838 952 ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:461] 0 24 104 163 258 309 343 451 490 532 ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 460
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:32918] 1 1 3 1 2 5 1 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:32918] 23 177 207 341 451 470 663 834 838 952 ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:461] 0 24 104 163 258 309 343 451 490 532 ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 460
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:32918] 1 1 3 1 2 5 1 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ scale.data   : num[0 , 0 ] 
        .. .. .. .. .. .. ..@ key          : chr "rna_"
        .. .. .. .. .. .. ..@ assay.orig   : NULL
        .. .. .. .. .. .. ..@ var.features : logi(0) 
        .. .. .. .. .. .. ..@ meta.features:'data.frame':	2000 obs. of  0 variables
        .. .. .. .. .. .. ..@ misc         : list()
        .. .. .. ..@ meta.data   :'data.frame':	460 obs. of  13 variables:
        .. .. .. .. ..$ barcode           : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" ...
        .. .. .. .. ..$ orig.ident        : Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..$ nCount_RNA        : num [1:460] 46 289 249 263 193 93 472 101 112 41 ...
        .. .. .. .. ..$ nFeature_RNA      : int [1:460] 24 80 59 95 51 34 108 39 42 14 ...
        .. .. .. .. ..$ samples           : chr [1:460] "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" ...
        .. .. .. .. ..$ percent.mt        : num [1:460] 19.57 12.8 25.3 30.8 7.77 ...
        .. .. .. .. ..$ doublet_scores    : num [1:460] NA 0.13 0.295 0.923 NA ...
        .. .. .. .. ..$ doublet_class     : Factor w/ 2 levels "singlet","doublet": NA 1 1 2 NA NA 2 NA NA NA ...
        .. .. .. .. ..$ emptyDrops_Total  : int [1:460] 46 289 249 263 193 93 472 101 112 41 ...
        .. .. .. .. ..$ emptyDrops_LogProb: num [1:460] NA -367 -268 -404 -268 ...
        .. .. .. .. ..$ emptyDrops_PValue : num [1:460] NA 1e-04 9e-04 1e-04 1e-04 ...
        .. .. .. .. ..$ emptyDrops_Limited: logi [1:460] NA TRUE FALSE TRUE TRUE NA ...
        .. .. .. .. ..$ emptyDrops_FDR    : num [1:460] NA 0 0 0 0.000124 ...
        .. .. .. ..@ active.assay: chr "RNA"
        .. .. .. ..@ active.ident: Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..- attr(*, "names")= chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" ...
        .. .. .. ..@ graphs      : list()
        .. .. .. ..@ neighbors   : list()
        .. .. .. ..@ reductions  : list()
        .. .. .. ..@ images      : list()
        .. .. .. ..@ project.name: chr "mock_experiment"
        .. .. .. ..@ misc        : list()
        .. .. .. ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        .. .. .. .. ..$ : int [1:3] 4 1 0
        .. .. .. ..@ commands    : list()
        .. .. .. ..@ tools       :List of 1
        .. .. .. .. ..$ flag_filtered: logi FALSE
        ..$ disable_qc_filters: logi FALSE

---

    Code
      task_name
    Output
      [1] "prepareExperiment"
    Code
      rlang::hash(res)
    Output
      [1] "2e95b4be5dd6efe24f65415243a6496f"
    Code
      str(res)
    Output
      List of 2
       $ data  : list()
       $ output:List of 8
        ..$ config            :List of 5
        .. ..$ name         : chr "mock_experiment"
        .. ..$ samples      :List of 2
        .. .. ..$ : chr "mock_sample_1_id"
        .. .. ..$ : chr "mock_sample_2_id"
        .. ..$ organism     : Named list()
        .. ..$ input        :List of 1
        .. .. ..$ type: chr "10x"
        .. ..$ sampleOptions: NULL
        ..$ counts_list       :List of 2
        .. ..$ mock_sample_1_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:36762] 66 454 1662 1694 1964 70 194 290 307 435 ...
        .. .. .. ..@ p       : int [1:501] 0 5 41 179 201 204 270 421 503 699 ...
        .. .. .. ..@ Dim     : int [1:2] 2000 500
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" ...
        .. .. .. ..@ x       : num [1:36762] 1 1 2 10 1 1 2 1 1 1 ...
        .. .. .. ..@ factors : list()
        .. ..$ mock_sample_2_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:33141] 23 177 207 341 451 470 663 834 838 952 ...
        .. .. .. ..@ p       : int [1:502] 0 24 30 110 169 264 315 349 457 496 ...
        .. .. .. ..@ Dim     : int [1:2] 2000 501
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" ...
        .. .. .. ..@ x       : num [1:33141] 1 1 3 1 2 5 1 1 1 1 ...
        .. .. .. ..@ factors : list()
        ..$ annot             :'data.frame':	2000 obs. of  3 variables:
        .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        ..$ edrops            :List of 2
        .. ..$ mock_sample_1_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" ...
        .. .. .. ..@ nrows          : int 500
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:963, 1] 6.16e-05 1.75e-04 6.16e-05 6.16e-05 1.75e-04 ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:963] "CERNA1" "GPR174" "DCLK2" "AC001226.2" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 68.2
        .. .. .. .. ..$ retain : num 515
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:500] 15 208 422 49 6 520 578 329 637 442 ...
        .. .. .. .. ..$ LogProb: num [1:500] NA -173 -587 NA NA ...
        .. .. .. .. ..$ PValue : num [1:500] NA 0.9235 0.0001 NA NA ...
        .. .. .. .. ..$ Limited: logi [1:500] NA FALSE TRUE NA NA FALSE ...
        .. .. .. .. ..$ FDR    : num [1:500] NA 0.988069 0.000141 NA NA ...
        .. ..$ mock_sample_2_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" ...
        .. .. .. ..@ nrows          : int 501
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:893, 1] 8.09e-05 8.09e-05 8.09e-05 8.09e-05 2.60e-04 ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:893] "CERNA1" "GPR174" "AC001226.2" "AC080129.2" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 85.4
        .. .. .. .. ..$ retain : num 195
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:501] 46 21 289 249 263 193 93 472 101 112 ...
        .. .. .. .. ..$ LogProb: num [1:501] NA NA -367 -268 -404 ...
        .. .. .. .. ..$ PValue : num [1:501] NA NA 1e-04 9e-04 1e-04 ...
        .. .. .. .. ..$ Limited: logi [1:501] NA NA TRUE FALSE TRUE TRUE ...
        .. .. .. .. ..$ FDR    : num [1:501] NA NA 0 0 0 ...
        ..$ doublet_scores    :List of 2
        .. ..$ mock_sample_1_id:'data.frame':	216 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:216] "GGTCACGAGTGAGCCA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 1 1 1 1 1 1 1 1 ...
        .. .. ..$ doublet_scores: num [1:216] 0.337 0.379 0.28 0.741 0.189 ...
        .. ..$ mock_sample_2_id:'data.frame':	242 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:242] "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "CAGCGTGGTACGATCT-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 2 2 1 2 1 2 2 1 ...
        .. .. ..$ doublet_scores: num [1:242] 0.13 0.295 0.923 0.984 0.406 ...
        ..$ scdata_list       :List of 2
        .. ..$ mock_sample_2_id:Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        .. .. .. ..@ assays      :List of 1
        .. .. .. .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
        .. .. .. .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:32918] 23 177 207 341 451 470 663 834 838 952 ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:461] 0 24 104 163 258 309 343 451 490 532 ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 460
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:32918] 1 1 3 1 2 5 1 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:32918] 23 177 207 341 451 470 663 834 838 952 ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:461] 0 24 104 163 258 309 343 451 490 532 ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 460
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:32918] 1 1 3 1 2 5 1 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ scale.data   : num[0 , 0 ] 
        .. .. .. .. .. .. ..@ key          : chr "rna_"
        .. .. .. .. .. .. ..@ assay.orig   : NULL
        .. .. .. .. .. .. ..@ var.features : logi(0) 
        .. .. .. .. .. .. ..@ meta.features:'data.frame':	2000 obs. of  0 variables
        .. .. .. .. .. .. ..@ misc         : list()
        .. .. .. ..@ meta.data   :'data.frame':	460 obs. of  14 variables:
        .. .. .. .. ..$ barcode           : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" ...
        .. .. .. .. ..$ orig.ident        : Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..$ nCount_RNA        : num [1:460] 46 289 249 263 193 93 472 101 112 41 ...
        .. .. .. .. ..$ nFeature_RNA      : int [1:460] 24 80 59 95 51 34 108 39 42 14 ...
        .. .. .. .. ..$ samples           : chr [1:460] "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" ...
        .. .. .. .. ..$ percent.mt        : num [1:460] 19.57 12.8 25.3 30.8 7.77 ...
        .. .. .. .. ..$ doublet_scores    : num [1:460] NA 0.13 0.295 0.923 NA ...
        .. .. .. .. ..$ doublet_class     : Factor w/ 2 levels "singlet","doublet": NA 1 1 2 NA NA 2 NA NA NA ...
        .. .. .. .. ..$ emptyDrops_Total  : int [1:460] 46 289 249 263 193 93 472 101 112 41 ...
        .. .. .. .. ..$ emptyDrops_LogProb: num [1:460] NA -367 -268 -404 -268 ...
        .. .. .. .. ..$ emptyDrops_PValue : num [1:460] NA 1e-04 9e-04 1e-04 1e-04 ...
        .. .. .. .. ..$ emptyDrops_Limited: logi [1:460] NA TRUE FALSE TRUE TRUE NA ...
        .. .. .. .. ..$ emptyDrops_FDR    : num [1:460] NA 0 0 0 0.000124 ...
        .. .. .. .. ..$ cells_id          : int [1:460] 560 320 152 73 227 145 633 48 127 302 ...
        .. .. .. ..@ active.assay: chr "RNA"
        .. .. .. ..@ active.ident: Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..- attr(*, "names")= chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" ...
        .. .. .. ..@ graphs      : list()
        .. .. .. ..@ neighbors   : list()
        .. .. .. ..@ reductions  : list()
        .. .. .. ..@ images      : list()
        .. .. .. ..@ project.name: chr "mock_experiment"
        .. .. .. ..@ misc        :List of 2
        .. .. .. .. ..$ gene_annotations:'data.frame':	2000 obs. of  3 variables:
        .. .. .. .. .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ experimentId    : chr "mock_experiment_id"
        .. .. .. ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        .. .. .. .. ..$ : int [1:3] 4 1 0
        .. .. .. ..@ commands    : list()
        .. .. .. ..@ tools       :List of 1
        .. .. .. .. ..$ flag_filtered: logi FALSE
        .. ..$ mock_sample_1_id:Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        .. .. .. ..@ assays      :List of 1
        .. .. .. .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
        .. .. .. .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:36426] 70 194 290 307 435 467 506 512 573 618 ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:437] 0 36 174 196 262 413 495 691 821 877 ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 436
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:36426] 1 2 1 1 1 1 3 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:36426] 70 194 290 307 435 467 506 512 573 618 ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:437] 0 36 174 196 262 413 495 691 821 877 ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 436
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:36426] 1 2 1 1 1 1 3 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ scale.data   : num[0 , 0 ] 
        .. .. .. .. .. .. ..@ key          : chr "rna_"
        .. .. .. .. .. .. ..@ assay.orig   : NULL
        .. .. .. .. .. .. ..@ var.features : logi(0) 
        .. .. .. .. .. .. ..@ meta.features:'data.frame':	2000 obs. of  0 variables
        .. .. .. .. .. .. ..@ misc         : list()
        .. .. .. ..@ meta.data   :'data.frame':	436 obs. of  14 variables:
        .. .. .. .. ..$ barcode           : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" ...
        .. .. .. .. ..$ orig.ident        : Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..$ nCount_RNA        : num [1:436] 208 422 49 520 578 329 637 442 229 521 ...
        .. .. .. .. ..$ nFeature_RNA      : int [1:436] 36 138 22 66 151 82 196 130 56 57 ...
        .. .. .. .. ..$ samples           : chr [1:436] "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" ...
        .. .. .. .. ..$ percent.mt        : num [1:436] 76.4 25.4 53.1 64.2 11.2 ...
        .. .. .. .. ..$ doublet_scores    : num [1:436] NA 0.337 NA 0.379 0.28 ...
        .. .. .. .. ..$ doublet_class     : Factor w/ 2 levels "singlet","doublet": NA 1 NA 1 1 1 1 1 NA 1 ...
        .. .. .. .. ..$ emptyDrops_Total  : int [1:436] 208 422 49 520 578 329 637 442 229 521 ...
        .. .. .. .. ..$ emptyDrops_LogProb: num [1:436] -173 -587 NA -333 -667 ...
        .. .. .. .. ..$ emptyDrops_PValue : num [1:436] 0.9235 0.0001 NA 0.1191 0.0001 ...
        .. .. .. .. ..$ emptyDrops_Limited: logi [1:436] FALSE TRUE NA FALSE TRUE TRUE ...
        .. .. .. .. ..$ emptyDrops_FDR    : num [1:436] 0.988069 0.000141 NA 0 0 ...
        .. .. .. .. ..$ cells_id          : int [1:436] 216 894 758 202 851 518 475 14 761 399 ...
        .. .. .. ..@ active.assay: chr "RNA"
        .. .. .. ..@ active.ident: Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..- attr(*, "names")= chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" ...
        .. .. .. ..@ graphs      : list()
        .. .. .. ..@ neighbors   : list()
        .. .. .. ..@ reductions  : list()
        .. .. .. ..@ images      : list()
        .. .. .. ..@ project.name: chr "mock_experiment"
        .. .. .. ..@ misc        :List of 2
        .. .. .. .. ..$ gene_annotations:'data.frame':	2000 obs. of  3 variables:
        .. .. .. .. .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" ...
        .. .. .. .. ..$ experimentId    : chr "mock_experiment_id"
        .. .. .. ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        .. .. .. .. ..$ : int [1:3] 4 1 0
        .. .. .. ..@ commands    : list()
        .. .. .. ..@ tools       :List of 1
        .. .. .. .. ..$ flag_filtered: logi FALSE
        ..$ disable_qc_filters: logi FALSE
        ..$ qc_config         :List of 7
        .. ..$ cellSizeDistribution:List of 2
        .. .. ..$ mock_sample_2_id:List of 3
        .. .. .. ..$ enabled       : logi FALSE
        .. .. .. ..$ auto          : logi TRUE
        .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. ..$ minCellSize: num 17
        .. .. .. .. ..$ binStep    : num 200
        .. .. ..$ mock_sample_1_id:List of 3
        .. .. .. ..$ enabled       : logi FALSE
        .. .. .. ..$ auto          : logi TRUE
        .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. ..$ minCellSize: num 27
        .. .. .. .. ..$ binStep    : num 200
        .. ..$ mitochondrialContent:List of 2
        .. .. ..$ mock_sample_2_id:List of 3
        .. .. .. ..$ enabled       : logi TRUE
        .. .. .. ..$ auto          : logi TRUE
        .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. ..$ maxFraction: num 0.521
        .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. .. ..$ mock_sample_1_id:List of 3
        .. .. .. ..$ enabled       : logi TRUE
        .. .. .. ..$ auto          : logi TRUE
        .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. ..$ maxFraction: num 0.603
        .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. ..$ classifier          :List of 2
        .. .. ..$ mock_sample_2_id:List of 4
        .. .. .. ..$ enabled       : logi TRUE
        .. .. .. ..$ prefiltered   : logi FALSE
        .. .. .. ..$ auto          : logi TRUE
        .. .. .. ..$ filterSettings:List of 1
        .. .. .. .. ..$ FDR: num 0.01
        .. .. ..$ mock_sample_1_id:List of 4
        .. .. .. ..$ enabled       : logi TRUE
        .. .. .. ..$ prefiltered   : logi FALSE
        .. .. .. ..$ auto          : logi TRUE
        .. .. .. ..$ filterSettings:List of 1
        .. .. .. .. ..$ FDR: num 0.01
        .. ..$ numGenesVsNumUmis   :List of 2
        .. .. ..$ mock_sample_2_id:List of 3
        .. .. .. ..$ enabled       : logi TRUE
        .. .. .. ..$ auto          : logi TRUE
        .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. ..$ mock_sample_1_id:List of 3
        .. .. .. ..$ enabled       : logi TRUE
        .. .. .. ..$ auto          : logi TRUE
        .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. ..$ doubletScores       :List of 2
        .. .. ..$ mock_sample_2_id:List of 3
        .. .. .. ..$ enabled       : logi TRUE
        .. .. .. ..$ auto          : logi TRUE
        .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. ..$ probabilityThreshold: num 0.805
        .. .. .. .. ..$ binStep             : num 0.02
        .. .. ..$ mock_sample_1_id:List of 3
        .. .. .. ..$ enabled       : logi TRUE
        .. .. .. ..$ auto          : logi TRUE
        .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. ..$ probabilityThreshold: num 0.865
        .. .. .. .. ..$ binStep             : num 0.02
        .. ..$ dataIntegration     :List of 2
        .. .. ..$ dataIntegration        :List of 2
        .. .. .. ..$ method        : chr "harmony"
        .. .. .. ..$ methodSettings:List of 4
        .. .. .. .. ..$ seuratv4 :List of 2
        .. .. .. .. .. ..$ numGenes     : num 2000
        .. .. .. .. .. ..$ normalisation: chr "logNormalize"
        .. .. .. .. ..$ unisample:List of 2
        .. .. .. .. .. ..$ numGenes     : num 2000
        .. .. .. .. .. ..$ normalisation: chr "logNormalize"
        .. .. .. .. ..$ harmony  :List of 2
        .. .. .. .. .. ..$ numGenes     : num 2000
        .. .. .. .. .. ..$ normalisation: chr "logNormalize"
        .. .. .. .. ..$ fastmnn  :List of 2
        .. .. .. .. .. ..$ numGenes     : num 2000
        .. .. .. .. .. ..$ normalisation: chr "logNormalize"
        .. .. ..$ dimensionalityReduction:List of 3
        .. .. .. ..$ method               : chr "rpca"
        .. .. .. ..$ numPCs               : NULL
        .. .. .. ..$ excludeGeneCategories: list()
        .. ..$ configureEmbedding  :List of 2
        .. .. ..$ embeddingSettings :List of 2
        .. .. .. ..$ method        : chr "umap"
        .. .. .. ..$ methodSettings:List of 2
        .. .. .. .. ..$ umap:List of 2
        .. .. .. .. .. ..$ minimumDistance: num 0.3
        .. .. .. .. .. ..$ distanceMetric : chr "cosine"
        .. .. .. .. ..$ tsne:List of 2
        .. .. .. .. .. ..$ perplexity  : num 30
        .. .. .. .. .. ..$ learningRate: num 200
        .. .. ..$ clusteringSettings:List of 2
        .. .. .. ..$ method        : chr "louvain"
        .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. ..$ louvain:List of 1
        .. .. .. .. .. ..$ resolution: num 0.8

---

    Code
      task_name
    Output
      [1] "uploadToAWS"
    Code
      rlang::hash(res)
    Output
      [1] "0379510c6ca5436e14b53fad6d4ad7de"
    Code
      str(res)
    Output
      List of 2
       $ data  :List of 2
        ..$ item :List of 5
        .. ..$ apiVersion      : chr "2.0.0-data-ingest-seurat-rds-automated"
        .. ..$ experimentId    : chr "mock_experiment_id"
        .. ..$ experimentName  : chr "mock_experiment"
        .. ..$ meta            :List of 2
        .. .. ..$ organism: Named list()
        .. .. ..$ type    : chr "10x"
        .. ..$ processingConfig:List of 7
        .. .. ..$ cellSizeDistribution:List of 2
        .. .. .. ..$ mock_sample_2_id:List of 3
        .. .. .. .. ..$ enabled       : logi FALSE
        .. .. .. .. ..$ auto          : logi TRUE
        .. .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. .. ..$ minCellSize: num 17
        .. .. .. .. .. ..$ binStep    : num 200
        .. .. .. ..$ mock_sample_1_id:List of 3
        .. .. .. .. ..$ enabled       : logi FALSE
        .. .. .. .. ..$ auto          : logi TRUE
        .. .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. .. ..$ minCellSize: num 27
        .. .. .. .. .. ..$ binStep    : num 200
        .. .. ..$ mitochondrialContent:List of 2
        .. .. .. ..$ mock_sample_2_id:List of 3
        .. .. .. .. ..$ enabled       : logi TRUE
        .. .. .. .. ..$ auto          : logi TRUE
        .. .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. .. ..$ maxFraction: num 0.521
        .. .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. .. .. ..$ mock_sample_1_id:List of 3
        .. .. .. .. ..$ enabled       : logi TRUE
        .. .. .. .. ..$ auto          : logi TRUE
        .. .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. .. ..$ maxFraction: num 0.603
        .. .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. .. ..$ classifier          :List of 2
        .. .. .. ..$ mock_sample_2_id:List of 4
        .. .. .. .. ..$ enabled       : logi TRUE
        .. .. .. .. ..$ prefiltered   : logi FALSE
        .. .. .. .. ..$ auto          : logi TRUE
        .. .. .. .. ..$ filterSettings:List of 1
        .. .. .. .. .. ..$ FDR: num 0.01
        .. .. .. ..$ mock_sample_1_id:List of 4
        .. .. .. .. ..$ enabled       : logi TRUE
        .. .. .. .. ..$ prefiltered   : logi FALSE
        .. .. .. .. ..$ auto          : logi TRUE
        .. .. .. .. ..$ filterSettings:List of 1
        .. .. .. .. .. ..$ FDR: num 0.01
        .. .. ..$ numGenesVsNumUmis   :List of 2
        .. .. .. ..$ mock_sample_2_id:List of 3
        .. .. .. .. ..$ enabled       : logi TRUE
        .. .. .. .. ..$ auto          : logi TRUE
        .. .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. ..$ mock_sample_1_id:List of 3
        .. .. .. .. ..$ enabled       : logi TRUE
        .. .. .. .. ..$ auto          : logi TRUE
        .. .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. ..$ doubletScores       :List of 2
        .. .. .. ..$ mock_sample_2_id:List of 3
        .. .. .. .. ..$ enabled       : logi TRUE
        .. .. .. .. ..$ auto          : logi TRUE
        .. .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. .. ..$ probabilityThreshold: num 0.805
        .. .. .. .. .. ..$ binStep             : num 0.02
        .. .. .. ..$ mock_sample_1_id:List of 3
        .. .. .. .. ..$ enabled       : logi TRUE
        .. .. .. .. ..$ auto          : logi TRUE
        .. .. .. .. ..$ filterSettings:List of 2
        .. .. .. .. .. ..$ probabilityThreshold: num 0.865
        .. .. .. .. .. ..$ binStep             : num 0.02
        .. .. ..$ dataIntegration     :List of 2
        .. .. .. ..$ dataIntegration        :List of 2
        .. .. .. .. ..$ method        : chr "harmony"
        .. .. .. .. ..$ methodSettings:List of 4
        .. .. .. .. .. ..$ seuratv4 :List of 2
        .. .. .. .. .. .. ..$ numGenes     : num 2000
        .. .. .. .. .. .. ..$ normalisation: chr "logNormalize"
        .. .. .. .. .. ..$ unisample:List of 2
        .. .. .. .. .. .. ..$ numGenes     : num 2000
        .. .. .. .. .. .. ..$ normalisation: chr "logNormalize"
        .. .. .. .. .. ..$ harmony  :List of 2
        .. .. .. .. .. .. ..$ numGenes     : num 2000
        .. .. .. .. .. .. ..$ normalisation: chr "logNormalize"
        .. .. .. .. .. ..$ fastmnn  :List of 2
        .. .. .. .. .. .. ..$ numGenes     : num 2000
        .. .. .. .. .. .. ..$ normalisation: chr "logNormalize"
        .. .. .. ..$ dimensionalityReduction:List of 3
        .. .. .. .. ..$ method               : chr "rpca"
        .. .. .. .. ..$ numPCs               : NULL
        .. .. .. .. ..$ excludeGeneCategories: list()
        .. .. ..$ configureEmbedding  :List of 2
        .. .. .. ..$ embeddingSettings :List of 2
        .. .. .. .. ..$ method        : chr "umap"
        .. .. .. .. ..$ methodSettings:List of 2
        .. .. .. .. .. ..$ umap:List of 2
        .. .. .. .. .. .. ..$ minimumDistance: num 0.3
        .. .. .. .. .. .. ..$ distanceMetric : chr "cosine"
        .. .. .. .. .. ..$ tsne:List of 2
        .. .. .. .. .. .. ..$ perplexity  : num 30
        .. .. .. .. .. .. ..$ learningRate: num 200
        .. .. .. ..$ clusteringSettings:List of 2
        .. .. .. .. ..$ method        : chr "louvain"
        .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. ..$ louvain:List of 1
        .. .. .. .. .. .. ..$ resolution: num 0.8
        ..$ table: NULL
       $ output: list()

