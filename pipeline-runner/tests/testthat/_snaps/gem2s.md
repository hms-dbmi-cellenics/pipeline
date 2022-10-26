# some gem2s steps work

    Code
      str(res, vec.len = 20)
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
      str(res, vec.len = 20)
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
        .. .. .. ..@ i       : int [1:36762] 66 454 1662 1694 1964 70 194 290 307 435 467 506 512 573 618 641 666 674 762 812 819 934 947 1058 1091 1103 1367 | __truncated__ ...
        .. .. .. ..@ p       : int [1:501] 0 5 41 179 201 204 270 421 503 699 829 836 841 897 954 1001 1022 1052 1057 1252 1392 1441 1453 1461 1486 1565 157| __truncated__ ...
        .. .. .. ..@ Dim     : int [1:2] 2000 500
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. ..$ : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GCAACCGTCCTTTGAT-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "CGTCAAATCGGAAGGT-1" "GCACGGTGTGGCTGAA-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "TGGATGTCATCGCTCT-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" ...
        .. .. .. ..@ x       : num [1:36762] 1 1 2 10 1 1 2 1 1 1 1 3 1 1 1 1 1 3 1 1 1 1 1 1 1 1 3 1 1 1 2 159 3 1 1 1 1 1 3 1 3 1 1 2 1 4 1 1 1 26 ...
        .. .. .. ..@ factors : list()
        .. ..$ mock_sample_2_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:33141] 23 177 207 341 451 470 663 834 838 952 1016 1129 1316 1386 1434 1448 1479 1676 1694 1745 1871 1964 1969 1995 425 | __truncated__ ...
        .. .. .. ..@ p       : int [1:502] 0 24 30 110 169 264 315 349 457 496 538 552 565 640 679 798 873 967 1049 1112 1263 1353 1396 1477 1553 1647 1649 | __truncated__ ...
        .. .. .. ..@ Dim     : int [1:2] 2000 501
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. ..$ : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" ...
        .. .. .. ..@ x       : num [1:33141] 1 1 3 1 2 5 1 1 1 1 1 1 1 1 1 1 1 2 9 1 1 6 1 2 3 4 1 1 1 11 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. ..@ factors : list()
        ..$ annot      :'data.frame':	2000 obs. of  3 variables:
        .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...

---

    Code
      str(res, vec.len = 20)
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
        .. .. .. ..@ i       : int [1:36762] 66 454 1662 1694 1964 70 194 290 307 435 467 506 512 573 618 641 666 674 762 812 819 934 947 1058 1091 1103 1367 | __truncated__ ...
        .. .. .. ..@ p       : int [1:501] 0 5 41 179 201 204 270 421 503 699 829 836 841 897 954 1001 1022 1052 1057 1252 1392 1441 1453 1461 1486 1565 157| __truncated__ ...
        .. .. .. ..@ Dim     : int [1:2] 2000 500
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. ..$ : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GCAACCGTCCTTTGAT-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "CGTCAAATCGGAAGGT-1" "GCACGGTGTGGCTGAA-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "TGGATGTCATCGCTCT-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" ...
        .. .. .. ..@ x       : num [1:36762] 1 1 2 10 1 1 2 1 1 1 1 3 1 1 1 1 1 3 1 1 1 1 1 1 1 1 3 1 1 1 2 159 3 1 1 1 1 1 3 1 3 1 1 2 1 4 1 1 1 26 ...
        .. .. .. ..@ factors : list()
        .. ..$ mock_sample_2_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:33141] 23 177 207 341 451 470 663 834 838 952 1016 1129 1316 1386 1434 1448 1479 1676 1694 1745 1871 1964 1969 1995 425 | __truncated__ ...
        .. .. .. ..@ p       : int [1:502] 0 24 30 110 169 264 315 349 457 496 538 552 565 640 679 798 873 967 1049 1112 1263 1353 1396 1477 1553 1647 1649 | __truncated__ ...
        .. .. .. ..@ Dim     : int [1:2] 2000 501
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. ..$ : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" ...
        .. .. .. ..@ x       : num [1:33141] 1 1 3 1 2 5 1 1 1 1 1 1 1 1 1 1 1 2 9 1 1 6 1 2 3 4 1 1 1 11 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. ..@ factors : list()
        ..$ annot      :'data.frame':	2000 obs. of  3 variables:
        .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        ..$ edrops     :List of 2
        .. ..$ mock_sample_1_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GCAACCGTCCTTTGAT-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "CGTCAAATCGGAAGGT-1" "GCACGGTGTGGCTGAA-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "TGGATGTCATCGCTCT-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" ...
        .. .. .. ..@ nrows          : int 500
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:963, 1] 6.16e-05 1.75e-04 6.16e-05 6.16e-05 1.75e-04 1.04e-03 1.37e-03 3.90e-04 1.75e-04 3.90e-04 6.16e-05 1.75e-04 6.16e| __truncated__ ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:963] "CERNA1" "GPR174" "DCLK2" "AC001226.2" "PAX5" "SMG6" "EDF1" "ZNF605" "ANKS4B" "ACTR3B" "AC004067.1" "HSD17B4" "AC007881.3" "KCNAB3" "PSMA3" "TMSB15A" "LIPE" "PTX3" "DCBLD1" "FUNDC1" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 68.2
        .. .. .. .. ..$ retain : num 515
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:500] 15 208 422 49 6 520 578 329 637 442 13 436 229 521 119 59 42 206 594 441 80 26 10 58 223 120 52 262 149 133 164 3| __truncated__ ...
        .. .. .. .. ..$ LogProb: num [1:500] NA -173.1 -586.5 NA NA -332.5 -666.9 -375.6 -843.2 -548.1 NA -118.3 -268.8 -288.8 -204.3 NA NA -89.2 -842.1 -608.| __truncated__ ...
        .. .. .. .. ..$ PValue : num [1:500] NA 0.9235 0.0001 NA NA 0.1191 0.0001 0.0001 0.0001 0.0001 NA 1 0.0028 0.7043 0.0011 NA NA 1 0.0001 0.0001 NA NA NA NA 0.0001 ...
        .. .. .. .. ..$ Limited: logi [1:500] NA FALSE TRUE NA NA FALSE TRUE TRUE TRUE TRUE NA FALSE FALSE FALSE FALSE NA NA FALSE TRUE TRUE NA NA NA NA TRUE F| __truncated__ ...
        .. .. .. .. ..$ FDR    : num [1:500] NA 0.988069 0.000141 NA NA 0 0 0.000141 0 0.000141 NA 1 0.003363 0 0.001344 NA NA 1 0 0.000141 NA NA NA NA 0.000141 ...
        .. ..$ mock_sample_2_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" ...
        .. .. .. ..@ nrows          : int 501
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:893, 1] 8.09e-05 8.09e-05 8.09e-05 8.09e-05 2.60e-04 6.13e-04 2.09e-03 1.58e-04 8.09e-05 1.58e-04 1.58e-04 1.58e-04 2.60e| __truncated__ ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:893] "CERNA1" "GPR174" "AC001226.2" "AC080129.2" "PAX5" "SMG6" "EDF1" "ZNF605" "LINC02413" "ACTR3B" "HSD17B4" "AC007881.3" "KCNAB3" "PSMA3" "LIPE" "FOXJ1" "PTX3" "DCBLD1" "FUNDC1" "KCNE1" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 85.4
        .. .. .. .. ..$ retain : num 195
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:501] 46 21 289 249 263 193 93 472 101 112 41 38 297 81 372 272 174 406 107 529 238 93 232 165 300 64 411 373 533 155 3| __truncated__ ...
        .. .. .. .. ..$ LogProb: num [1:501] NA NA -367 -268 -404 -268 NA -490 -172 -198 NA NA -336 NA -500 -312 -381 -365 -274 -640 -372 NA -325 -311 -403 ...
        .. .. .. .. ..$ PValue : num [1:501] NA NA 0.0001 0.0009 0.0001 0.0001 NA 0.0001 0.0016 0.0002 NA NA 0.0001 NA 0.0001 0.0001 0.0001 0.0001 0.0001 0.00| __truncated__ ...
        .. .. .. .. ..$ Limited: logi [1:501] NA NA TRUE FALSE TRUE TRUE NA TRUE FALSE FALSE NA NA TRUE NA TRUE TRUE TRUE TRUE TRUE TRUE TRUE NA TRUE TRUE TRUE| __truncated__ ...
        .. .. .. .. ..$ FDR    : num [1:501] NA NA 0 0 0 0.000124 NA 0 0.001789 0.000239 NA NA 0 NA 0 0 0.000124 0 0.000124 0 0 NA 0 0.000124 0 ...

---

    Code
      str(res, vec.len = 20)
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
        .. .. .. ..@ i       : int [1:36762] 66 454 1662 1694 1964 70 194 290 307 435 467 506 512 573 618 641 666 674 762 812 819 934 947 1058 1091 1103 1367 | __truncated__ ...
        .. .. .. ..@ p       : int [1:501] 0 5 41 179 201 204 270 421 503 699 829 836 841 897 954 1001 1022 1052 1057 1252 1392 1441 1453 1461 1486 1565 157| __truncated__ ...
        .. .. .. ..@ Dim     : int [1:2] 2000 500
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. ..$ : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GCAACCGTCCTTTGAT-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "CGTCAAATCGGAAGGT-1" "GCACGGTGTGGCTGAA-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "TGGATGTCATCGCTCT-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" ...
        .. .. .. ..@ x       : num [1:36762] 1 1 2 10 1 1 2 1 1 1 1 3 1 1 1 1 1 3 1 1 1 1 1 1 1 1 3 1 1 1 2 159 3 1 1 1 1 1 3 1 3 1 1 2 1 4 1 1 1 26 ...
        .. .. .. ..@ factors : list()
        .. ..$ mock_sample_2_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:33141] 23 177 207 341 451 470 663 834 838 952 1016 1129 1316 1386 1434 1448 1479 1676 1694 1745 1871 1964 1969 1995 425 | __truncated__ ...
        .. .. .. ..@ p       : int [1:502] 0 24 30 110 169 264 315 349 457 496 538 552 565 640 679 798 873 967 1049 1112 1263 1353 1396 1477 1553 1647 1649 | __truncated__ ...
        .. .. .. ..@ Dim     : int [1:2] 2000 501
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. ..$ : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" ...
        .. .. .. ..@ x       : num [1:33141] 1 1 3 1 2 5 1 1 1 1 1 1 1 1 1 1 1 2 9 1 1 6 1 2 3 4 1 1 1 11 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. ..@ factors : list()
        ..$ annot         :'data.frame':	2000 obs. of  3 variables:
        .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        ..$ edrops        :List of 2
        .. ..$ mock_sample_1_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GCAACCGTCCTTTGAT-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "CGTCAAATCGGAAGGT-1" "GCACGGTGTGGCTGAA-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "TGGATGTCATCGCTCT-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" ...
        .. .. .. ..@ nrows          : int 500
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:963, 1] 6.16e-05 1.75e-04 6.16e-05 6.16e-05 1.75e-04 1.04e-03 1.37e-03 3.90e-04 1.75e-04 3.90e-04 6.16e-05 1.75e-04 6.16e| __truncated__ ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:963] "CERNA1" "GPR174" "DCLK2" "AC001226.2" "PAX5" "SMG6" "EDF1" "ZNF605" "ANKS4B" "ACTR3B" "AC004067.1" "HSD17B4" "AC007881.3" "KCNAB3" "PSMA3" "TMSB15A" "LIPE" "PTX3" "DCBLD1" "FUNDC1" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 68.2
        .. .. .. .. ..$ retain : num 515
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:500] 15 208 422 49 6 520 578 329 637 442 13 436 229 521 119 59 42 206 594 441 80 26 10 58 223 120 52 262 149 133 164 3| __truncated__ ...
        .. .. .. .. ..$ LogProb: num [1:500] NA -173.1 -586.5 NA NA -332.5 -666.9 -375.6 -843.2 -548.1 NA -118.3 -268.8 -288.8 -204.3 NA NA -89.2 -842.1 -608.| __truncated__ ...
        .. .. .. .. ..$ PValue : num [1:500] NA 0.9235 0.0001 NA NA 0.1191 0.0001 0.0001 0.0001 0.0001 NA 1 0.0028 0.7043 0.0011 NA NA 1 0.0001 0.0001 NA NA NA NA 0.0001 ...
        .. .. .. .. ..$ Limited: logi [1:500] NA FALSE TRUE NA NA FALSE TRUE TRUE TRUE TRUE NA FALSE FALSE FALSE FALSE NA NA FALSE TRUE TRUE NA NA NA NA TRUE F| __truncated__ ...
        .. .. .. .. ..$ FDR    : num [1:500] NA 0.988069 0.000141 NA NA 0 0 0.000141 0 0.000141 NA 1 0.003363 0 0.001344 NA NA 1 0 0.000141 NA NA NA NA 0.000141 ...
        .. ..$ mock_sample_2_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" ...
        .. .. .. ..@ nrows          : int 501
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:893, 1] 8.09e-05 8.09e-05 8.09e-05 8.09e-05 2.60e-04 6.13e-04 2.09e-03 1.58e-04 8.09e-05 1.58e-04 1.58e-04 1.58e-04 2.60e| __truncated__ ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:893] "CERNA1" "GPR174" "AC001226.2" "AC080129.2" "PAX5" "SMG6" "EDF1" "ZNF605" "LINC02413" "ACTR3B" "HSD17B4" "AC007881.3" "KCNAB3" "PSMA3" "LIPE" "FOXJ1" "PTX3" "DCBLD1" "FUNDC1" "KCNE1" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 85.4
        .. .. .. .. ..$ retain : num 195
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:501] 46 21 289 249 263 193 93 472 101 112 41 38 297 81 372 272 174 406 107 529 238 93 232 165 300 64 411 373 533 155 3| __truncated__ ...
        .. .. .. .. ..$ LogProb: num [1:501] NA NA -367 -268 -404 -268 NA -490 -172 -198 NA NA -336 NA -500 -312 -381 -365 -274 -640 -372 NA -325 -311 -403 ...
        .. .. .. .. ..$ PValue : num [1:501] NA NA 0.0001 0.0009 0.0001 0.0001 NA 0.0001 0.0016 0.0002 NA NA 0.0001 NA 0.0001 0.0001 0.0001 0.0001 0.0001 0.00| __truncated__ ...
        .. .. .. .. ..$ Limited: logi [1:501] NA NA TRUE FALSE TRUE TRUE NA TRUE FALSE FALSE NA NA TRUE NA TRUE TRUE TRUE TRUE TRUE TRUE TRUE NA TRUE TRUE TRUE| __truncated__ ...
        .. .. .. .. ..$ FDR    : num [1:501] NA NA 0 0 0 0.000124 NA 0 0.001789 0.000239 NA NA 0 NA 0 0 0.000124 0 0.000124 0 0 NA 0 0.000124 0 ...
        ..$ doublet_scores:List of 2
        .. ..$ mock_sample_1_id:'data.frame':	216 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:216] "GGTCACGAGTGAGCCA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "TTACCATCACTGTTCC-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" "TCATGGAAGAAGATCT-1" "CTCCAACGTGATTAGA-1" "AGTAACCAGTATTGCC-1" "TGTGCGGTCATTACTC-1" "AGATGCTTCGTTGCCT-1" "TGATGGTCAACCGATT-1" "AACCTTTCAGTGGTGA-1" "TATCCTAGTGATACTC-1" "TTTACGTAGCGATCGA-1" "GTGCTTCCAAATCGTC-1" "CCAATTTCAGGCTTGC-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 ...
        .. .. ..$ doublet_scores: num [1:216] 0.2368 0.3751 0.152 0.6659 0.1344 0.2483 0.4523 0.5696 0.2508 0.3544 0.4421 0.2415 0.1291 0.078 0.3669 0.0889 0.8| __truncated__ ...
        .. ..$ mock_sample_2_id:'data.frame':	242 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:242] "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "CAGCGTGGTACGATCT-1" "GCAGGCTCAAACTGCT-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "TGACAGTAGGCATGCA-1" "CAAGGGACAGAGCGTA-1" "TTTATGCAGAGGATCC-1" "AGCTCAAAGGCCCGTT-1" "TCGGGACGTAACATGA-1" "AATGACCCAGAGTCAG-1" "TCGTAGAGTGAGTCAG-1" "TCTTAGTAGGGATGTC-1" "ATCCTATCAAGACCTT-1" "CGATCGGGTGTGATGG-1" "CCTCTCCAGCGTGTTT-1" "TTCGATTGTTCCTACC-1" "CACGTTCTCGACCTAA-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 2 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. ..$ doublet_scores: num [1:242] 0.0458 0.5706 0.966 0.5052 0.3721 0.664 0.1622 0.564 0.9349 0.1361 0.0261 0.0479 0.1194 0.8091 0.9956 0.2044 0.99| __truncated__ ...

---

    Code
      str(res, vec.len = 20)
    Output
      List of 2
       $ data  : list()
       $ output:List of 6
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
        .. .. .. ..@ i       : int [1:36762] 66 454 1662 1694 1964 70 194 290 307 435 467 506 512 573 618 641 666 674 762 812 819 934 947 1058 1091 1103 1367 | __truncated__ ...
        .. .. .. ..@ p       : int [1:501] 0 5 41 179 201 204 270 421 503 699 829 836 841 897 954 1001 1022 1052 1057 1252 1392 1441 1453 1461 1486 1565 157| __truncated__ ...
        .. .. .. ..@ Dim     : int [1:2] 2000 500
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. ..$ : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GCAACCGTCCTTTGAT-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "CGTCAAATCGGAAGGT-1" "GCACGGTGTGGCTGAA-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "TGGATGTCATCGCTCT-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" ...
        .. .. .. ..@ x       : num [1:36762] 1 1 2 10 1 1 2 1 1 1 1 3 1 1 1 1 1 3 1 1 1 1 1 1 1 1 3 1 1 1 2 159 3 1 1 1 1 1 3 1 3 1 1 2 1 4 1 1 1 26 ...
        .. .. .. ..@ factors : list()
        .. ..$ mock_sample_2_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:33141] 23 177 207 341 451 470 663 834 838 952 1016 1129 1316 1386 1434 1448 1479 1676 1694 1745 1871 1964 1969 1995 425 | __truncated__ ...
        .. .. .. ..@ p       : int [1:502] 0 24 30 110 169 264 315 349 457 496 538 552 565 640 679 798 873 967 1049 1112 1263 1353 1396 1477 1553 1647 1649 | __truncated__ ...
        .. .. .. ..@ Dim     : int [1:2] 2000 501
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. ..$ : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" ...
        .. .. .. ..@ x       : num [1:33141] 1 1 3 1 2 5 1 1 1 1 1 1 1 1 1 1 1 2 9 1 1 6 1 2 3 4 1 1 1 11 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. ..@ factors : list()
        ..$ annot         :'data.frame':	2000 obs. of  3 variables:
        .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        ..$ edrops        :List of 2
        .. ..$ mock_sample_1_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GCAACCGTCCTTTGAT-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "CGTCAAATCGGAAGGT-1" "GCACGGTGTGGCTGAA-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "TGGATGTCATCGCTCT-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" ...
        .. .. .. ..@ nrows          : int 500
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:963, 1] 6.16e-05 1.75e-04 6.16e-05 6.16e-05 1.75e-04 1.04e-03 1.37e-03 3.90e-04 1.75e-04 3.90e-04 6.16e-05 1.75e-04 6.16e| __truncated__ ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:963] "CERNA1" "GPR174" "DCLK2" "AC001226.2" "PAX5" "SMG6" "EDF1" "ZNF605" "ANKS4B" "ACTR3B" "AC004067.1" "HSD17B4" "AC007881.3" "KCNAB3" "PSMA3" "TMSB15A" "LIPE" "PTX3" "DCBLD1" "FUNDC1" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 68.2
        .. .. .. .. ..$ retain : num 515
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:500] 15 208 422 49 6 520 578 329 637 442 13 436 229 521 119 59 42 206 594 441 80 26 10 58 223 120 52 262 149 133 164 3| __truncated__ ...
        .. .. .. .. ..$ LogProb: num [1:500] NA -173.1 -586.5 NA NA -332.5 -666.9 -375.6 -843.2 -548.1 NA -118.3 -268.8 -288.8 -204.3 NA NA -89.2 -842.1 -608.| __truncated__ ...
        .. .. .. .. ..$ PValue : num [1:500] NA 0.9235 0.0001 NA NA 0.1191 0.0001 0.0001 0.0001 0.0001 NA 1 0.0028 0.7043 0.0011 NA NA 1 0.0001 0.0001 NA NA NA NA 0.0001 ...
        .. .. .. .. ..$ Limited: logi [1:500] NA FALSE TRUE NA NA FALSE TRUE TRUE TRUE TRUE NA FALSE FALSE FALSE FALSE NA NA FALSE TRUE TRUE NA NA NA NA TRUE F| __truncated__ ...
        .. .. .. .. ..$ FDR    : num [1:500] NA 0.988069 0.000141 NA NA 0 0 0.000141 0 0.000141 NA 1 0.003363 0 0.001344 NA NA 1 0 0.000141 NA NA NA NA 0.000141 ...
        .. ..$ mock_sample_2_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" ...
        .. .. .. ..@ nrows          : int 501
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:893, 1] 8.09e-05 8.09e-05 8.09e-05 8.09e-05 2.60e-04 6.13e-04 2.09e-03 1.58e-04 8.09e-05 1.58e-04 1.58e-04 1.58e-04 2.60e| __truncated__ ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:893] "CERNA1" "GPR174" "AC001226.2" "AC080129.2" "PAX5" "SMG6" "EDF1" "ZNF605" "LINC02413" "ACTR3B" "HSD17B4" "AC007881.3" "KCNAB3" "PSMA3" "LIPE" "FOXJ1" "PTX3" "DCBLD1" "FUNDC1" "KCNE1" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 85.4
        .. .. .. .. ..$ retain : num 195
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:501] 46 21 289 249 263 193 93 472 101 112 41 38 297 81 372 272 174 406 107 529 238 93 232 165 300 64 411 373 533 155 3| __truncated__ ...
        .. .. .. .. ..$ LogProb: num [1:501] NA NA -367 -268 -404 -268 NA -490 -172 -198 NA NA -336 NA -500 -312 -381 -365 -274 -640 -372 NA -325 -311 -403 ...
        .. .. .. .. ..$ PValue : num [1:501] NA NA 0.0001 0.0009 0.0001 0.0001 NA 0.0001 0.0016 0.0002 NA NA 0.0001 NA 0.0001 0.0001 0.0001 0.0001 0.0001 0.00| __truncated__ ...
        .. .. .. .. ..$ Limited: logi [1:501] NA NA TRUE FALSE TRUE TRUE NA TRUE FALSE FALSE NA NA TRUE NA TRUE TRUE TRUE TRUE TRUE TRUE TRUE NA TRUE TRUE TRUE| __truncated__ ...
        .. .. .. .. ..$ FDR    : num [1:501] NA NA 0 0 0 0.000124 NA 0 0.001789 0.000239 NA NA 0 NA 0 0 0.000124 0 0.000124 0 0 NA 0 0.000124 0 ...
        ..$ doublet_scores:List of 2
        .. ..$ mock_sample_1_id:'data.frame':	216 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:216] "GGTCACGAGTGAGCCA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "TTACCATCACTGTTCC-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" "TCATGGAAGAAGATCT-1" "CTCCAACGTGATTAGA-1" "AGTAACCAGTATTGCC-1" "TGTGCGGTCATTACTC-1" "AGATGCTTCGTTGCCT-1" "TGATGGTCAACCGATT-1" "AACCTTTCAGTGGTGA-1" "TATCCTAGTGATACTC-1" "TTTACGTAGCGATCGA-1" "GTGCTTCCAAATCGTC-1" "CCAATTTCAGGCTTGC-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 ...
        .. .. ..$ doublet_scores: num [1:216] 0.2368 0.3751 0.152 0.6659 0.1344 0.2483 0.4523 0.5696 0.2508 0.3544 0.4421 0.2415 0.1291 0.078 0.3669 0.0889 0.8| __truncated__ ...
        .. ..$ mock_sample_2_id:'data.frame':	242 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:242] "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "CAGCGTGGTACGATCT-1" "GCAGGCTCAAACTGCT-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "TGACAGTAGGCATGCA-1" "CAAGGGACAGAGCGTA-1" "TTTATGCAGAGGATCC-1" "AGCTCAAAGGCCCGTT-1" "TCGGGACGTAACATGA-1" "AATGACCCAGAGTCAG-1" "TCGTAGAGTGAGTCAG-1" "TCTTAGTAGGGATGTC-1" "ATCCTATCAAGACCTT-1" "CGATCGGGTGTGATGG-1" "CCTCTCCAGCGTGTTT-1" "TTCGATTGTTCCTACC-1" "CACGTTCTCGACCTAA-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 2 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. ..$ doublet_scores: num [1:242] 0.0458 0.5706 0.966 0.5052 0.3721 0.664 0.1622 0.564 0.9349 0.1361 0.0261 0.0479 0.1194 0.8091 0.9956 0.2044 0.99| __truncated__ ...
        ..$ scdata_list   :List of 2
        .. ..$ mock_sample_1_id:Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        .. .. .. ..@ assays      :List of 1
        .. .. .. .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
        .. .. .. .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:36426] 70 194 290 307 435 467 506 512 573 618 641 666 674 762 812 819 934 947 1058 1091 1103 1367 1400 1571 1624 1665 16| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:437] 0 36 174 196 262 413 495 691 821 877 934 981 1002 1032 1227 1367 1416 1428 1453 1532 1559 1637 1709 1758 1819 190| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 436
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" "GGGCTCACAATGTCTG-1" "TATCTTGCAGGCGTTC-1" "TTGGGCGCATTGACCA-1" "TCATGGAAGAAGATCT-1" "GGGTTATAGAATTGTG-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:36426] 1 2 1 1 1 1 3 1 1 1 1 1 3 1 1 1 1 1 1 1 1 3 1 1 1 2 159 3 1 1 1 1 1 3 1 3 1 1 2 1 4 1 1 1 26 7 1 1 2 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:36426] 70 194 290 307 435 467 506 512 573 618 641 666 674 762 812 819 934 947 1058 1091 1103 1367 1400 1571 1624 1665 16| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:437] 0 36 174 196 262 413 495 691 821 877 934 981 1002 1032 1227 1367 1416 1428 1453 1532 1559 1637 1709 1758 1819 190| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 436
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" "GGGCTCACAATGTCTG-1" "TATCTTGCAGGCGTTC-1" "TTGGGCGCATTGACCA-1" "TCATGGAAGAAGATCT-1" "GGGTTATAGAATTGTG-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:36426] 1 2 1 1 1 1 3 1 1 1 1 1 3 1 1 1 1 1 1 1 1 3 1 1 1 2 159 3 1 1 1 1 1 3 1 3 1 1 2 1 4 1 1 1 26 7 1 1 2 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ scale.data   : num[0 , 0 ] 
        .. .. .. .. .. .. ..@ key          : chr "rna_"
        .. .. .. .. .. .. ..@ assay.orig   : NULL
        .. .. .. .. .. .. ..@ var.features : logi(0) 
        .. .. .. .. .. .. ..@ meta.features:'data.frame':	2000 obs. of  0 variables
        .. .. .. .. .. .. ..@ misc         : list()
        .. .. .. ..@ meta.data   :'data.frame':	436 obs. of  13 variables:
        .. .. .. .. ..$ barcode           : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" "GGGCTCACAATGTCTG-1" "TATCTTGCAGGCGTTC-1" "TTGGGCGCATTGACCA-1" "TCATGGAAGAAGATCT-1" "GGGTTATAGAATTGTG-1" ...
        .. .. .. .. ..$ orig.ident        : Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..$ nCount_RNA        : num [1:436] 208 422 49 520 578 329 637 442 229 521 119 59 42 594 441 80 26 58 223 52 262 149 133 164 314 743 68 345 593 83 35| __truncated__ ...
        .. .. .. .. ..$ nFeature_RNA      : int [1:436] 36 138 22 66 151 82 196 130 56 57 47 21 30 195 140 49 12 25 79 27 78 72 49 61 88 204 17 116 168 26 125 26 56 57 1| __truncated__ ...
        .. .. .. .. ..$ samples           : chr [1:436] "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" ...
        .. .. .. .. ..$ percent.mt        : num [1:436] 76.44 25.36 53.06 64.23 11.25 8.81 27.16 21.04 66.81 85.6 20.17 59.32 28.57 20.88 14.74 20 50 6.9 10.76 44.23 13.| __truncated__ ...
        .. .. .. .. ..$ doublet_scores    : num [1:436] NA 0.237 NA 0.375 0.152 0.666 0.134 0.248 NA 0.452 NA NA NA 0.57 0.251 NA NA NA 0.354 NA 0.442 NA NA NA 0.241 ...
        .. .. .. .. ..$ doublet_class     : Factor w/ 2 levels "singlet","doublet": NA 1 NA 1 1 1 1 1 NA 1 NA NA NA 1 1 NA NA NA 1 NA 1 NA NA NA 1 1 NA 1 1 NA 1 NA NA NA NA NA 1 NA NA 1 NA 1 1 1 1 1 NA NA NA NA ...
        .. .. .. .. ..$ emptyDrops_Total  : int [1:436] 208 422 49 520 578 329 637 442 229 521 119 59 42 594 441 80 26 58 223 52 262 149 133 164 314 743 68 345 593 83 35| __truncated__ ...
        .. .. .. .. ..$ emptyDrops_LogProb: num [1:436] -173 -587 NA -333 -667 -376 -843 -548 -269 -289 -204 NA NA -842 -609 NA NA NA -351 NA -348 -300 -237 -275 -409 ...
        .. .. .. .. ..$ emptyDrops_PValue : num [1:436] 0.9235 0.0001 NA 0.1191 0.0001 0.0001 0.0001 0.0001 0.0028 0.7043 0.0011 NA NA 0.0001 0.0001 NA NA NA 0.0001 NA 0| __truncated__ ...
        .. .. .. .. ..$ emptyDrops_Limited: logi [1:436] FALSE TRUE NA FALSE TRUE TRUE TRUE TRUE FALSE FALSE FALSE NA NA TRUE TRUE NA NA NA TRUE NA TRUE TRUE TRUE TRUE TR| __truncated__ ...
        .. .. .. .. ..$ emptyDrops_FDR    : num [1:436] 0.988069 0.000141 NA 0 0 0.000141 0 0.000141 0.003363 0 0.001344 NA NA 0 0.000141 NA NA NA 0.000141 NA 0.000141 0| __truncated__ ...
        .. .. .. ..@ active.assay: chr "RNA"
        .. .. .. ..@ active.ident: Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..- attr(*, "names")= chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" "GGGCTCACAATGTCTG-1" "TATCTTGCAGGCGTTC-1" "TTGGGCGCATTGACCA-1" "TCATGGAAGAAGATCT-1" "GGGTTATAGAATTGTG-1" ...
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
        .. .. .. .. .. .. .. .. ..@ i       : int [1:32918] 23 177 207 341 451 470 663 834 838 952 1016 1129 1316 1386 1434 1448 1479 1676 1694 1745 1871 1964 1969 1995 19 2| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:461] 0 24 104 163 258 309 343 451 490 532 546 559 634 673 792 867 961 1043 1106 1257 1347 1390 1471 1547 1641 1755 186| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 460
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" "TTTATGCAGAGGATCC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:32918] 1 1 3 1 2 5 1 1 1 1 1 1 1 1 1 1 1 2 9 1 1 6 1 2 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:32918] 23 177 207 341 451 470 663 834 838 952 1016 1129 1316 1386 1434 1448 1479 1676 1694 1745 1871 1964 1969 1995 19 2| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:461] 0 24 104 163 258 309 343 451 490 532 546 559 634 673 792 867 961 1043 1106 1257 1347 1390 1471 1547 1641 1755 186| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 460
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" "TTTATGCAGAGGATCC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:32918] 1 1 3 1 2 5 1 1 1 1 1 1 1 1 1 1 1 2 9 1 1 6 1 2 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ scale.data   : num[0 , 0 ] 
        .. .. .. .. .. .. ..@ key          : chr "rna_"
        .. .. .. .. .. .. ..@ assay.orig   : NULL
        .. .. .. .. .. .. ..@ var.features : logi(0) 
        .. .. .. .. .. .. ..@ meta.features:'data.frame':	2000 obs. of  0 variables
        .. .. .. .. .. .. ..@ misc         : list()
        .. .. .. ..@ meta.data   :'data.frame':	460 obs. of  13 variables:
        .. .. .. .. ..$ barcode           : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" "TTTATGCAGAGGATCC-1" ...
        .. .. .. .. ..$ orig.ident        : Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..$ nCount_RNA        : num [1:460] 46 289 249 263 193 93 472 101 112 41 38 297 81 372 272 174 406 107 529 238 93 232 165 300 411 373 533 155 310 553| __truncated__ ...
        .. .. .. .. ..$ nFeature_RNA      : int [1:460] 24 80 59 95 51 34 108 39 42 14 13 75 39 119 75 94 82 63 151 90 43 81 76 94 114 108 121 42 85 135 46 12 68 97 47 4| __truncated__ ...
        .. .. .. .. ..$ samples           : chr [1:460] "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" ...
        .. .. .. .. ..$ percent.mt        : num [1:460] 19.57 12.8 25.3 30.8 7.77 20.43 18.22 20.79 14.29 31.71 36.84 9.09 38.27 19.09 18.01 18.39 22.91 14.02 31.57 28.1| __truncated__ ...
        .. .. .. .. ..$ doublet_scores    : num [1:460] NA 0.0458 0.5706 0.966 NA NA 0.5052 NA NA NA NA 0.3721 NA 0.664 0.1622 NA 0.564 NA 0.9349 0.1361 NA 0.0261 NA 0.0479 0.1194 ...
        .. .. .. .. ..$ doublet_class     : Factor w/ 2 levels "singlet","doublet": NA 1 1 1 NA NA 1 NA NA NA NA 1 NA 1 1 NA 1 NA 1 1 NA 1 NA 1 1 1 2 NA 1 2 NA NA 1 1 NA NA 1 1 NA 1 2 NA 1 1 1 NA 1 NA 1 NA ...
        .. .. .. .. ..$ emptyDrops_Total  : int [1:460] 46 289 249 263 193 93 472 101 112 41 38 297 81 372 272 174 406 107 529 238 93 232 165 300 411 373 533 155 310 553| __truncated__ ...
        .. .. .. .. ..$ emptyDrops_LogProb: num [1:460] NA -367 -268 -404 -268 NA -490 -172 -198 NA NA -336 NA -500 -312 -381 -365 -274 -640 -372 NA -325 -311 -403 -476 ...
        .. .. .. .. ..$ emptyDrops_PValue : num [1:460] NA 0.0001 0.0009 0.0001 0.0001 NA 0.0001 0.0016 0.0002 NA NA 0.0001 NA 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 | __truncated__ ...
        .. .. .. .. ..$ emptyDrops_Limited: logi [1:460] NA TRUE FALSE TRUE TRUE NA TRUE FALSE FALSE NA NA TRUE NA TRUE TRUE TRUE TRUE TRUE TRUE TRUE NA TRUE TRUE TRUE TR| __truncated__ ...
        .. .. .. .. ..$ emptyDrops_FDR    : num [1:460] NA 0 0 0 0.000124 NA 0 0.001789 0.000239 NA NA 0 NA 0 0 0.000124 0 0.000124 0 0 NA 0 0.000124 0 0 ...
        .. .. .. ..@ active.assay: chr "RNA"
        .. .. .. ..@ active.ident: Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..- attr(*, "names")= chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" "TTTATGCAGAGGATCC-1" ...
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

---

    Code
      str(res, vec.len = 20)
    Output
      List of 2
       $ data  : list()
       $ output:List of 7
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
        .. .. .. ..@ i       : int [1:36762] 66 454 1662 1694 1964 70 194 290 307 435 467 506 512 573 618 641 666 674 762 812 819 934 947 1058 1091 1103 1367 | __truncated__ ...
        .. .. .. ..@ p       : int [1:501] 0 5 41 179 201 204 270 421 503 699 829 836 841 897 954 1001 1022 1052 1057 1252 1392 1441 1453 1461 1486 1565 157| __truncated__ ...
        .. .. .. ..@ Dim     : int [1:2] 2000 500
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. ..$ : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GCAACCGTCCTTTGAT-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "CGTCAAATCGGAAGGT-1" "GCACGGTGTGGCTGAA-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "TGGATGTCATCGCTCT-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" ...
        .. .. .. ..@ x       : num [1:36762] 1 1 2 10 1 1 2 1 1 1 1 3 1 1 1 1 1 3 1 1 1 1 1 1 1 1 3 1 1 1 2 159 3 1 1 1 1 1 3 1 3 1 1 2 1 4 1 1 1 26 ...
        .. .. .. ..@ factors : list()
        .. ..$ mock_sample_2_id:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. ..@ i       : int [1:33141] 23 177 207 341 451 470 663 834 838 952 1016 1129 1316 1386 1434 1448 1479 1676 1694 1745 1871 1964 1969 1995 425 | __truncated__ ...
        .. .. .. ..@ p       : int [1:502] 0 24 30 110 169 264 315 349 457 496 538 552 565 640 679 798 873 967 1049 1112 1263 1353 1396 1477 1553 1647 1649 | __truncated__ ...
        .. .. .. ..@ Dim     : int [1:2] 2000 501
        .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. ..$ : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" ...
        .. .. .. ..@ x       : num [1:33141] 1 1 3 1 2 5 1 1 1 1 1 1 1 1 1 1 1 2 9 1 1 6 1 2 3 4 1 1 1 11 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. ..@ factors : list()
        ..$ annot         :'data.frame':	2000 obs. of  3 variables:
        .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        ..$ edrops        :List of 2
        .. ..$ mock_sample_1_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:500] "GGGAGTAAGGCTAACG-1" "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GCAACCGTCCTTTGAT-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "CGTCAAATCGGAAGGT-1" "GCACGGTGTGGCTGAA-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "TGGATGTCATCGCTCT-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" ...
        .. .. .. ..@ nrows          : int 500
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:963, 1] 6.16e-05 1.75e-04 6.16e-05 6.16e-05 1.75e-04 1.04e-03 1.37e-03 3.90e-04 1.75e-04 3.90e-04 6.16e-05 1.75e-04 6.16e| __truncated__ ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:963] "CERNA1" "GPR174" "DCLK2" "AC001226.2" "PAX5" "SMG6" "EDF1" "ZNF605" "ANKS4B" "ACTR3B" "AC004067.1" "HSD17B4" "AC007881.3" "KCNAB3" "PSMA3" "TMSB15A" "LIPE" "PTX3" "DCBLD1" "FUNDC1" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 68.2
        .. .. .. .. ..$ retain : num 515
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:500] 15 208 422 49 6 520 578 329 637 442 13 436 229 521 119 59 42 206 594 441 80 26 10 58 223 120 52 262 149 133 164 3| __truncated__ ...
        .. .. .. .. ..$ LogProb: num [1:500] NA -173.1 -586.5 NA NA -332.5 -666.9 -375.6 -843.2 -548.1 NA -118.3 -268.8 -288.8 -204.3 NA NA -89.2 -842.1 -608.| __truncated__ ...
        .. .. .. .. ..$ PValue : num [1:500] NA 0.9235 0.0001 NA NA 0.1191 0.0001 0.0001 0.0001 0.0001 NA 1 0.0028 0.7043 0.0011 NA NA 1 0.0001 0.0001 NA NA NA NA 0.0001 ...
        .. .. .. .. ..$ Limited: logi [1:500] NA FALSE TRUE NA NA FALSE TRUE TRUE TRUE TRUE NA FALSE FALSE FALSE FALSE NA NA FALSE TRUE TRUE NA NA NA NA TRUE F| __truncated__ ...
        .. .. .. .. ..$ FDR    : num [1:500] NA 0.988069 0.000141 NA NA 0 0 0.000141 0 0.000141 NA 1 0.003363 0 0.001344 NA NA 1 0 0.000141 NA NA NA NA 0.000141 ...
        .. ..$ mock_sample_2_id:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
        .. .. .. ..@ rownames       : chr [1:501] "AATGGCTCAGTCCGTG-1" "GCAGCTGAGCTTAAGA-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" ...
        .. .. .. ..@ nrows          : int 501
        .. .. .. ..@ elementType    : chr "ANY"
        .. .. .. ..@ elementMetadata: NULL
        .. .. .. ..@ metadata       :List of 5
        .. .. .. .. ..$ lower  : num 100
        .. .. .. .. ..$ niters : num 10000
        .. .. .. .. ..$ ambient: num [1:893, 1] 8.09e-05 8.09e-05 8.09e-05 8.09e-05 2.60e-04 6.13e-04 2.09e-03 1.58e-04 8.09e-05 1.58e-04 1.58e-04 1.58e-04 2.60e| __truncated__ ...
        .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
        .. .. .. .. .. .. ..$ : chr [1:893] "CERNA1" "GPR174" "AC001226.2" "AC080129.2" "PAX5" "SMG6" "EDF1" "ZNF605" "LINC02413" "ACTR3B" "HSD17B4" "AC007881.3" "KCNAB3" "PSMA3" "LIPE" "FOXJ1" "PTX3" "DCBLD1" "FUNDC1" "KCNE1" ...
        .. .. .. .. .. .. ..$ : NULL
        .. .. .. .. ..$ alpha  : num 85.4
        .. .. .. .. ..$ retain : num 195
        .. .. .. ..@ listData       :List of 5
        .. .. .. .. ..$ Total  : int [1:501] 46 21 289 249 263 193 93 472 101 112 41 38 297 81 372 272 174 406 107 529 238 93 232 165 300 64 411 373 533 155 3| __truncated__ ...
        .. .. .. .. ..$ LogProb: num [1:501] NA NA -367 -268 -404 -268 NA -490 -172 -198 NA NA -336 NA -500 -312 -381 -365 -274 -640 -372 NA -325 -311 -403 ...
        .. .. .. .. ..$ PValue : num [1:501] NA NA 0.0001 0.0009 0.0001 0.0001 NA 0.0001 0.0016 0.0002 NA NA 0.0001 NA 0.0001 0.0001 0.0001 0.0001 0.0001 0.00| __truncated__ ...
        .. .. .. .. ..$ Limited: logi [1:501] NA NA TRUE FALSE TRUE TRUE NA TRUE FALSE FALSE NA NA TRUE NA TRUE TRUE TRUE TRUE TRUE TRUE TRUE NA TRUE TRUE TRUE| __truncated__ ...
        .. .. .. .. ..$ FDR    : num [1:501] NA NA 0 0 0 0.000124 NA 0 0.001789 0.000239 NA NA 0 NA 0 0 0.000124 0 0.000124 0 0 NA 0 0.000124 0 ...
        ..$ doublet_scores:List of 2
        .. ..$ mock_sample_1_id:'data.frame':	216 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:216] "GGTCACGAGTGAGCCA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "TTACCATCACTGTTCC-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" "TCATGGAAGAAGATCT-1" "CTCCAACGTGATTAGA-1" "AGTAACCAGTATTGCC-1" "TGTGCGGTCATTACTC-1" "AGATGCTTCGTTGCCT-1" "TGATGGTCAACCGATT-1" "AACCTTTCAGTGGTGA-1" "TATCCTAGTGATACTC-1" "TTTACGTAGCGATCGA-1" "GTGCTTCCAAATCGTC-1" "CCAATTTCAGGCTTGC-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 ...
        .. .. ..$ doublet_scores: num [1:216] 0.2368 0.3751 0.152 0.6659 0.1344 0.2483 0.4523 0.5696 0.2508 0.3544 0.4421 0.2415 0.1291 0.078 0.3669 0.0889 0.8| __truncated__ ...
        .. ..$ mock_sample_2_id:'data.frame':	242 obs. of  3 variables:
        .. .. ..$ barcodes      : chr [1:242] "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "CAGCGTGGTACGATCT-1" "GCAGGCTCAAACTGCT-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "TGACAGTAGGCATGCA-1" "CAAGGGACAGAGCGTA-1" "TTTATGCAGAGGATCC-1" "AGCTCAAAGGCCCGTT-1" "TCGGGACGTAACATGA-1" "AATGACCCAGAGTCAG-1" "TCGTAGAGTGAGTCAG-1" "TCTTAGTAGGGATGTC-1" "ATCCTATCAAGACCTT-1" "CGATCGGGTGTGATGG-1" "CCTCTCCAGCGTGTTT-1" "TTCGATTGTTCCTACC-1" "CACGTTCTCGACCTAA-1" ...
        .. .. ..$ doublet_class : Factor w/ 2 levels "singlet","doublet": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 2 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. ..$ doublet_scores: num [1:242] 0.0458 0.5706 0.966 0.5052 0.3721 0.664 0.1622 0.564 0.9349 0.1361 0.0261 0.0479 0.1194 0.8091 0.9956 0.2044 0.99| __truncated__ ...
        ..$ scdata_list   :List of 2
        .. ..$ mock_sample_2_id:Formal class 'Seurat' [package "SeuratObject"] with 13 slots
        .. .. .. ..@ assays      :List of 1
        .. .. .. .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
        .. .. .. .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:32918] 23 177 207 341 451 470 663 834 838 952 1016 1129 1316 1386 1434 1448 1479 1676 1694 1745 1871 1964 1969 1995 19 2| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:461] 0 24 104 163 258 309 343 451 490 532 546 559 634 673 792 867 961 1043 1106 1257 1347 1390 1471 1547 1641 1755 186| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 460
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" "TTTATGCAGAGGATCC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:32918] 1 1 3 1 2 5 1 1 1 1 1 1 1 1 1 1 1 2 9 1 1 6 1 2 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:32918] 23 177 207 341 451 470 663 834 838 952 1016 1129 1316 1386 1434 1448 1479 1676 1694 1745 1871 1964 1969 1995 19 2| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:461] 0 24 104 163 258 309 343 451 490 532 546 559 634 673 792 867 961 1043 1106 1257 1347 1390 1471 1547 1641 1755 186| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 460
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" "TTTATGCAGAGGATCC-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:32918] 1 1 3 1 2 5 1 1 1 1 1 1 1 1 1 1 1 2 9 1 1 6 1 2 2 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ scale.data   : num[0 , 0 ] 
        .. .. .. .. .. .. ..@ key          : chr "rna_"
        .. .. .. .. .. .. ..@ assay.orig   : NULL
        .. .. .. .. .. .. ..@ var.features : logi(0) 
        .. .. .. .. .. .. ..@ meta.features:'data.frame':	2000 obs. of  0 variables
        .. .. .. .. .. .. ..@ misc         : list()
        .. .. .. ..@ meta.data   :'data.frame':	460 obs. of  14 variables:
        .. .. .. .. ..$ barcode           : chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" "TTTATGCAGAGGATCC-1" ...
        .. .. .. .. ..$ orig.ident        : Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..$ nCount_RNA        : num [1:460] 46 289 249 263 193 93 472 101 112 41 38 297 81 372 272 174 406 107 529 238 93 232 165 300 411 373 533 155 310 553| __truncated__ ...
        .. .. .. .. ..$ nFeature_RNA      : int [1:460] 24 80 59 95 51 34 108 39 42 14 13 75 39 119 75 94 82 63 151 90 43 81 76 94 114 108 121 42 85 135 46 12 68 97 47 4| __truncated__ ...
        .. .. .. .. ..$ samples           : chr [1:460] "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" "mock_sample_2_id" ...
        .. .. .. .. ..$ percent.mt        : num [1:460] 19.57 12.8 25.3 30.8 7.77 20.43 18.22 20.79 14.29 31.71 36.84 9.09 38.27 19.09 18.01 18.39 22.91 14.02 31.57 28.1| __truncated__ ...
        .. .. .. .. ..$ doublet_scores    : num [1:460] NA 0.0458 0.5706 0.966 NA NA 0.5052 NA NA NA NA 0.3721 NA 0.664 0.1622 NA 0.564 NA 0.9349 0.1361 NA 0.0261 NA 0.0479 0.1194 ...
        .. .. .. .. ..$ doublet_class     : Factor w/ 2 levels "singlet","doublet": NA 1 1 1 NA NA 1 NA NA NA NA 1 NA 1 1 NA 1 NA 1 1 NA 1 NA 1 1 1 2 NA 1 2 NA NA 1 1 NA NA 1 1 NA 1 2 NA 1 1 1 NA 1 NA 1 NA ...
        .. .. .. .. ..$ emptyDrops_Total  : int [1:460] 46 289 249 263 193 93 472 101 112 41 38 297 81 372 272 174 406 107 529 238 93 232 165 300 411 373 533 155 310 553| __truncated__ ...
        .. .. .. .. ..$ emptyDrops_LogProb: num [1:460] NA -367 -268 -404 -268 NA -490 -172 -198 NA NA -336 NA -500 -312 -381 -365 -274 -640 -372 NA -325 -311 -403 -476 ...
        .. .. .. .. ..$ emptyDrops_PValue : num [1:460] NA 0.0001 0.0009 0.0001 0.0001 NA 0.0001 0.0016 0.0002 NA NA 0.0001 NA 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 | __truncated__ ...
        .. .. .. .. ..$ emptyDrops_Limited: logi [1:460] NA TRUE FALSE TRUE TRUE NA TRUE FALSE FALSE NA NA TRUE NA TRUE TRUE TRUE TRUE TRUE TRUE TRUE NA TRUE TRUE TRUE TR| __truncated__ ...
        .. .. .. .. ..$ emptyDrops_FDR    : num [1:460] NA 0 0 0 0.000124 NA 0 0.001789 0.000239 NA NA 0 NA 0 0 0.000124 0 0.000124 0 0 NA 0 0.000124 0 0 ...
        .. .. .. .. ..$ cells_id          : int [1:460] 560 320 152 73 227 145 633 48 127 302 23 838 355 600 164 621 531 409 296 882 282 620 516 211 859 258 313 480 297 | __truncated__ ...
        .. .. .. ..@ active.assay: chr "RNA"
        .. .. .. ..@ active.ident: Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..- attr(*, "names")= chr [1:460] "AATGGCTCAGTCCGTG-1" "TACGGGCGTTATAGCC-1" "AAGACAACATCAACCA-1" "GGACGTCCACTGCTTC-1" "GTAAGTCCACGAGAAC-1" "GTCTCACAGTGTAGTA-1" "CAGCGTGGTACGATCT-1" "TCGTGGGTCGTGTCAA-1" "AGAAGCGAGTGGCAGT-1" "CACAACAAGTTGGCTT-1" "CCGTAGGCATAGTCGT-1" "GCAGGCTCAAACTGCT-1" "AACCTGAGTAACTTCG-1" "GGTGGCTAGAAGCCTG-1" "AGGCCACTCATGAGTC-1" "GTCGAATCACAGACGA-1" "TGACAGTAGGCATGCA-1" "ATTCTACGTGTCCAAT-1" "CAAGGGACAGAGCGTA-1" "TTTATGCAGAGGATCC-1" ...
        .. .. .. ..@ graphs      : list()
        .. .. .. ..@ neighbors   : list()
        .. .. .. ..@ reductions  : list()
        .. .. .. ..@ images      : list()
        .. .. .. ..@ project.name: chr "mock_experiment"
        .. .. .. ..@ misc        :List of 2
        .. .. .. .. ..$ gene_annotations:'data.frame':	2000 obs. of  3 variables:
        .. .. .. .. .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
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
        .. .. .. .. .. .. .. .. ..@ i       : int [1:36426] 70 194 290 307 435 467 506 512 573 618 641 666 674 762 812 819 934 947 1058 1091 1103 1367 1400 1571 1624 1665 16| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:437] 0 36 174 196 262 413 495 691 821 877 934 981 1002 1032 1227 1367 1416 1428 1453 1532 1559 1637 1709 1758 1819 190| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 436
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" "GGGCTCACAATGTCTG-1" "TATCTTGCAGGCGTTC-1" "TTGGGCGCATTGACCA-1" "TCATGGAAGAAGATCT-1" "GGGTTATAGAATTGTG-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:36426] 1 2 1 1 1 1 3 1 1 1 1 1 3 1 1 1 1 1 1 1 1 3 1 1 1 2 159 3 1 1 1 1 1 3 1 3 1 1 2 1 4 1 1 1 26 7 1 1 2 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
        .. .. .. .. .. .. .. .. ..@ i       : int [1:36426] 70 194 290 307 435 467 506 512 573 618 641 666 674 762 812 819 934 947 1058 1091 1103 1367 1400 1571 1624 1665 16| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ p       : int [1:437] 0 36 174 196 262 413 495 691 821 877 934 981 1002 1032 1227 1367 1416 1428 1453 1532 1559 1637 1709 1758 1819 190| __truncated__ ...
        .. .. .. .. .. .. .. .. ..@ Dim     : int [1:2] 2000 436
        .. .. .. .. .. .. .. .. ..@ Dimnames:List of 2
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. .. .. .. .. ..$ : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" "GGGCTCACAATGTCTG-1" "TATCTTGCAGGCGTTC-1" "TTGGGCGCATTGACCA-1" "TCATGGAAGAAGATCT-1" "GGGTTATAGAATTGTG-1" ...
        .. .. .. .. .. .. .. .. ..@ x       : num [1:36426] 1 2 1 1 1 1 3 1 1 1 1 1 3 1 1 1 1 1 1 1 1 3 1 1 1 2 159 3 1 1 1 1 1 3 1 3 1 1 2 1 4 1 1 1 26 7 1 1 2 1 ...
        .. .. .. .. .. .. .. .. ..@ factors : list()
        .. .. .. .. .. .. ..@ scale.data   : num[0 , 0 ] 
        .. .. .. .. .. .. ..@ key          : chr "rna_"
        .. .. .. .. .. .. ..@ assay.orig   : NULL
        .. .. .. .. .. .. ..@ var.features : logi(0) 
        .. .. .. .. .. .. ..@ meta.features:'data.frame':	2000 obs. of  0 variables
        .. .. .. .. .. .. ..@ misc         : list()
        .. .. .. ..@ meta.data   :'data.frame':	436 obs. of  14 variables:
        .. .. .. .. ..$ barcode           : chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" "GGGCTCACAATGTCTG-1" "TATCTTGCAGGCGTTC-1" "TTGGGCGCATTGACCA-1" "TCATGGAAGAAGATCT-1" "GGGTTATAGAATTGTG-1" ...
        .. .. .. .. ..$ orig.ident        : Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..$ nCount_RNA        : num [1:436] 208 422 49 520 578 329 637 442 229 521 119 59 42 594 441 80 26 58 223 52 262 149 133 164 314 743 68 345 593 83 35| __truncated__ ...
        .. .. .. .. ..$ nFeature_RNA      : int [1:436] 36 138 22 66 151 82 196 130 56 57 47 21 30 195 140 49 12 25 79 27 78 72 49 61 88 204 17 116 168 26 125 26 56 57 1| __truncated__ ...
        .. .. .. .. ..$ samples           : chr [1:436] "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" "mock_sample_1_id" ...
        .. .. .. .. ..$ percent.mt        : num [1:436] 76.44 25.36 53.06 64.23 11.25 8.81 27.16 21.04 66.81 85.6 20.17 59.32 28.57 20.88 14.74 20 50 6.9 10.76 44.23 13.| __truncated__ ...
        .. .. .. .. ..$ doublet_scores    : num [1:436] NA 0.237 NA 0.375 0.152 0.666 0.134 0.248 NA 0.452 NA NA NA 0.57 0.251 NA NA NA 0.354 NA 0.442 NA NA NA 0.241 ...
        .. .. .. .. ..$ doublet_class     : Factor w/ 2 levels "singlet","doublet": NA 1 NA 1 1 1 1 1 NA 1 NA NA NA 1 1 NA NA NA 1 NA 1 NA NA NA 1 1 NA 1 1 NA 1 NA NA NA NA NA 1 NA NA 1 NA 1 1 1 1 1 NA NA NA NA ...
        .. .. .. .. ..$ emptyDrops_Total  : int [1:436] 208 422 49 520 578 329 637 442 229 521 119 59 42 594 441 80 26 58 223 52 262 149 133 164 314 743 68 345 593 83 35| __truncated__ ...
        .. .. .. .. ..$ emptyDrops_LogProb: num [1:436] -173 -587 NA -333 -667 -376 -843 -548 -269 -289 -204 NA NA -842 -609 NA NA NA -351 NA -348 -300 -237 -275 -409 ...
        .. .. .. .. ..$ emptyDrops_PValue : num [1:436] 0.9235 0.0001 NA 0.1191 0.0001 0.0001 0.0001 0.0001 0.0028 0.7043 0.0011 NA NA 0.0001 0.0001 NA NA NA 0.0001 NA 0| __truncated__ ...
        .. .. .. .. ..$ emptyDrops_Limited: logi [1:436] FALSE TRUE NA FALSE TRUE TRUE TRUE TRUE FALSE FALSE FALSE NA NA TRUE TRUE NA NA NA TRUE NA TRUE TRUE TRUE TRUE TR| __truncated__ ...
        .. .. .. .. ..$ emptyDrops_FDR    : num [1:436] 0.988069 0.000141 NA 0 0 0.000141 0 0.000141 0.003363 0 0.001344 NA NA 0 0.000141 NA NA NA 0.000141 NA 0.000141 0| __truncated__ ...
        .. .. .. .. ..$ cells_id          : int [1:436] 216 894 758 202 851 518 475 14 761 399 671 822 341 446 542 141 631 218 294 423 455 718 370 176 391 259 380 10 710| __truncated__ ...
        .. .. .. ..@ active.assay: chr "RNA"
        .. .. .. ..@ active.ident: Factor w/ 1 level "mock_experiment": 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
        .. .. .. .. ..- attr(*, "names")= chr [1:436] "GGCAGTCAGGCCTTGC-1" "GGTCACGAGTGAGCCA-1" "ATCCCTGTCACCATGA-1" "GGTTAACTCATATGGC-1" "TACCTCGTCGACCCAG-1" "GATGAGGCAGGCATTT-1" "CGTAATGTCATACGAC-1" "AGTGACTCAGGGTCTC-1" "TTCACCGAGCAAATCA-1" "TTACCATCACTGTTCC-1" "TCAGTGACAGGCTATT-1" "CATCGTCGTTTCGCTC-1" "CTCCCTCAGCACTGGA-1" "GTTCGCTAGGGACCAT-1" "TATCGCCAGTCCCAAT-1" "GGGCTCACAATGTCTG-1" "TATCTTGCAGGCGTTC-1" "TTGGGCGCATTGACCA-1" "TCATGGAAGAAGATCT-1" "GGGTTATAGAATTGTG-1" ...
        .. .. .. ..@ graphs      : list()
        .. .. .. ..@ neighbors   : list()
        .. .. .. ..@ reductions  : list()
        .. .. .. ..@ images      : list()
        .. .. .. ..@ project.name: chr "mock_experiment"
        .. .. .. ..@ misc        :List of 2
        .. .. .. .. ..$ gene_annotations:'data.frame':	2000 obs. of  3 variables:
        .. .. .. .. .. ..$ input        : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. ..$ name         : chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. .. ..$ original_name: chr [1:2000] "CERNA1" "IGKV2-18" "LUADT1" "AL133297.2" "AL021937.1" "GPR174" "DCLK2" "C16orf90" "PCSK6-AS1" "AC001226.2" "POU6F2" "AC080129.2" "AC015908.3" "CYP2A6" "PAX5" "COL6A5" "SMG6" "AL645634.1" "LINC02372" "EDF1" ...
        .. .. .. .. ..$ experimentId    : chr "mock_experiment_id"
        .. .. .. ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
        .. .. .. .. ..$ : int [1:3] 4 1 0
        .. .. .. ..@ commands    : list()
        .. .. .. ..@ tools       :List of 1
        .. .. .. .. ..$ flag_filtered: logi FALSE
        ..$ qc_config     :List of 7
        .. ..$ cellSizeDistribution:List of 2
        .. .. ..$ mock_sample_2_id:List of 4
        .. .. .. ..$ enabled              : logi FALSE
        .. .. .. ..$ auto                 : logi TRUE
        .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. ..$ minCellSize: num 17
        .. .. .. .. ..$ binStep    : num 200
        .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. ..$ minCellSize: num 17
        .. .. .. .. ..$ binStep    : num 200
        .. .. ..$ mock_sample_1_id:List of 4
        .. .. .. ..$ enabled              : logi FALSE
        .. .. .. ..$ auto                 : logi TRUE
        .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. ..$ minCellSize: num 27
        .. .. .. .. ..$ binStep    : num 200
        .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. ..$ minCellSize: num 27
        .. .. .. .. ..$ binStep    : num 200
        .. ..$ mitochondrialContent:List of 2
        .. .. ..$ mock_sample_2_id:List of 4
        .. .. .. ..$ enabled              : logi TRUE
        .. .. .. ..$ auto                 : logi TRUE
        .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. ..$ maxFraction: num 0.521
        .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. ..$ maxFraction: num 0.521
        .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. .. ..$ mock_sample_1_id:List of 4
        .. .. .. ..$ enabled              : logi TRUE
        .. .. .. ..$ auto                 : logi TRUE
        .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. ..$ maxFraction: num 0.603
        .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. ..$ maxFraction: num 0.603
        .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. ..$ classifier          :List of 2
        .. .. ..$ mock_sample_2_id:List of 5
        .. .. .. ..$ enabled              : logi TRUE
        .. .. .. ..$ prefiltered          : logi FALSE
        .. .. .. ..$ auto                 : logi TRUE
        .. .. .. ..$ filterSettings       :List of 1
        .. .. .. .. ..$ FDR: num 0.01
        .. .. .. ..$ defaultFilterSettings:List of 1
        .. .. .. .. ..$ FDR: num 0.01
        .. .. ..$ mock_sample_1_id:List of 5
        .. .. .. ..$ enabled              : logi TRUE
        .. .. .. ..$ prefiltered          : logi FALSE
        .. .. .. ..$ auto                 : logi TRUE
        .. .. .. ..$ filterSettings       :List of 1
        .. .. .. .. ..$ FDR: num 0.01
        .. .. .. ..$ defaultFilterSettings:List of 1
        .. .. .. .. ..$ FDR: num 0.01
        .. ..$ numGenesVsNumUmis   :List of 2
        .. .. ..$ mock_sample_2_id:List of 4
        .. .. .. ..$ enabled              : logi TRUE
        .. .. .. ..$ auto                 : logi TRUE
        .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. ..$ mock_sample_1_id:List of 4
        .. .. .. ..$ enabled              : logi TRUE
        .. .. .. ..$ auto                 : logi TRUE
        .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. ..$ p.level: num 0.001
        .. ..$ doubletScores       :List of 2
        .. .. ..$ mock_sample_2_id:List of 4
        .. .. .. ..$ enabled              : logi TRUE
        .. .. .. ..$ auto                 : logi TRUE
        .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. ..$ probabilityThreshold: num 0.979
        .. .. .. .. ..$ binStep             : num 0.02
        .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. ..$ probabilityThreshold: num 0.979
        .. .. .. .. ..$ binStep             : num 0.02
        .. .. ..$ mock_sample_1_id:List of 4
        .. .. .. ..$ enabled              : logi TRUE
        .. .. .. ..$ auto                 : logi TRUE
        .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. ..$ probabilityThreshold: num 0.84
        .. .. .. .. ..$ binStep             : num 0.02
        .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. ..$ probabilityThreshold: num 0.84
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
      str(res, vec.len = 20)
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
        .. .. .. ..$ mock_sample_2_id:List of 4
        .. .. .. .. ..$ enabled              : logi FALSE
        .. .. .. .. ..$ auto                 : logi TRUE
        .. .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. .. ..$ minCellSize: num 17
        .. .. .. .. .. ..$ binStep    : num 200
        .. .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. .. ..$ minCellSize: num 17
        .. .. .. .. .. ..$ binStep    : num 200
        .. .. .. ..$ mock_sample_1_id:List of 4
        .. .. .. .. ..$ enabled              : logi FALSE
        .. .. .. .. ..$ auto                 : logi TRUE
        .. .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. .. ..$ minCellSize: num 27
        .. .. .. .. .. ..$ binStep    : num 200
        .. .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. .. ..$ minCellSize: num 27
        .. .. .. .. .. ..$ binStep    : num 200
        .. .. ..$ mitochondrialContent:List of 2
        .. .. .. ..$ mock_sample_2_id:List of 4
        .. .. .. .. ..$ enabled              : logi TRUE
        .. .. .. .. ..$ auto                 : logi TRUE
        .. .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. .. ..$ maxFraction: num 0.521
        .. .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. .. ..$ maxFraction: num 0.521
        .. .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. .. .. ..$ mock_sample_1_id:List of 4
        .. .. .. .. ..$ enabled              : logi TRUE
        .. .. .. .. ..$ auto                 : logi TRUE
        .. .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. .. ..$ maxFraction: num 0.603
        .. .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. .. ..$ method        : chr "absoluteThreshold"
        .. .. .. .. .. ..$ methodSettings:List of 1
        .. .. .. .. .. .. ..$ absoluteThreshold:List of 2
        .. .. .. .. .. .. .. ..$ maxFraction: num 0.603
        .. .. .. .. .. .. .. ..$ binStep    : num 0.3
        .. .. ..$ classifier          :List of 2
        .. .. .. ..$ mock_sample_2_id:List of 5
        .. .. .. .. ..$ enabled              : logi TRUE
        .. .. .. .. ..$ prefiltered          : logi FALSE
        .. .. .. .. ..$ auto                 : logi TRUE
        .. .. .. .. ..$ filterSettings       :List of 1
        .. .. .. .. .. ..$ FDR: num 0.01
        .. .. .. .. ..$ defaultFilterSettings:List of 1
        .. .. .. .. .. ..$ FDR: num 0.01
        .. .. .. ..$ mock_sample_1_id:List of 5
        .. .. .. .. ..$ enabled              : logi TRUE
        .. .. .. .. ..$ prefiltered          : logi FALSE
        .. .. .. .. ..$ auto                 : logi TRUE
        .. .. .. .. ..$ filterSettings       :List of 1
        .. .. .. .. .. ..$ FDR: num 0.01
        .. .. .. .. ..$ defaultFilterSettings:List of 1
        .. .. .. .. .. ..$ FDR: num 0.01
        .. .. ..$ numGenesVsNumUmis   :List of 2
        .. .. .. ..$ mock_sample_2_id:List of 4
        .. .. .. .. ..$ enabled              : logi TRUE
        .. .. .. .. ..$ auto                 : logi TRUE
        .. .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. ..$ mock_sample_1_id:List of 4
        .. .. .. .. ..$ enabled              : logi TRUE
        .. .. .. .. ..$ auto                 : logi TRUE
        .. .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. .. ..$ regressionType        : chr "linear"
        .. .. .. .. .. ..$ regressionTypeSettings:List of 2
        .. .. .. .. .. .. ..$ linear:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. .. .. .. .. ..$ spline:List of 1
        .. .. .. .. .. .. .. ..$ p.level: num 0.001
        .. .. ..$ doubletScores       :List of 2
        .. .. .. ..$ mock_sample_2_id:List of 4
        .. .. .. .. ..$ enabled              : logi TRUE
        .. .. .. .. ..$ auto                 : logi TRUE
        .. .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. .. ..$ probabilityThreshold: num 0.979
        .. .. .. .. .. ..$ binStep             : num 0.02
        .. .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. .. ..$ probabilityThreshold: num 0.979
        .. .. .. .. .. ..$ binStep             : num 0.02
        .. .. .. ..$ mock_sample_1_id:List of 4
        .. .. .. .. ..$ enabled              : logi TRUE
        .. .. .. .. ..$ auto                 : logi TRUE
        .. .. .. .. ..$ filterSettings       :List of 2
        .. .. .. .. .. ..$ probabilityThreshold: num 0.84
        .. .. .. .. .. ..$ binStep             : num 0.02
        .. .. .. .. ..$ defaultFilterSettings:List of 2
        .. .. .. .. .. ..$ probabilityThreshold: num 0.84
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

