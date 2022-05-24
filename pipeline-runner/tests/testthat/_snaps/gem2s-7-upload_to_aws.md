# get_cell_sets without metadata matches snapshot

    Code
      str(cell_sets)
    Output
      List of 1
       $ cellSets:List of 2
        ..$ :List of 5
        .. ..$ key     : chr "scratchpad"
        .. ..$ name    : chr "Custom cell sets"
        .. ..$ rootNode: logi TRUE
        .. ..$ children: list()
        .. ..$ type    : chr "cellSets"
        ..$ :List of 5
        .. ..$ key     : chr "sample"
        .. ..$ name    : chr "Samples"
        .. ..$ rootNode: logi TRUE
        .. ..$ children:List of 2
        .. .. ..$ :List of 4
        .. .. .. ..$ key    : chr "123abc"
        .. .. .. ..$ name   : chr "WT1"
        .. .. .. ..$ color  : chr "#77aadd"
        .. .. .. ..$ cellIds: int [1:40] 0 1 2 3 4 5 6 7 8 9 ...
        .. .. ..$ :List of 4
        .. .. .. ..$ key    : chr "123def"
        .. .. .. ..$ name   : chr "WT2"
        .. .. .. ..$ color  : chr "#ee8866"
        .. .. .. ..$ cellIds: int [1:40] 40 41 42 43 44 45 46 47 48 49 ...
        .. ..$ type    : chr "metadataCategorical"

# get_cell_sets with two metadata groups matches snapshot

    Code
      str(cell_sets)
    Output
      List of 1
       $ cellSets:List of 4
        ..$ :List of 5
        .. ..$ key     : chr "scratchpad"
        .. ..$ name    : chr "Custom cell sets"
        .. ..$ rootNode: logi TRUE
        .. ..$ children: list()
        .. ..$ type    : chr "cellSets"
        ..$ :List of 5
        .. ..$ key     : chr "sample"
        .. ..$ name    : chr "Samples"
        .. ..$ rootNode: logi TRUE
        .. ..$ children:List of 2
        .. .. ..$ :List of 4
        .. .. .. ..$ key    : chr "123abc"
        .. .. .. ..$ name   : chr "WT1"
        .. .. .. ..$ color  : chr "#77aadd"
        .. .. .. ..$ cellIds: int [1:40] 0 1 2 3 4 5 6 7 8 9 ...
        .. .. ..$ :List of 4
        .. .. .. ..$ key    : chr "123def"
        .. .. .. ..$ name   : chr "WT2"
        .. .. .. ..$ color  : chr "#ee8866"
        .. .. .. ..$ cellIds: int [1:40] 40 41 42 43 44 45 46 47 48 49 ...
        .. ..$ type    : chr "metadataCategorical"
        ..$ :List of 5
        .. ..$ key     : chr "Group1"
        .. ..$ name    : chr "Group1"
        .. ..$ rootNode: logi TRUE
        .. ..$ children:List of 2
        .. .. ..$ :List of 4
        .. .. .. ..$ key    : chr "Group1-Hello"
        .. .. .. ..$ name   : chr "Hello"
        .. .. .. ..$ color  : chr "#eedd88"
        .. .. .. ..$ cellIds: int [1:40] 0 1 2 3 4 5 6 7 8 9 ...
        .. .. ..$ :List of 4
        .. .. .. ..$ key    : chr "Group1-WT2"
        .. .. .. ..$ name   : chr "WT2"
        .. .. .. ..$ color  : chr "#ffaabb"
        .. .. .. ..$ cellIds: int [1:40] 40 41 42 43 44 45 46 47 48 49 ...
        .. ..$ type    : chr "metadataCategorical"
        ..$ :List of 5
        .. ..$ key     : chr "Group2"
        .. ..$ name    : chr "Group2"
        .. ..$ rootNode: logi TRUE
        .. ..$ children:List of 1
        .. .. ..$ :List of 4
        .. .. .. ..$ key    : chr "Group2-WT"
        .. .. .. ..$ name   : chr "WT"
        .. .. .. ..$ color  : chr "#99ddff"
        .. .. .. ..$ cellIds: int [1:80] 0 1 2 3 4 5 6 7 8 9 ...
        .. ..$ type    : chr "metadataCategorical"

