# send_pipeline_fail_update handles a gem2s call successefully

    Code
      str(mockery::mock_args(mock_publish))
    Output
      List of 1
       $ :List of 3
        ..$ Message          : chr "[\n \"test_experiment\",\n\"test_task\",\n{\n \"processName\": \"gem2s\",\n\"experimentId\": \"test_experiment\"| __truncated__
        ..$ TopicArn         : chr "test_topic"
        ..$ MessageAttributes:List of 1
        .. ..$ type:List of 3
        .. .. ..$ DataType   : chr "String"
        .. .. ..$ StringValue: chr "GEM2SResponse"
        .. .. ..$ BinaryValue: NULL

# send_pipeline_fail_update handles a qc call successfully with no global_env config

    Code
      mockery::mock_args(mock_publish)
    Output
      [[1]]
      [[1]]$Message
      [1] "{\n \"experimentId\": null,\n\"taskName\": null,\n\"input\": {\n \"processName\": \"qc\" \n},\n\"response\": {\n \"error\": \"qc\" \n},\n\"pipelineVersion\":        2,\n\"apiUrl\": \"test_url\" \n}"
      
      [[1]]$TopicArn
      [1] "test_topic"
      
      [[1]]$MessageAttributes
      [[1]]$MessageAttributes$type
      [[1]]$MessageAttributes$type$DataType
      [1] "String"
      
      [[1]]$MessageAttributes$type$StringValue
      [1] "PipelineResponse"
      
      [[1]]$MessageAttributes$type$BinaryValue
      NULL
      
      
      
      

# send_pipeline_fail_update handles a qc call successfully with global_env config

    Code
      mockery::mock_args(mock_publish)
    Output
      [[1]]
      [[1]]$Message
      [1] "{\n \"experimentId\": null,\n\"taskName\": \"mitochondrialContent\",\n\"input\": {\n \"processName\": \"qc\",\n\"sampleUuid\": \"00000000-0000-0000-000000000111\",\n\"taskName\": \"mitochondrialContent\" \n},\n\"response\": {\n \"error\": \"qc\" \n},\n\"pipelineVersion\":        2,\n\"apiUrl\": \"test_url\",\n\"output\": {\n \"bucket\": \"test_bucket\",\n\"key\": \"mock-uuid\" \n} \n}"
      
      [[1]]$TopicArn
      [1] "test_topic"
      
      [[1]]$MessageAttributes
      [[1]]$MessageAttributes$type
      [[1]]$MessageAttributes$type$DataType
      [1] "String"
      
      [[1]]$MessageAttributes$type$StringValue
      [1] "PipelineResponse"
      
      [[1]]$MessageAttributes$type$BinaryValue
      NULL
      
      
      
      

---

    Code
      mockery::mock_args(mock_put_object_in_s3)
    Output
      [[1]]
      [[1]][[1]]
      [[1]][[1]]$aws_config
      list()
      
      [[1]][[1]]$results_bucket
      [1] "test_bucket"
      
      [[1]][[1]]$api_url
      [1] "test_url"
      
      [[1]][[1]]$sns_topic
      [1] "test_topic"
      
      
      [[1]][[2]]
      [1] "test_bucket"
      
      [[1]][[3]]
        [1] 7b 0a 20 22 63 6f 6e 66 69 67 22 3a 20 7b 0a 20 22 65 6d 62 65 64 64 69 6e
       [26] 67 53 65 74 74 69 6e 67 73 22 3a 20 7b 0a 20 22 6d 65 74 68 6f 64 22 3a 20
       [51] 22 75 6d 61 70 22 2c 0a 22 6d 65 74 68 6f 64 53 65 74 74 69 6e 67 73 22 3a
       [76] 20 7b 0a 20 22 74 73 6e 65 22 3a 20 7b 0a 20 22 70 65 72 70 6c 65 78 69 74
      [101] 79 22 3a 20 20 20 20 20 39 2e 31 38 2c 0a 22 6c 65 61 72 6e 69 6e 67 52 61
      [126] 74 65 22 3a 20 20 20 20 20 20 32 30 30 20 0a 7d 2c 0a 22 75 6d 61 70 22 3a
      [151] 20 7b 0a 20 22 64 69 73 74 61 6e 63 65 4d 65 74 72 69 63 22 3a 20 22 63 6f
      [176] 73 69 6e 65 22 2c 0a 22 6d 69 6e 69 6d 75 6d 44 69 73 74 61 6e 63 65 22 3a
      [201] 20 20 20 20 20 20 30 2e 33 20 0a 7d 20 0a 7d 20 0a 7d 2c 0a 22 63 6c 75 73
      [226] 74 65 72 69 6e 67 53 65 74 74 69 6e 67 73 22 3a 20 7b 0a 20 22 6d 65 74 68
      [251] 6f 64 22 3a 20 22 6c 6f 75 76 61 69 6e 22 2c 0a 22 6d 65 74 68 6f 64 53 65
      [276] 74 74 69 6e 67 73 22 3a 20 7b 0a 20 22 6c 6f 75 76 61 69 6e 22 3a 20 7b 0a
      [301] 20 22 72 65 73 6f 6c 75 74 69 6f 6e 22 3a 20 20 20 20 20 20 30 2e 38 20 0a
      [326] 7d 20 0a 7d 20 0a 7d 20 0a 7d 20 0a 7d
      
      [[1]][[4]]
      [1] "mock-uuid"
      
      

# unflatten_cell_sets works

    Code
      cell_sets
    Output
      [[1]]
      [[1]]$key
      [1] "louvain"
      
      [[1]]$name
      [1] "fake louvain clusters"
      
      [[1]]$rootNode
      [1] TRUE
      
      [[1]]$type
      [1] "cellSets"
      
      [[1]]$children
      [[1]]$children[[1]]
      [[1]]$children[[1]]$key
      [1] "louvain-0"
      
      [[1]]$children[[1]]$name
      [1] "Cluster 0"
      
      [[1]]$children[[1]]$rootNode
      [1] FALSE
      
      [[1]]$children[[1]]$color
      [1] "#e377c2"
      
      [[1]]$children[[1]]$type
      [1] "cellSets"
      
      [[1]]$children[[1]]$cellIds
       [1]  1  2  3  4  5  6  7  8  9 10 11 12 16 17 18 19 20 21 22 23 24 25 26 27 28
      [26] 29 30
      
      
      [[1]]$children[[2]]
      [[1]]$children[[2]]$key
      [1] "louvain-1"
      
      [[1]]$children[[2]]$name
      [1] "Cluster 1"
      
      [[1]]$children[[2]]$rootNode
      [1] FALSE
      
      [[1]]$children[[2]]$color
      [1] "#8c564b"
      
      [[1]]$children[[2]]$type
      [1] "cellSets"
      
      [[1]]$children[[2]]$cellIds
       [1] 31 32 33 34 35 36 37 38 39 40 41 42 43 47 48 49 50 51 52 53 54 55
      
      
      [[1]]$children[[3]]
      [[1]]$children[[3]]$key
      [1] "louvain-2"
      
      [[1]]$children[[3]]$name
      [1] "Cluster 2"
      
      [[1]]$children[[3]]$rootNode
      [1] FALSE
      
      [[1]]$children[[3]]$color
      [1] "#d62728"
      
      [[1]]$children[[3]]$type
      [1] "cellSets"
      
      [[1]]$children[[3]]$cellIds
      [1] 57 58 59 60
      
      
      [[1]]$children[[4]]
      [[1]]$children[[4]]$key
      [1] "louvain-3"
      
      [[1]]$children[[4]]$name
      [1] "Cluster 3"
      
      [[1]]$children[[4]]$rootNode
      [1] FALSE
      
      [[1]]$children[[4]]$color
      [1] "#2ca02c"
      
      [[1]]$children[[4]]$type
      [1] "cellSets"
      
      [[1]]$children[[4]]$cellIds
      [1] 61 62 64 65
      
      
      [[1]]$children[[5]]
      [[1]]$children[[5]]$key
      [1] "louvain-4"
      
      [[1]]$children[[5]]$name
      [1] "Cluster 4"
      
      [[1]]$children[[5]]$rootNode
      [1] FALSE
      
      [[1]]$children[[5]]$color
      [1] "#ff7f0e"
      
      [[1]]$children[[5]]$type
      [1] "cellSets"
      
      [[1]]$children[[5]]$cellIds
      [1] 66 68 69 70
      
      
      [[1]]$children[[6]]
      [[1]]$children[[6]]$key
      [1] "louvain-5"
      
      [[1]]$children[[6]]$name
      [1] "Cluster 5"
      
      [[1]]$children[[6]]$rootNode
      [1] FALSE
      
      [[1]]$children[[6]]$color
      [1] "#1f77b4"
      
      [[1]]$children[[6]]$type
      [1] "cellSets"
      
      [[1]]$children[[6]]$cellIds
      [1] 71 72 73 75
      
      
      [[1]]$children[[7]]
      [[1]]$children[[7]]$key
      [1] "louvain-6"
      
      [[1]]$children[[7]]$name
      [1] "Cluster 6"
      
      [[1]]$children[[7]]$rootNode
      [1] FALSE
      
      [[1]]$children[[7]]$color
      [1] "#f8e71c"
      
      [[1]]$children[[7]]$type
      [1] "cellSets"
      
      [[1]]$children[[7]]$cellIds
      [1] 76 77 79 80
      
      
      
      
      [[2]]
      [[2]]$key
      [1] "scratchpad"
      
      [[2]]$name
      [1] "Custom cell sets"
      
      [[2]]$rootNode
      [1] TRUE
      
      [[2]]$type
      [1] "cellSets"
      
      [[2]]$children
      list()
      
      
      [[3]]
      [[3]]$key
      [1] "sample"
      
      [[3]]$name
      [1] "Samples"
      
      [[3]]$rootNode
      [1] TRUE
      
      [[3]]$type
      [1] "metadataCategorical"
      
      [[3]]$children
      [[3]]$children[[1]]
      [[3]]$children[[1]]$key
      [1] "sample-id-1"
      
      [[3]]$children[[1]]$name
      [1] "KO"
      
      [[3]]$children[[1]]$rootNode
      NULL
      
      [[3]]$children[[1]]$color
      [1] "#8c564b"
      
      [[3]]$children[[1]]$type
      NULL
      
      [[3]]$children[[1]]$cellIds
       [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
      [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
      
      
      [[3]]$children[[2]]
      [[3]]$children[[2]]$key
      [1] "sample-id-2"
      
      [[3]]$children[[2]]$name
      [1] "WT1"
      
      [[3]]$children[[2]]$rootNode
      NULL
      
      [[3]]$children[[2]]$color
      [1] "#d62728"
      
      [[3]]$children[[2]]$type
      NULL
      
      [[3]]$children[[2]]$cellIds
       [1] 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65
      [26] 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80
      
      
      
      
      [[4]]
      [[4]]$key
      [1] "Track_1"
      
      [[4]]$name
      [1] "Track_1"
      
      [[4]]$rootNode
      [1] TRUE
      
      [[4]]$type
      [1] "metadataCategorical"
      
      [[4]]$children
      [[4]]$children[[1]]
      [[4]]$children[[1]]$key
      [1] "Track_1-KMeta"
      
      [[4]]$children[[1]]$name
      [1] "KMeta"
      
      [[4]]$children[[1]]$rootNode
      NULL
      
      [[4]]$children[[1]]$color
      [1] "#8c564b"
      
      [[4]]$children[[1]]$type
      NULL
      
      [[4]]$children[[1]]$cellIds
       [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
      [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
      
      
      [[4]]$children[[2]]
      [[4]]$children[[2]]$key
      [1] "Track_1-WMetaT"
      
      [[4]]$children[[2]]$name
      [1] "WMetaT"
      
      [[4]]$children[[2]]$rootNode
      NULL
      
      [[4]]$children[[2]]$color
      [1] "#d62728"
      
      [[4]]$children[[2]]$type
      NULL
      
      [[4]]$children[[2]]$cellIds
       [1]  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59
      [20]  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78
      [39]  79  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97
      [58]  98  99 100
      
      
      
      

