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
      
      

