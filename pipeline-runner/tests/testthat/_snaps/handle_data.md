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

