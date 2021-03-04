remoter-client
===============

This is a simple container running a remoter client to run R script.

How to run (individually)
-------------------------

You can launch an individual client instance with details of a particular task if
you do not want to run the entire pipeline.

Make sure you have the server running locally. Refer to the server README on how to do this.

Then, you can launch the client with the argument being the JSON-ified string
of the task that is to be executed. A sample task is supplied in this folder. You can check
the relevant API code under [here](https://github.com/biomage-ltd/api/tree/master/src/api/general-services/pipeline-manage/constructors) to see the relevant code for generating the schema. In addition, the schema
is also available in the API repository spec folder.

For example:

    docker run --rm -ti remoter-client "$(cat sample_task_input.json)"
