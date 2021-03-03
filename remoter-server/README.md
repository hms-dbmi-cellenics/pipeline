remoter-server
===============

This is a simple container running a remoter server for a persistent R kernel.

How to run (individually)
-------------------------

You can launch an individual client instance with details of a particular task if
you do not want to run the entire pipeline.

Make sure you have the server running locally. Then, you can launch the client
with the argument being the JSON-ified string of the task that is to be executed.

For example:

    docker run --rm -ti remoter-client "$(cat test.json)"
