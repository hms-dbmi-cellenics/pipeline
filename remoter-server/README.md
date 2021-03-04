remoter-server
===============

This is a simple container running a remoter server for a persistent R kernel.

Task definitions for individual pipeline items are located under `src/`.

How to run
----------

Simply run:

    docker-compose up --build

from this directory.

You can also open the pipeline project in a development container through VS Code,
which will also start the server.