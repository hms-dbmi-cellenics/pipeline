# Biomage Pipeline

The Cellscope pipeline project for dependency-managed work processing.

## Getting started

The steps of the pipeline that are run through this project are started
spontaneously on your machine as Docker containers, simulating Kubernetes
in local development.

We have included a utility so you can automatically monitor containers spawned
and read their logs as they are executing.

For local development, you should already have Docker and Node.js installed, as well as
Inframock running.

Afterwards, you can install and run local-runner.

```bash
# cd ./local-runner
# npm install
# npm run build
# npm start
```

A similar message should appear:

```
> node src/app.js

Loading CloudFormation for local container launcher...
Creating mock Lambda function on InfraMock...
No previous stack found on InfraMock.
Stack with ARN arn:aws:cloudformation:eu-west-1:000000000000:stack/local-container-launcher/106d1df9 successfully created.
Waiting for Docker events...
```

Logs from pipelines run through the API will apear here.

## Rebuiling the docker images

```bash
> npm run build
```

## Docker issues on Linux

The remoter client and server do not have their own private network, and they communicate through the host
using the special address `host.docker.internal`. This is not recognised by Docker on Linux. This address
can be overriden by the environment variable `DOCKER_GATEWAY_HOST`:

```bash
EXPORT DOCKER_GATEWAY_HOST=`docker network inspect --format='{{range .IPAM.Config}}{{.Gateway}}{{end}}'
npm start
```

## Debugging locally

To save the arguments (`config`, `scdata`, etc) to a task function, specify DEBUG_STEP and DEBUG_PATH.
Available tasks include all those defined in `run_step` of  [init.r](qc-runner/src/init.r) as well as `DEBUG_STEP=all` 
to save the arguments to all task functions:

```bash
# e.g. DEBUG_STEP=dataIntegration
DEBUG_STEP=task_name DEBUG_PATH=/path/to/debug/folder npm start
```
Note: `DEBUG_PATH` needs to be an absolute path (rather than relative), i.e. being in `pipeline/local-runner` this can be used 
to populate the local subfolder `mydebug`:

```bash
DEBUG_STEP=all DEBUG_PATH=${PWD}/mydebug npm restart
```

When the pipeline is run, it will save the arguments to the specified `task_name` in `DEBUG_PATH`. You
can load these into your R environment:

```R
# clicking the file in RStudio does this for you
load('/path/to/debug/folder/{task_name}_{sample_id}.RData')

# if you need to load multiple tasks, you can load each into a seperate environment
# you would when access objects using e.g. task_env$scdata
task_env <- new.env()
load('/path/to/debug/folder/{task_name}_{sample_id}.RData', envir = task_env)
```
