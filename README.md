# Biomage Pipeline
The Cellscope pipeline project for dependency-managed work processing.

Getting started
---------------

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
$ EXPORT DOCKER_GATEWAY_HOST=`docker network inspect --format='{{range .IPAM.Config}}{{.Gateway}}{{end}}'
$ npm start
``` 

**TODO**: Another OSX/Linux incompatibility is related to the aws development configuration, where OSX expects
an url while Linux expects a host name.