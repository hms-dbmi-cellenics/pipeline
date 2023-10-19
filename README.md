[![codecov](https://codecov.io/gh/hms-dbmi-cellenics/pipeline/branch/master/graph/badge.svg?token=kQ19q1EenW)](https://codecov.io/gh/hms-dbmi-cellenics/pipeline)
# Cellenics Pipeline

The Cellenics pipeline project for dependency-managed work processing.

## Getting started

The steps of the pipeline that are run through this project are started
spontaneously on your machine as Docker containers, simulating Kubernetes
in local development.

We have included a utility so you can automatically monitor containers spawned
and read their logs as they are executing.

For local development, you should already have Docker and Node.js installed, as well as
Inframock running.

Afterwards, you can install the pipeline dependencies with:

```bash
make install
```

To build and run the pipeline:

```bash
make build && make run
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
make build
```

## Local development and adding dependencies

First make sure the project library is synchronized with the lockfile:

```R
# inside pipeline-runner folder
renv::restore()
```

**NOTE**: To restore Bioconductor packages your R version needs to be the same as in the [Dockerfile](pipeline-runner/Dockerfile) (4.2.0).

`install.packages(...)` and use them (e.g. `dplyr::left_join(...)`) as you normally would. Then, update the lockfile:

```R
renv::snapshot()
```

commit the changes to the lockfile (used to install dependencies in the Dockerfile). See [renv docs](https://rstudio.github.io/renv/) for more info.

### Development dependencies

Packages used for interactive development, such as `devtools`, `usethis`, `roxygen2`,
`styler` and the R `languageserver` (to develop R in vscode!) and their dependencies
should not be added to the lockfile, since they are not required at runtime. 
`renv` has been configured to ignore them. 

To install them, run the following block, with no arguments. This installs all
packages in the DESCRIPTION file, which includes the development dependencies in
the Suggests section.

```R
renv::install()
```

### Running tests locally

There are several ways to run tests locally. The easiest one being using the Rstudio
shortcut `Cmd + shift + T`.

Other ways to run tests locally:

```R
devtools::test()
```

```R
testthat::test_local()
```

### snapshot tests

There are two special test files, `test-gem2s.R` and `test-qc.R`. They run sequential
snapshot tests for gem2s and qc respectively. They output, for each step, a summarised
text representation of the current pipeline state, and a hash of the complete R
object. For more information on the workings of snapshot testing, see [testthat snapshots](https://testthat.r-lib.org/articles/snapshotting.html).

In case of a test failure, you need to manually examine the difference between
the new output and the original snapshot. To do so, either use your favorite `diff`
tool, or the following `testthat` function, which opens an interactive diff tool:

```R
testthat::snapshot_review()
```

To accept the new changes, you can either click on the `Accept` button in the interactive
diff tool opened by [`snapshot_review()`](https://testthat.r-lib.org/reference/snapshot_accept.html),
or run the following command:

```R
testthat::snapshot_accept()
```

#### Adding test datasets to the snapshot tests

The tests were written trying to simplify the way to add new test datasets.
Basically, two things are required:

1. Add sample files to the `pipeline-runner/tests/testthat/mock_data/originals_bucket` folder
2. Create an `input.json` file, named like `<experiment_id>-input.json` and put it
   in the `pipeline-runner/tests/testthatmock_data/input` folder. It should 
   contain the correct sample file ids (the name of the sample files in the
   `originals_bucket` folder) and `experiment_id`


## Debugging locally

**TLDR**: save something inside `/debug` in a data processing or gem2s step to
 access it later from `./local-runner/debug`.

 **TLDR2**: if the pipeline throws an error, `tryCatchLog` will save a [dump file](https://github.com/aryoda/tryCatchLog#how-do-i-perform-a-post-mortem-analysis-of-my-crashed-r-script) in  `./local-runner/debug` that can be used for inspecting the workspace and object values along the call stack.

To save the parameters (`config`, `seurat_obj`, etc) to a **data processing** task function, specify `DEBUG_STEP`.
Available tasks include all task names listed in `run_processing_step` [init.R](pipeline-runner/init.R#L69) as well as `DEBUG_STEP=all` 
to save the parameters to all data processing task functions:

```bash
# e.g. DEBUG_STEP=dataIntegration
DEBUG_STEP=task_name make run
```

When the pipeline is run, it will save the parameters to the specified `task_name` in `$(pwd)/debug`. You
can load these into your R environment:

```R
# clicking the file in RStudio does this for you
load('{task_name}_{sample_id}.RData')

# if you need to load multiple tasks, you can load each into a seperate environment
# you would when access objects using e.g. task_env$scdata
task_env <- new.env()
load('{task_name}_{sample_id}.RData', envir = task_env)
```

## Troubleshooting

#### Linux Mint 20.3 Cinnamon
```bash
Error in curl::curl_fetch_memory(url, handle = handle) : 
Timeout was reached: [172.17.0.1:4566] Connection timeout after 60001 ms
Calls: init ... request_fetch -> request_fetch.write_memory -> <Anonymous>
Execution halted
```
Turn off firewall or allow incoming traffic. This would allow AWS to send packages to the pipeline, which would otherwise be blocked by the firewall.

1. Open *Firewall Configuration* from the Start Menu.
2. Select **Allow** in the **Outgoing** dropdown menu (Alternatively, set **Status** to OFF).

#### Javascript not connecting to the docker daemon
When running the pipeline locally you may get the following error: 
```bash
Creating mock Lambda function on InfraMock...
Stack with ARN arn:aws:cloudformation:eu-west-1:000000000000:stack/local-container-launcher/f693a088 successfully created.
Waiting for Docker events...
pipeline/local-runner/src/app.js:88
    containers.forEach((info) => {
               ^

TypeError: Cannot read property 'forEach' of null
    at Object.callback (pipeline/local-runner/src/app.js:88:16)
    at pipeline/local-runner/node_modules/dockerode/lib/docker.js:428:12
    at Modem.buildPayload (pipeline/local-runner/node_modules/docker-modem/lib/modem.js:297:19)
    at ClientRequest.<anonymous> (pipeline/local-runner/node_modules/docker-modem/lib/modem.js:282:10)
```

To solve this you need to Allow the default Docker socket to be used by third parties:

1. Open 'Settings' in the Docker desktop app (GUI)
2. Click Advanced
3. Make sure the setting 'Allow the default Docker socket to be used' is turned on
