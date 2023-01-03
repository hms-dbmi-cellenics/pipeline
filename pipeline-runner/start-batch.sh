#!/bin/bash -e
python3 ./datadog-batch/monitor.py &
Rscript init.R