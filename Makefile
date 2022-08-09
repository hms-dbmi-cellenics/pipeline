#!make
#----------------------------------------
# Settings
#----------------------------------------
.DEFAULT_GOAL := help
#--------------------------------------------------
# Variables
#--------------------------------------------------
# If unix name is not Darwin assume we are on Linux and add this needed env variable
# to be picked up by the R worker.
# See README.md/Running on Docker issues, for more info
ifneq ($(shell uname -s), Darwin)
	# Get the gateway address of the default bridge network
	export HOST_IP=$(shell docker network inspect bridge --format='{{(index .IPAM.Config 0).Gateway}}')
endif
#--------------------------------------------------
# Targets
#--------------------------------------------------
install: 
	@echo "Installing local runner"
	@(cd ./local-runner && npm install)
	@echo "Installing R env packages"
	@(cd ./pipeline-runner && R -e "renv::restore()")
build: 
    # regenerate sysdata.rda env file
	@(cd ./pipeline-runner && Rscript data-raw/sysdata.R)
	@(cd ./local-runner && npm run build)
test:
	@(cd ./pipeline-runner && R -e "devtools::test()")
run: build run-only
run-only:
	@(cd ./local-runner && npm start)
.PHONY: install build run run-only help
help: ## Shows available targets
	@fgrep -h "## " $(MAKEFILE_LIST) | fgrep -v fgrep | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-13s\033[0m %s\n", $$1, $$2}'
