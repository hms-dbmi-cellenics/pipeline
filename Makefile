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
	export HOST_IP=$(shell hostname -I | awk '{print $$1}')
endif
#--------------------------------------------------
# Targets
#--------------------------------------------------
install: 
	@(cd ./local-runner && npm install)
build: 
	@(cd ./local-runner && npm run build)
run:
	@(cd ./local-runner && npm start)
.PHONY: install build run help
help: ## Shows available targets
	@fgrep -h "## " $(MAKEFILE_LIST) | fgrep -v fgrep | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-13s\033[0m %s\n", $$1, $$2}'
