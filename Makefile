#!make
#----------------------------------------
# Settings
#----------------------------------------
.DEFAULT_GOAL := help
#--------------------------------------------------
# Variables
#--------------------------------------------------
# If unix name is not Darwin assume we are on Linux and add this needed env variable
# See README.md/Running on Docker issues, for more info
ifneq ($(shell uname -s), Darwin)
ifndef DOCKER_GATEWAY_HOST
	export DOCKER_GATEWAY_HOST=`docker network inspect --format='{{range .IPAM.Config}}{{.Gateway}}{{end}}'
endif
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
