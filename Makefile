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
	@echo "Installing renv packages"
	@(cd ./pipeline-runner && R -e "renv::restore()")
update-sysdata: 
    # regenerate sysdata.rda env file
    # this step depends on your local R installation to run 
	@(cd ./pipeline-runner && Rscript data-raw/sysdata.R)	
build: 
	@(cd ./local-runner && npm run build)
build-batch-staging: 
	@(cd pipeline-runner && docker build --target batch --tag 242905224710.dkr.ecr.eu-west-1.amazonaws.com/pipeline:batch-staging --build-arg GITHUB_PAT=${GITHUB_API_TOKEN} .)
	@aws ecr get-login-password --region 'eu-west-1' | docker login --username AWS --password-stdin 242905224710.dkr.ecr.eu-west-1.amazonaws.com
	@docker push 242905224710.dkr.ecr.eu-west-1.amazonaws.com/pipeline:batch-staging
test: build ## Executes unit tests (usage: make test)
	@docker run \
		--entrypoint /bin/bash \
		-v $(PWD)/pipeline-runner/tests/testthat/_snaps:/src/pipeline-runner/tests/testthat/_snaps \
		biomage-pipeline-runner \
		-c "R -e 'testthat::test_local(stop_on_failure = FALSE)'"
test-file: build ## Tests a specific test file (usage: make test-file FILE=test-file.R)
	@docker run \
		--entrypoint /bin/bash \
		-v $(PWD)/pipeline-runner/tests/testthat/_snaps:/src/pipeline-runner/tests/testthat/_snaps \
		biomage-pipeline-runner \
		-c "R -e \"pkgload::load_all(); testthat::test_file('tests/testthat/$(FILE)')\""
snap-accept: build ## Accept updated snaps (usage: make snap-accept)
	@docker run \
		--entrypoint /bin/bash \
		-v $(PWD)/pipeline-runner/tests/testthat/_snaps:/src/pipeline-runner/tests/testthat/_snaps \
		biomage-pipeline-runner \
		-c "R -e \"testthat::snapshot_accept()\""
hooks: ## Configures path to git hooks
	@git config core.hooksPath .githooks
run: build run-only
run-only:
	@(cd ./local-runner && npm start)
.PHONY: install build run run-only help
help: ## Shows available targets
	@fgrep -h "## " $(MAKEFILE_LIST) | fgrep -v fgrep | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-13s\033[0m %s\n", $$1, $$2}'
