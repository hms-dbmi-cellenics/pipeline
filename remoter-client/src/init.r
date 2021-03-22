require("RJSONIO")
require("remoter")
require("ids")

PORT=6969

# get request from the arguments of the docker run command, and escape all quotes.
request <- commandArgs(trailingOnly = TRUE)[1]
run_id <- ids::random_id()
message("Got request with ID", run_id, "...")
message(request)
message("")
parsed = RJSONIO::fromJSON(request)
if (parsed$server == "host.docker.internal") {
    parsed$server = Sys.getenv("DOCKER_GATEWAY_HOST")
    if (parsed$server == "") {
        parsed$server = "host.docker.internal"
    }
}

# Get sample ids
sample_id = Sys.getenv("SAMPLE_ID", "")

# load wrapper in case it changed from last run
message("Loading wrapper for server ", parsed$server, "...")
remoter::batch(addr = parsed$server, port = PORT, file = "./wrapper.r")

message('')
message('Copying request...')
message(sprintf("c2s(request, 'request_%s')", run_id))
remoter::batch(addr = parsed$server, port = PORT, script = sprintf("c2s(request, 'request_%s')", run_id))

message(sprintf('Launching work for sample %s...', sample_id))
message(sprintf("wrapper(request_%s, '%s')", run_id, sample_id))
remoter::batch(addr = parsed$server, port = PORT, script = sprintf("wrapper(request_%s, '%s')", run_id, sample_id))
message('Exiting...')
