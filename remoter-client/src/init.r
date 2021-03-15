require("remoter")
require("ids")

# get request from the arguments of the docker run command, and escape all quotes.
request <- commandArgs(trailingOnly = TRUE)[1]
run_id <- ids::random_id()
message("Got request with ID", run_id, "...")
message(request)
message("")
parsed = RJSONIO::fromJSON(request)

# load wrapper in case it changed from last run
message("Loading wrapper...")
remoter::batch(addr = parsed$server, port = 6969, file = "./wrapper.r")

message('')
message('Copying request...')
message(sprintf("c2s(request, 'request_%s')", run_id))
remoter::batch(addr = "host.docker.internal", port = 6969, script = sprintf("c2s(request, 'request_%s')", run_id))

message('Launching work...')
message(sprintf("wrapper(request_%s)", run_id))
remoter::batch(addr = "host.docker.internal", port = 6969, script = sprintf("wrapper(request_%s)", run_id))
