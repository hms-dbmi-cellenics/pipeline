require("remoter")
try(
    remoter::server(port=6969)
)

sleep = Sys.getenv("REMOTER_DEBUG_SLEEP", "");
if (sleep != "") {
    sleep = as.integer(sleep)
    message('About to sleep ', sleep, ' seconds before exiting.')
    Sys.sleep(sleep)
}
message('Exiting...')
