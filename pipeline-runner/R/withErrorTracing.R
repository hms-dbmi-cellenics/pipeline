#' Add Error Tracing to tryCatch
#'
#' @param expr expression to evaluate
#' @param skip_pattern grep pattern. Will not print calls including or after
#'  first match of \code{skip_pattern}
#'
withErrorTracing = function(expr, silentSuccess=FALSE, skip_pattern = '^multi_assign') {
    hasFailed = FALSE
    messages = list()
    warnings = list()

    errorTracer = function(obj) {

        # Storing the call stack
        calls = sys.calls()
        calls = calls[1:length(calls)-1]
        # Keeping the calls only
        trace = limitedLabels(c(calls, attr(obj, "calls")))

        # print what's useful from the calls
        trace <- rev(trace)
        is.trace <- grep(skip_pattern, trace)[1]-1

        cat('\n-------\n')
        cat('🚩 Backtrack:\n')
        for (i in 2:is.trace) message(trace[i])
        cat('-------\n\n')

        # Muffle any redundant output of the same message
        optionalRestart = function(r) { res = findRestart(r); if (!is.null(res)) invokeRestart(res) }
        optionalRestart("muffleMessage")
        optionalRestart("muffleWarning")
    }

    vexpr = withCallingHandlers(withVisible(expr),  error=errorTracer)
    if (silentSuccess && !hasFailed) {
        cat(paste(warnings, collapse=""))
    }
    if (vexpr$visible) vexpr$value else invisible(vexpr$value)
}
