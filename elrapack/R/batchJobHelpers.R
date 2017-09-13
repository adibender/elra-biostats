#' Extracts string from simulation log files, following a specified pattern.
#'
#' @import BatchJobs
#' @param reg A registry object (see BatchJobs documentation)
#' @param ids A subset of job ids
#' @param pattern Greps Log files for Lines with this pattern.
#' @keywords internal
#' @rdname get_convergence_info
#' @export
get_pattern <- function(
    reg,
    ids,
    pattern=".*Algorithm parameters.*form=*") {

    all.files       <- getLogFiles(reg, ids=ids)
    logs.read       <- lapply(all.files, readLines)
    lines.with.name <- lapply(logs.read, function(z)
        unlist(lapply(z, function(x) grep(pattern, x, value = TRUE))))

    # return
    sub("(.*)\\=", "", unlist(lines.with.name))

}


#' Extract number of iterations from Log files
#'
#' @rdname get_convergence_info
#' @inherit get_pattern
get_maxIterations_perNode <- function(reg, ids) {

  all.files       <- getLogFiles(reg, ids=ids)
  logs.read       <- lapply(all.files, readLines)
  iteration.lines <- lapply(logs.read, function(z)
    unlist(lapply(z, function(x) grep("*Iterations*", x, value = TRUE))))

  n.iterations <- sapply(iteration.lines, length)

  data.frame(name=get_pattern(reg, ids=ids), iterations = n.iterations)

}

#' Extract runtime per job
#'
#' @rdname get_convergence_info
#' @inherit get_pattern
#' @importFrom stringr str_split
get_runtime_perNode <- function(reg, ids) {

    all.files <- getLogFiles(reg, ids=ids)
    logs.read <- lapply(all.files, readLines)
    last.line <- sapply(logs.read, function(z) z[length(z)])
    times     <- lapply(str_split(last.line, " "), as.numeric)
    elapsed   <- sapply(times, max, na.rm = TRUE)

    # return
    data.frame(name = get_pattern(reg, ids=ids), elapsed = elapsed)

}


#' Extract Info like run-time, convergence etc. for each job
#'
#' Wrapper function that obtains iterations, convergence, runtime
#' and node name for specific Jobs (\code{ids}).
#'
#' @rdname get_convergence_info
#' @inherit get_pattern
get_convergence_info <- function(reg, ids) {

    not.converged <- grepLogs(
        reg,
        ids     = ids,
        pattern = "*algorithm did not converge*")

    conv <- data.frame(
        iterations = get_maxIterations_perNode(reg, ids=ids)$iterations,
        converged = factor("yes", levels = c("yes", "no")),
        time      = getJobInfo(reg, ids  = ids)$time.running,
        node      = getJobInfo(reg, ids  = ids)$nodename)

    conv$converged[ids %in% not.converged] <- "no"
    conv$ids <- ids

    # return
    conv

}


