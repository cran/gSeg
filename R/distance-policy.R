# Fast paths are enabled only for exact, verified stats-library/R-wrapper
# identities. Unknown builds use stats::dist(); no architecture-only guess.
.gseg_dist_state <- new.env(parent = emptyenv())
.gseg_dist_build_info <- function() {
  tryCatch({
    path <- getLoadedDLLs()[["stats"]][["path"]]
    list(r = R.version.string, platform = R.version$platform,
         stats_md5 = unname(suppressWarnings(tools::md5sum(path))),
         body_text = paste(deparse(body(stats::dist), width.cutoff = 500), collapse = "\n"),
         formals_text = paste(deparse(formals(stats::dist), width.cutoff = 500), collapse = "\n"))
  }, error = function(e) NULL)
}
.gseg_dist_match_build <- function(info, builds = .gseg_dist_builds) {
  fields <- c("r", "platform", "stats_md5", "body_text", "formals_text")
  if (is.null(info)) return(NA_integer_)
  for (build in builds) {
    if (all(vapply(fields, function(k) identical(info[[k]], build[[k]]), logical(1))))
      return(if (identical(build$mode, "fused")) 1L else 0L)
  }
  NA_integer_
}

# The default matrix conversion is an identity only while its generic,
# default method and applicable method tables remain standard. Do not invoke
# custom/active bindings while deciding whether to use the reference path.
.gseg_dist_conversion_ok <- function() {
  tryCatch({
    b <- asNamespace("base")
    for (name in c("as.matrix", "as.matrix.default", "is.matrix", "UseMethod", ".__S3MethodsTable__."))
      if (bindingIsActive(name, b)) return(FALSE)
    generic <- get("as.matrix", b, inherits = FALSE)
    default <- get("as.matrix.default", b, inherits = FALSE)
    template_generic <- function(x, ...) UseMethod("as.matrix")
    template_default <- function(x, ...) {
      if (is.matrix(x)) x else array(x, c(length(x), 1L),
                                    if (!is.null(names(x))) list(names(x), NULL) else NULL)
    }
    if (!identical(environment(generic), b) || !identical(environment(default), b) ||
        !identical(body(generic), body(template_generic)) ||
        !identical(formals(generic), formals(template_generic)) ||
        !identical(body(default), body(template_default)) ||
        !identical(formals(default), formals(template_default)) ||
        !identical(get("is.matrix", b, inherits = FALSE), .Primitive("is.matrix")) ||
        !identical(get("UseMethod", b, inherits = FALSE), .Primitive("UseMethod")) ||
        isdebugged(generic) || isdebugged(default)) return(FALSE)
    names <- paste0("as.matrix.", c("matrix", "array", "double", "integer", "numeric"))
    scopes <- list(.GlobalEnv, asNamespace("stats"))
    for (scope in scopes) {
      for (name in names)
        if (exists(name, envir = scope, inherits = TRUE)) return(FALSE)
      env <- scope
      repeat {
        if (exists("as.matrix.default", envir = env, inherits = FALSE)) {
          if (bindingIsActive("as.matrix.default", env) ||
              !identical(get("as.matrix.default", env, inherits = FALSE), default)) return(FALSE)
          break
        }
        if (identical(env, emptyenv())) return(FALSE)
        env <- parent.env(env)
      }
    }
    table <- get(".__S3MethodsTable__.", b, inherits = FALSE)
    if (!is.environment(table)) return(FALSE)
    for (name in names)
      if (exists(name, table, inherits = FALSE)) return(FALSE)
    if (exists("as.matrix.default", table, inherits = FALSE) &&
        (bindingIsActive("as.matrix.default", table) ||
         !identical(get("as.matrix.default", table, inherits = FALSE), default))) return(FALSE)
    TRUE
  }, error = function(e) FALSE)
}

.gseg_dist_policy <- function() {
  if (exists("mode", envir = .gseg_dist_state, inherits = FALSE)) {
    mode <- .gseg_dist_state$mode
    if (is.na(mode) || !identical(stats::dist, .gseg_dist_state$reference_function) ||
        !.gseg_dist_conversion_ok()) return(NA_integer_)
    return(mode)
  }
  mode <- .gseg_dist_match_build(.gseg_dist_build_info())
  if (!is.na(mode) && !.gseg_dist_conversion_ok()) mode <- NA_integer_
  if (!is.na(mode)) {
    x <- matrix(sin(seq_len(8L * 31L)), 8L)
    value <- tryCatch(.Call(C_gseg_dist_reference, x, 8L, 31L, mode, TRUE),
                      error = function(e) NULL)
    reference <- stats::dist(x)
    attributes(reference) <- NULL
    if (!identical(reference, value, num.eq = FALSE, single.NA = FALSE))
      mode <- NA_integer_
  }
  .gseg_dist_state$reference_function <- stats::dist
  .gseg_dist_state$mode <- mode
  mode
}
.gseg_reference_dist <- function(x) {
  # Preserve all user-defined conversion and class behavior in stats::dist.
  if (is.object(x) || !(typeof(x) %in% c("integer", "double")))
    return(stats::dist(x))
  dims <- attr(x, "dim", exact = TRUE)
  if (length(dims) != 2L) return(stats::dist(x))
  n <- dims[1L]
  p <- dims[2L]
  if (n < 16L || p < 1L || length(x) > .Machine$integer.max - n)
    return(stats::dist(x))
  mode <- .gseg_dist_policy()
  if (is.na(mode)) return(stats::dist(x))
  # Bound the additional layout copy by the unavoidable distance output.
  rowmajor <- length(x) <= n * (n - 1.0) / 2.0
  value <- .Call(C_gseg_dist_reference, x, n, p, mode, rowmajor)
  if (is.null(value)) return(stats::dist(x))
  structure(value, Size = n, Labels = rownames(x), Diag = FALSE,
            Upper = FALSE, method = "euclidean", call = quote(stats::dist(x = x)),
            class = "dist")
}
