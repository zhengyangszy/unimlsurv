# Small internal helpers -----------------------------------------------------
#
# `unimlsurv` mirrors the layering of `uniml` but keeps its own assertions, so a
# change in one package cannot silently change the other one's error behaviour.

.abort <- function(...) {
  stop(..., call. = FALSE)
}

.assert_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .abort(sprintf("`%s` must be a single non-missing character value.", name))
  }
  invisible(TRUE)
}

.assert_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .abort(sprintf("`%s` must be a single non-missing logical value.", name))
  }
  invisible(TRUE)
}

.assert_scalar_number <- function(x, name, integer = FALSE, min = -Inf,
                                  max = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    .abort(sprintf("`%s` must be a single finite numeric value.", name))
  }
  if (integer && x != round(x)) {
    .abort(sprintf("`%s` must be an integer value.", name))
  }
  if (x < min || x > max) {
    .abort(sprintf("`%s` must be between %s and %s.", name, format(min), format(max)))
  }
  invisible(TRUE)
}

.assert_character_vector <- function(x, name, allow_empty = TRUE) {
  if (!is.character(x) || anyNA(x)) {
    .abort(sprintf("`%s` must be a character vector without missing values.", name))
  }
  if (!allow_empty && length(x) == 0L) {
    .abort(sprintf("`%s` must not be empty.", name))
  }
  invisible(TRUE)
}

.assert_choice <- function(x, name, choices) {
  .assert_scalar_character(x, name)
  if (!x %in% choices) {
    .abort(sprintf(
      "`%s` must be one of: %s.",
      name,
      paste(sprintf("\"%s\"", choices), collapse = ", ")
    ))
  }
  invisible(TRUE)
}

.assert_numeric_vector <- function(x, name, min = -Inf, max = Inf,
                                   allow_na = FALSE, allow_empty = FALSE) {
  if (!is.numeric(x)) {
    .abort(sprintf("`%s` must be a numeric vector.", name))
  }
  if (!allow_empty && length(x) == 0L) {
    .abort(sprintf("`%s` must not be empty.", name))
  }
  if (!allow_na && anyNA(x)) {
    .abort(sprintf("`%s` must not contain missing values.", name))
  }
  finite <- x[!is.na(x)]
  if (any(!is.finite(finite))) {
    .abort(sprintf("`%s` must contain finite values.", name))
  }
  if (any(finite < min | finite > max)) {
    .abort(sprintf("`%s` must lie between %s and %s.", name, format(min), format(max)))
  }
  invisible(TRUE)
}

.assert_positive_integer_vector <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0L || anyNA(x) ||
      any(x != round(x)) || any(x <= 0)) {
    .abort(sprintf("`%s` must be a non-empty vector of positive integers.", name))
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.format_number <- function(x, digits = 4L) {
  if (length(x) == 0L || is.null(x)) {
    return(NA_character_)
  }
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}
