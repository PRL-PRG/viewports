viewports_set_debug_mode <- function(debug=TRUE) {
  if (typeof(debug) != "logical") {
    stop(paste0("Argument is of type", typeof(debug),
                ", must be a logical vector."))
  }
  if (length(debug) == 0) {
    stop("Argument is a 0-length vector")
  }
  if (length(debug) > 1) {
    warning(paste0("Argument is a vector containing multiple values, picking the",
                   "first one, ignoring the rest"))
  }
  invisible(.Call("set_debug_mode", debug))
}

slice <- function(vector, start, size) {
  .Call("create_slice",
        vector,
        .expect_in_range(.expect_exactly_one(.expect_types(start, c("integer", "double"))), 1, length(vector)),
        .expect_in_range(.expect_exactly_one(.expect_types(size,  c("integer", "double"))), 0, length(vector) - start))
}

mosaic <- function(vector, indices_or_mask) {
  if (typeof(indices_or_mask) == "logical") .expect_same_length(indices_or_mask, vector)
  .Call("create_mosaic",
        vector,
        .expect_types(indices_or_mask, c("integer", "double", "logical")))
}

prism <- function(vector, indices) {
  .Call("create_prism",
        vector,
        .expect_types(indices, c("integer", "double")))
}

.expect_exactly_one <- function(vector, name=substitute(vector)) {
  if (length(vector) > 1) {
    warning(paste0("`", name, "` ",
                   "is a vector containing multiple values, ",
                   "picking the first one, ignoring the rest"))
  }
  if (length(vector) == 0) {
    stop(paste0("`", name, "` ",
                "is a zero-length vector, ",
                "but it should be a vector containing a single value"))
  }
  vector
}

.expect_types <- function(vector, expected_types, name=substitute(vector)) {
  if (!(typeof(vector) %in% expected_types)) {
    stop(paste0("`", name, "` ",
                "is a vector of type `", typeof(vector), "` ",
                "but a vector of one of the following types was expected: ", paste(expected_types, collapse=" ")))
  }
  vector
}

.expect_in_range <- function(value, lower_bound, upper_bound, name=substitute(value)) {
    if (!(value >= lower_bound)) {
        stop(paste0("`", name, "`", " should be greater or equal to ", lower_bound))
    }
    if (!(value <= upper_bound)) {
        stop(paste0("`", name, "`", " should be less or equal to ", upper_bound))
    }
    value
}

.expect_same_length <- function(vector, other_vector, name=substitute(vector), other_name=substitute(vector)) {
    if (length(vector) != length(other_vector)) {
        stop(paste0("`", name, "`", " should be the same length as `", other_name, "`"))
    }
}