#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Altrep.h>

#include "debug.h"
#include "helpers.h"

#include "slices.h"
#include "common.h"
#include "mosaics.h"
#include "prisms.h"

#define MAKE_SURE
#include "make_sure.h"

#ifndef NA_RAW
#define NA_RAW 0
#endif

static R_altrep_class_t slice_integer_altrep;
static R_altrep_class_t slice_numeric_altrep;
static R_altrep_class_t slice_logical_altrep;
static R_altrep_class_t slice_complex_altrep;
static R_altrep_class_t slice_raw_altrep;
//static R_altrep_class_t slice_vector_altrep;

static inline R_altrep_class_t class_from_sexp_type(SEXPTYPE type) { // @suppress("No return")
    switch (type) {
        case INTSXP:  return slice_integer_altrep;
        case REALSXP: return slice_numeric_altrep;
        case LGLSXP:  return slice_logical_altrep;
        case CPLXSXP: return slice_complex_altrep;
        case RAWSXP:  return slice_raw_altrep;
        default:      Rf_error("No ALTREP slice class for vector of type %s", type2str(type));
    }
}

# define how_many_ints_in_R_xlen_t (sizeof(R_xlen_t) / sizeof(int))
typedef union {
    int      integers[how_many_ints_in_R_xlen_t];
    R_xlen_t length;
} converter_t;

void read_start_and_size(SEXP/*INTSXP*/ window, R_xlen_t *start, R_xlen_t *size) {
    converter_t start_converter;
    for (int i = 0; i < how_many_ints_in_R_xlen_t; i++) {
        start_converter.integers[i] = INTEGER_ELT(window, i);
    }

    converter_t size_converter;
    for (int i = 0; i < how_many_ints_in_R_xlen_t; i++) {
        size_converter.integers[i] = INTEGER_ELT(window, how_many_ints_in_R_xlen_t + i);
    }

    *start = start_converter.length;
    *size = size_converter.length;
}

SEXP slice_new(SEXP source, R_xlen_t start, R_xlen_t size) {
    make_sure(TYPEOF(source) == INTSXP || TYPEOF(source) == REALSXP || TYPEOF(source) == CPLXSXP
           || TYPEOF(source) == LGLSXP || TYPEOF(source) == VECSXP  || TYPEOF(source) == STRSXP, Rf_error,
		      "type of source should be one of INTSXP, REALSXP, CPLXSXP, LGLSXP, VECSXP, or STRSXP");

    make_sure(start < XLENGTH(source), Rf_error, "start cannot be greater than length of source");
    make_sure(start + size <= XLENGTH(source), Rf_error, "viewport must fit within the length of source");

    if (get_debug_mode()) {
        Rprintf("slice_new\n");
        Rprintf("           SEXP: %p\n", source);
        Rprintf("          start: %li\n", start);
        Rprintf("           size: %li\n", size);
    }

    SEXP/*INTSXP*/ window = allocVector(INTSXP, 2 * how_many_ints_in_R_xlen_t);

    converter_t start_converter = { .length = start };
    converter_t size_converter  = { .length = size  };

    for (size_t i = 0; i < how_many_ints_in_R_xlen_t; i++) {
        SET_INTEGER_ELT(window, i, start_converter.integers[i]);
    }

    for (size_t i = 0; i < how_many_ints_in_R_xlen_t; i++) {
        SET_INTEGER_ELT(window, how_many_ints_in_R_xlen_t + i, size_converter.integers[i]);
    }

    SEXP/*LISTSXP*/ data = allocSExp(LISTSXP);
    SETCAR (data, source);     // The original vector
    SET_TAG(data, R_NilValue); // Starts as R_NilValue, becomes a vector if it the slice is written to
    SETCDR (data, R_NilValue); // Nothing here

    return R_new_altrep(class_from_sexp_type(TYPEOF(source)), window, data);
}

static inline SEXP/*INTSXP*/ get_window(SEXP x) {
    return R_altrep_data1(x);
}

static inline SEXP get_source(SEXP x) {
    SEXP/*LISTSXP*/ cell =  R_altrep_data2(x);
    return CAR(cell);
}

static inline SEXP get_materialized_data(SEXP x) {
    SEXP/*LISTSXP*/ cell =  R_altrep_data2(x);
    return TAG(cell);
}

static inline void set_materialized_data(SEXP x, SEXP data) {
    SEXP/*LISTSXP*/ cell =  R_altrep_data2(x);
    SET_TAG(cell, data);
}

static inline bool is_materialized(SEXP x) {
    SEXP/*LISTSXP*/ cell =  R_altrep_data2(x);
    return TAG(cell) != R_NilValue;
}

SEXP slice_duplicate(SEXP x, Rboolean deep) {

    if (get_debug_mode()) {
        Rprintf("slice_duplicate\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("           deep: %i\n", deep);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*INTSXP*/ window = get_window(x);
    SEXP           source = get_source(x);

    if (deep) {
        R_xlen_t start = 0;
        R_xlen_t size  = 0;
        read_start_and_size(window, &start, &size);
        SEXP slice = slice_new(source, start, size);
        if (is_materialized(x)) {
            SEXP data = get_materialized_data(x);
            set_materialized_data(slice, duplicate(data));
        }
        return slice;
    } else {
        SEXP/*LISTSXP*/ meta = allocSExp(LISTSXP);
        SETCAR (meta, source);     // The original vector
        SEXP data = is_materialized(x) ? get_materialized_data(x) : R_NilValue;
        SET_TAG(meta, data);       // Starts as R_NilValue, becomes a vector if it the slice is written to
        SETCDR (meta, R_NilValue); // Nothing here
        return R_new_altrep(class_from_sexp_type(TYPEOF(x)), window, meta);
    }
}

static Rboolean slice_inspect(SEXP x, int pre, int deep, int pvec, void (*inspect_subtree)(SEXP, int, int, int)) {

    Rprintf("slice_altrep %s\n", type2char(TYPEOF(x)));

    inspect_subtree(R_altrep_data1(x), pre, deep, pvec);
    inspect_subtree(R_altrep_data2(x), pre, deep, pvec);

    return FALSE;
}

static R_xlen_t slice_length(SEXP x) {
    if (get_debug_mode()) {
        Rprintf("slice_length\n");
        Rprintf("           SEXP: %p\n", x);
    }

    SEXP/*INTSXP*/ window = get_window(x);

    R_xlen_t start = 0;
    R_xlen_t size  = 0;
    read_start_and_size(window, &start, &size);

    return size;
}

const void *extract_read_only_data_pointer(SEXP x) {
    SEXP/*INTSXP*/ window = get_window(x);
    SEXP           source = get_source(x);

    R_xlen_t start = 0;
    R_xlen_t size  = 0;
    read_start_and_size(window, &start, &size);

    const void *data = DATAPTR_RO(source);

    SEXPTYPE type = TYPEOF(source);
    switch (type) {
        case INTSXP:  {
            const int *ints = (const int *) data;
            return (const void *) (&ints[start]);
        }
        case REALSXP: {
            const double *doubles = (const double *) data;
            return (const void *) (&doubles[start]);
        }
        case LGLSXP: {
            const Rboolean *booleans = (const Rboolean *) data;
            return (const void *) (&booleans[start]);
        }
        case RAWSXP: {
            const Rbyte *bytes = (const Rbyte *) data;
            return (const void *) (&bytes[start]);
        }
        case CPLXSXP: {
            const Rcomplex *trouble = (const Rcomplex *) data;
            return (const void *) (&trouble[start]);
        }
        case STRSXP:
        case VECSXP: {
            const SEXP *sexps = (const SEXP *) data;
            return (const void *) (&sexps[start]);
        }
        default: Rf_error("Illegal source type for a slice: %i.\n", type);
    }

    make_sure(false, Rf_error, "unreachable");
    return NULL;
}

static void *slice_dataptr(SEXP x, Rboolean writeable) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    if (get_debug_mode()) {
        Rprintf("slice_dataptr\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("      writeable: %i\n", writeable);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP data = get_materialized_data(x);
        return (writeable) ? DATAPTR(data) : ((void *) DATAPTR_RO(data));
    }

    if (writeable) {
        SEXP/*INTSXP*/ window = get_window(x);
        SEXP           source = get_source(x);

        R_xlen_t start = 0;
        R_xlen_t size  = 0;
        read_start_and_size(window, &start, &size);

        SEXP data = copy_data_in_range(source, start, size);
        set_materialized_data(x, data);
        return DATAPTR(data);
    }

    return (void *) extract_read_only_data_pointer(x);
}

static const void *slice_dataptr_or_null(SEXP x) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    if (get_debug_mode()) {
        Rprintf("slice_dataptr_or_null\n");
        Rprintf("           SEXP: %p\n", x);
    }

    return extract_read_only_data_pointer(x);
}

R_xlen_t project_index(SEXP/*INTSXP*/ window, R_xlen_t index) {
	make_sure(window != NULL && window != R_NilValue, Rf_error, "window cannot be null");

    R_xlen_t start = 0;
    R_xlen_t size  = 0;
    read_start_and_size(window, &start, &size);

    make_sure(index < size, Rf_error, "index out of range");
    R_xlen_t projected_index = start + index;

    if (get_debug_mode()) {
        Rprintf("    projecting index\n");
        Rprintf("         window SEXP: %p\n",  window);
        Rprintf("               start: %li\n", start);
        Rprintf("                size: %li\n", size);
        Rprintf("         input index: %li\n", index);
        Rprintf("     projected index: %li\n", projected_index);
    }

    return projected_index;
}

static inline bool index_within_bounds(SEXP window, R_xlen_t index) {
    R_xlen_t start = 0;
    R_xlen_t size  = 0;
    read_start_and_size(window, &start, &size);
    return index < size;
}

static int slice_integer_element(SEXP x, R_xlen_t i) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    SEXP/*INTSXP*/ window = get_window(x);
    SEXP/*INTSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == INTSXP, Rf_error, "type of source must be INTSXP");

    if (get_debug_mode()) {
        Rprintf("slice_integer_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          source: %p\n", source);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*INTSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == INTSXP, Rf_error, "type of data must be INTSXP");
        return INTEGER_ELT(data, i);
    }

    if (!index_within_bounds(window, i)) {
        return NA_INTEGER;
    }

    R_xlen_t projected_index = project_index(window, i);

    return INTEGER_ELT(source, projected_index);
}

static double slice_numeric_element(SEXP x, R_xlen_t i) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    SEXP/*INTSXP*/  window = get_window(x);
    SEXP/*CPLXSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == REALSXP, Rf_error, "type of source must be REALSXP");

    if (get_debug_mode()) {
        Rprintf("slice_numeric_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          source: %p\n", source);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*REALSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == REALSXP, Rf_error, "type of data must be REALSXP");
        return REAL_ELT(data, i);
    }

    if (!index_within_bounds(window, i)) {
        return NA_REAL;
    }

    R_xlen_t projected_index = project_index(window, i);

    return REAL_ELT(source, projected_index);
}

static Rbyte slice_raw_element(SEXP x, R_xlen_t i) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    SEXP/*INTSXP*/  window = get_window(x);
    SEXP/*CPLXSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == RAWSXP, Rf_error, "type of source must be RAWSXP");

    if (get_debug_mode()) {
        Rprintf("slice_raw_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          source: %p\n", source);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*RAWSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == RAWSXP, Rf_error, "type of data must be RAWSXP");
        return RAW_ELT(data, i);
    }

    if (!index_within_bounds(window, i)) {
        return NA_RAW;
    }

    R_xlen_t projected_index = project_index(window, i);

    return RAW_ELT(source, projected_index);
}

static Rcomplex slice_complex_element(SEXP x, R_xlen_t i) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    SEXP/*INTSXP*/  window = get_window(x);
    SEXP/*CPLXSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == CPLXSXP, Rf_error, "type of source must be CPLXSXP");

    if (get_debug_mode()) {
        Rprintf("slice_complex_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          source: %p\n", source);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*CPLXSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == CPLXSXP, Rf_error, "type of data must be CPLXSXP");
        return COMPLEX_ELT(data, i);
    }

    if (!index_within_bounds(window, i)) {
#ifndef NA_COMPLEX
        Rcomplex NA_COMPLEX;
        NA_COMPLEX.i = NA_REAL;
        NA_COMPLEX.r = NA_REAL;
#endif
        return NA_COMPLEX;
    }

    R_xlen_t projected_index = project_index(window, i);

    return COMPLEX_ELT(source, projected_index);
}

static int slice_logical_element(SEXP x, R_xlen_t i) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    SEXP/*INTSXP*/ window = get_window(x);
    SEXP/*LGLSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == LGLSXP, Rf_error, "type of source must be LGLSXP");

    if (get_debug_mode()) {
        Rprintf("slice_logical_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          source: %p\n", source);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*LGLSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == LGLSXP, Rf_error, "type of data must be LGLSXP");
        return LOGICAL_ELT(data, i);
    }

    if (!index_within_bounds(window, i)) {
        return NA_LOGICAL;
    }

    R_xlen_t projected_index = project_index(window, i);

    return LOGICAL_ELT(source, projected_index);
}

static R_xlen_t slice_integer_get_region(SEXP x, R_xlen_t i, R_xlen_t n, int *buf) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    if (get_debug_mode()) {
        Rprintf("slice_integer_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*INTSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == INTSXP, Rf_error, "type of data must be INTSXP");
        return INTEGER_GET_REGION(data, i, n, buf);
    }

    SEXP/*INTSXP*/ window = get_window(x);
    SEXP/*INTSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == INTSXP, Rf_error, "type of source must be INTSXP");

    R_xlen_t projected_index = project_index(window, i);

    return INTEGER_GET_REGION(source, projected_index, n, buf);
}

static R_xlen_t slice_numeric_get_region(SEXP x, R_xlen_t i, R_xlen_t n, double *buf) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    if (get_debug_mode()) {
        Rprintf("slice_numeric_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*REALSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == REALSXP, Rf_error, "type of data must be REALSXP");
        return REAL_GET_REGION(data, i, n, buf);
    }

    SEXP/*INTSXP*/ window = get_window(x);
    SEXP/*REALSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == REALSXP, Rf_error, "type of source must be REALSXP");

    R_xlen_t projected_index = project_index(window, i);

    return REAL_GET_REGION(source, projected_index, n, buf);
}

static R_xlen_t slice_raw_get_region(SEXP x, R_xlen_t i, R_xlen_t n, Rbyte *buf) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    if (get_debug_mode()) {
        Rprintf("slice_raw_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*RAWSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == RAWSXP, Rf_error, "type of data must be RAWSXP");
        return RAW_GET_REGION(data, i, n, buf);
    }

    SEXP/*INTSXP*/ window = get_window(x);
    SEXP/*RAWSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == RAWSXP, Rf_error, "type of source must be RAWSXP");

    R_xlen_t projected_index = project_index(window, i);

    return RAW_GET_REGION(source, projected_index, n, buf);
}

static R_xlen_t slice_complex_get_region(SEXP x, R_xlen_t i, R_xlen_t n, Rcomplex *buf) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    if (get_debug_mode()) {
        Rprintf("slice_complex_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*CPLXSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == CPLXSXP, Rf_error, "type of data must be CPLXSXP");
        return COMPLEX_GET_REGION(data, i, n, buf);
    }

    SEXP/*INTSXP*/  window = get_window(x);
    SEXP/*CPLXSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == CPLXSXP, Rf_error, "type of source must be CPLXSXP");

    R_xlen_t projected_index = project_index(window, i);

    return COMPLEX_GET_REGION(source, projected_index, n, buf);
}

static R_xlen_t slice_logical_get_region(SEXP x, R_xlen_t i, R_xlen_t n, int *buf) {
	make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    if (get_debug_mode()) {
        Rprintf("slice_logical_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*LGLSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == LGLSXP, Rf_error, "type of data must be LGLSXP");
        return LOGICAL_GET_REGION(data, i, n, buf);
    }

    SEXP/*INTSXP*/ window = get_window(x);
    SEXP/*LGLSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == LGLSXP, Rf_error, "type of source must be LGLSXP");

    R_xlen_t projected_index = project_index(window, i);

    return LOGICAL_GET_REGION(source, projected_index, n, buf);
}

SEXP translate_indices(SEXP original, R_xlen_t offset, R_xlen_t size) {
	make_sure(TYPEOF(original) == INTSXP || TYPEOF(original) == REALSXP, Rf_error,
			  "type of original must be either INTSXP or REALSXP");
	make_sure(sizeof(R_xlen_t) <= sizeof(double), Rf_error,
			  "a vector of doubles must be able to contain elements of type R_xlen_t");

	SEXP translated = allocVector(REALSXP, XLENGTH(original)); // FIXME translation altrep vector

	switch (TYPEOF(original)) {
	case INTSXP:
		for (R_xlen_t i = 0; i < XLENGTH(translated); i++) {
			int original_element = INTEGER_ELT(original, i);
			if (original_element == NA_INTEGER || ((R_xlen_t) original_element) > size) {
				SET_REAL_ELT(translated, i, NA_REAL);
			} else {
				SET_REAL_ELT(translated, i, ((R_xlen_t) original_element) + offset);
			}
		}
		break;

	case REALSXP:
		for (R_xlen_t i = 0; i < XLENGTH(translated); i++) {
			double original_element = REAL_ELT(original, i);
			if (ISNAN(original_element) || original_element > size) {
				SET_REAL_ELT(translated, i, NA_REAL);
			} else {
				SET_REAL_ELT(translated, i, ((R_xlen_t) original_element) + offset);
			}
		}
		break;

	default:
		Rf_error("Indices are expected to be either INTSXP or REALSXP");
	}

	return translated;
}

static SEXP slice_extract_subset(SEXP x, SEXP indices, SEXP call) {
    make_sure(x != NULL && x != R_NilValue, Rf_error, "x cannot be null");

    if (get_debug_mode()) {
        Rprintf("slice_extract_subset\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("        indices: %p\n", indices);
        Rprintf("           call: %p\n", call);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*INTSXP*/ window = get_window(x);
    SEXP           source = get_source(x);

    // No indices.
    R_xlen_t size = XLENGTH(indices);
    if (size == 0) {
        return allocVector(TYPEOF(source), 0);
    }

    R_xlen_t window_start = 0;
    R_xlen_t window_size  = 0;
    read_start_and_size(window, &window_start, &window_size);

    if (is_materialized(x)) {
            return copy_data_at_indices(source, indices);
        }

    if (!are_indices_in_range(indices, 1, window_size)) {
    	SEXP translated_indices = translate_indices(indices, window_start, window_size);
    	return copy_data_at_indices(source, translated_indices);
    }

    if (!are_indices_contiguous(indices)) {
    	SEXP translated_indices = translate_indices(indices, window_start, window_size);
        if (are_indices_monotonic(indices)) {
        	return create_mosaic(source, translated_indices);
        } else {
        	return create_prism(source, translated_indices);
        }
    }

    // Contiguous indices.
    R_xlen_t start = get_first_element_as_length(indices) - 1;
    R_xlen_t projected_start = project_index(window, start);
    return slice_new(source, projected_start, size);
}

// R_set_altstring_Set_elt_method
// static void string_set_elt(SEXP x, R_xlen_t i, SEXP v)

// R_set_altstring_Elt_method
// string_elt(SEXP x, R_xlen_t i)

void init_common(R_altrep_class_t cls) {
    R_set_altrep_Duplicate_method(cls, slice_duplicate);
    R_set_altrep_Inspect_method(cls, slice_inspect);
    R_set_altrep_Length_method(cls, slice_length);

    R_set_altvec_Dataptr_method(cls, slice_dataptr);
    R_set_altvec_Dataptr_or_null_method(cls, slice_dataptr_or_null);
    R_set_altvec_Extract_subset_method(cls, slice_extract_subset);
}

// UFO Inits
void init_slice_integer_altrep_class(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altinteger_class("slice_integer_altrep", "viewports", dll);
    slice_integer_altrep = cls;

    init_common(cls);

    R_set_altinteger_Elt_method(cls, slice_integer_element);
    R_set_altinteger_Get_region_method(cls, slice_integer_get_region);
}

void init_slice_logical_altrep_class(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altlogical_class("slice_logical_altrep", "viewports", dll);
    slice_logical_altrep = cls;

    init_common(cls);

    R_set_altlogical_Elt_method(cls, slice_logical_element);
    R_set_altlogical_Get_region_method(cls, slice_logical_get_region);
}

void init_slice_numeric_altrep_class(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altreal_class("slice_numeric_altrep", "viewports", dll);
    slice_numeric_altrep = cls;

    init_common(cls);

    R_set_altreal_Elt_method(cls, slice_numeric_element);
    R_set_altreal_Get_region_method(cls, slice_numeric_get_region);
}

void init_slice_complex_altrep_class(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altcomplex_class("slice_complex_altrep", "viewports", dll);
    slice_complex_altrep = cls;

    init_common(cls);

    R_set_altcomplex_Elt_method(cls, slice_complex_element);
    R_set_altcomplex_Get_region_method(cls, slice_complex_get_region);
}

void init_slice_raw_altrep_class(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altraw_class("slice_raw_altrep", "viewports", dll);
    slice_raw_altrep = cls;

    init_common(cls);

    R_set_altraw_Elt_method(cls, slice_raw_element);
    R_set_altraw_Get_region_method(cls, slice_raw_get_region);
}

void init_slice_altrep_class(DllInfo * dll) {
    init_slice_integer_altrep_class(dll);
    init_slice_numeric_altrep_class(dll);
    init_slice_logical_altrep_class(dll);
    init_slice_complex_altrep_class(dll);
    init_slice_raw_altrep_class(dll);
}

SEXP create_slice(SEXP source, SEXP/*INTSXP|REALSXP*/ start_sexp, SEXP/*INTSXP|REALSXP*/ size_sexp) {
    make_sure(TYPEOF(start_sexp) == INTSXP || TYPEOF(start_sexp) == REALSXP, Rf_error,
    		  "type of start must be either INTSXP or REALSXP");
    make_sure(TYPEOF(size_sexp) == INTSXP || TYPEOF(size_sexp) == REALSXP, Rf_error,
    		  "type of size must be either INTSXP or REALSXP");
    make_sure(XLENGTH(start_sexp) > 0, Rf_error, "start cannot be a zero-length vector");
    make_sure(XLENGTH(size_sexp) > 0, Rf_error, "size cannot be a zero-length vector");

    R_xlen_t start       = get_first_element_as_length(start_sexp) - 1;
    R_xlen_t size        = get_first_element_as_length(size_sexp);

    make_sure(start < XLENGTH(source), Rf_error, "start cannot be greater than length of source");
    make_sure(start + size <= XLENGTH(source), Rf_error, "viewport must fit within the length of source");

    if (get_debug_mode()) {
        Rprintf("create slice\n");
        Rprintf("           SEXP: %p\n", source);
        Rprintf("          start: %li\n", start);
        Rprintf("           size: %li\n", size);
    }

    return slice_new(source, start, size);
}
