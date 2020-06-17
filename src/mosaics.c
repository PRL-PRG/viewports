#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Altrep.h>

#include "debug.h"
#include "helpers.h"
#include "bitmap_sexp.h"

#include "common.h"
#include "mosaics.h"

#define MAKE_SURE
#include "make_sure.h"

static R_altrep_class_t mosaic_integer_altrep;
static R_altrep_class_t mosaic_numeric_altrep;
static R_altrep_class_t mosaic_logical_altrep;
static R_altrep_class_t mosaic_complex_altrep;
static R_altrep_class_t mosaic_raw_altrep;

static inline R_altrep_class_t class_from_sexp_type(SEXPTYPE type) { // @suppress("No return")
    switch (type) {
        case INTSXP:  return mosaic_integer_altrep;
        case REALSXP: return mosaic_numeric_altrep;
        case LGLSXP:  return mosaic_logical_altrep;
        case CPLXSXP: return mosaic_complex_altrep;
        case RAWSXP:  return mosaic_raw_altrep;
        default:      Rf_error("No ALTREP mosaic class for vector of type %s", type2str(type));
    }
}

#define how_many_ints_in_R_xlen_t (sizeof(R_xlen_t) / sizeof(int))

R_xlen_t convert_logical_mask_to_bitmap(SEXP/*LGLSXP*/ mask, SEXP/*INTSXP*/ bitmap) {
    make_sure(TYPEOF(mask) == LGLSXP, Rf_error, "mask must be a vector of type LGLSXP");

    R_xlen_t size = XLENGTH(mask);
    make_sure(XLENGTH(mask) == XTRUELENGTH(bitmap), Rf_error, "mask must the same length as the bitmap");

    R_xlen_t elements = 0;
    for (R_xlen_t i = 0; i < size; i++) {
        Rboolean current = LOGICAL_ELT(mask, i);
        if (current == NA_LOGICAL) {
            Rf_error("Mosaics cannot be created from a logical mask containing NA\n");
        }
        if (current == TRUE) {
            bitmap_set(bitmap, i);
            elements++;
        }
    }

    return elements;
}

R_xlen_t convert_integer_indices_to_bitmap(SEXP/*INTSXP*/ indices, SEXP/*INTSXP*/ bitmap) {
    make_sure(TYPEOF(indices) == INTSXP, Rf_error, "type of indices must be INTSXP");
    R_xlen_t size = XLENGTH(indices);

    int previous = NA_INTEGER;
    for (int i = 0; i < size; i++) {
        int current = INTEGER_ELT(indices, i);
        if (current == NA_INTEGER) {
            Rf_error("Mosaics cannot be created from an ordered index containing NA\n");
        }
        if (previous != NA_INTEGER && previous >= current) {
            Rf_error("Mosaics can only be created from an ordered index list, but %li >= %li\n", previous, current);
        }
        make_sure(current > 0, Rf_error, "index value must be greater than zero");
        bitmap_set(bitmap, (R_xlen_t) current - 1);
        previous = current;
    }

    return size;
}

R_xlen_t convert_numeric_indices_to_bitmap(SEXP/*REALSXP*/ indices, SEXP/*INTSXP*/ bitmap){
	make_sure(TYPEOF(indices) == REALSXP, Rf_error, "type of indices must be REALSXP");
    R_xlen_t size = XLENGTH(indices);

    double previous = NA_REAL;
    for (int i = 0; i < size; i++) {
        double current = REAL_ELT(indices, i);
        if (ISNAN(current)) {
            Rf_error("Mosaics cannot be created from an ordered index containing NA\n");
        }
        if (!ISNAN(previous) && previous >= current) {
            Rf_error("Mosaics can only be created from an ordered index list, but %li >= %li\n", previous, current);
        }
        make_sure(current > 0, Rf_error, "index value must be greater than zero");
        bitmap_set(bitmap, (R_xlen_t) current - 1);
        previous = current;
    }

    return size;
}


R_xlen_t convert_indices_to_bitmap(SEXP/*INTSXP | REALSXP | LGLSXP*/ indices, SEXP/*INTSXP*/ bitmap) {  // @suppress("No return")
    SEXPTYPE type = TYPEOF(indices);
    make_sure(type == INTSXP || type == REALSXP || type == LGLSXP, Rf_error,
    		  "type of indices must be one of INTSXP, REALSXP, or LGLSXP");

    switch (type) {
        case INTSXP:  return convert_integer_indices_to_bitmap(indices, bitmap);
        case REALSXP: return convert_numeric_indices_to_bitmap(indices, bitmap);
        case LGLSXP:  return convert_logical_mask_to_bitmap   (indices, bitmap);
        default:      Rf_error("Mosaics can be indexed by logical, integer, or numeric vectors but found: %d\n", type);
    }
}

SEXP/*A*/ mosaic_new(SEXP/*A*/ source, SEXP/*INTSXP*/ bitmap, R_xlen_t size) {
    make_sure(TYPEOF(source) == INTSXP
           || TYPEOF(source) == REALSXP
		   || TYPEOF(source) == RAWSXP
           || TYPEOF(source) == CPLXSXP
           || TYPEOF(source) == LGLSXP
           || TYPEOF(source) == VECSXP
           || TYPEOF(source) == STRSXP, Rf_error,
		      "type of indices must be one of INTSXP, REALSXP, RAWSXP, LGLSXP, CPLXSXP, VECSXP, or STRSXP");

    if (get_debug_mode()) {
        Rprintf("mosaic_new\n");
        Rprintf("           SEXP: %p\n", source);
        Rprintf("         bitmap: %p\n",  bitmap);
        Rprintf("           size: %li\n", size);
    }

    SEXP/*REALSXP*/ length_vector = allocVector(REALSXP, 1);
    SET_REAL_ELT(length_vector, 0, size);

    SEXP/*LISTSXP*/ data = allocSExp(LISTSXP);
    SETCAR (data, source);        // The original vector
    SET_TAG(data, R_NilValue);    // Starts as R_NilValue, becomes a vector if it the mosaic is written to
    SETCDR (data, length_vector); // The number of set bits in the mask goes here, aka length

    return R_new_altrep(class_from_sexp_type(TYPEOF(source)), bitmap, data);
}

static inline SEXP/*INTSXP*/ get_bitmap(SEXP x) {
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

static inline SEXP get_length_vector(SEXP x) {
    SEXP/*LISTSXP*/ cell =  R_altrep_data2(x);
    return CDR(cell);
}

static inline R_xlen_t get_length(SEXP x) {
    SEXP/*LISTSXP*/ cell =  R_altrep_data2(x);
    return REAL_ELT(CDR(cell), 0);
}

SEXP mosaic_duplicate(SEXP x, Rboolean deep) {

    if (get_debug_mode()) {
        Rprintf("mosaic_duplicate\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("           deep: %i\n", deep);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*INTSXP*/ bitmap = get_bitmap(x);
    SEXP           source = get_source(x);
    R_xlen_t       length = get_length(x);

    if (deep) {
        SEXP mosaic = mosaic_new(source, bitmap, length);
        if (is_materialized(x)) {
            SEXP data = get_materialized_data(x);
            set_materialized_data(mosaic, duplicate(data));
        }
        return mosaic;
    } else {
        SEXP/*LISTSXP*/ meta = allocSExp(LISTSXP);
        SETCAR (meta, source);               // The original vector
        SEXP data = is_materialized(x) ? get_materialized_data(x) : R_NilValue;
        SET_TAG(meta, data);                 // Starts as R_NilValue, becomes a vector if it the mosaic is written to
        SETCDR (meta, get_length_vector(x)); // Length vector here
        return R_new_altrep(class_from_sexp_type(TYPEOF(source)), bitmap, meta);
    }
}

static Rboolean mosaic_inspect(SEXP x, int pre, int deep, int pvec, void (*inspect_subtree)(SEXP, int, int, int)) {

    Rprintf("mosaic_altrep %s\n", type2char(TYPEOF(x)));

    inspect_subtree(R_altrep_data1(x), pre, deep, pvec);
    inspect_subtree(R_altrep_data2(x), pre, deep, pvec);

    return FALSE;
}

static R_xlen_t mosaic_length(SEXP x) {
    if (get_debug_mode()) {
        Rprintf("mosaic_length\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("         length: %p\n", get_length(x));
    }
    return get_length(x);
}

SEXP copy_from_source(SEXP x) {
    SEXP/*INTSXP*/ bitmap = get_bitmap(x);
    SEXP           source = get_source(x);
    R_xlen_t       length = get_length(x);

    make_sure(XTRUELENGTH(bitmap) == XLENGTH(source), Rf_error,
    		  "bitmap must be the same length as source");

    SEXP materialized = allocVector(TYPEOF(source), length);
    R_xlen_t cursor = 0;
    for (R_xlen_t index = 0; index < XLENGTH(source); index++) {
        if (bitmap_get(bitmap, index)) {                            // FIXME Conditional jump or move depends on uninitialised value(s)
            copy_element(source, index, materialized, cursor);
            cursor++;
        }
    }
    make_sure(cursor == length, Rf_error,
    		  "the number of copied elements is different than the length of the output vector");
    return materialized;
}

static void *mosaic_dataptr(SEXP x, Rboolean writeable) {
    make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_dataptr\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("      writeable: %i\n", writeable);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP data = get_materialized_data(x);
        return writeable ? DATAPTR(data) : (void *) DATAPTR_RO(data);
    }

    SEXP data = copy_from_source(x);
    set_materialized_data(x, data);
    return writeable ? DATAPTR(data) : (void *) DATAPTR_RO(data);
}

static const void *mosaic_dataptr_or_null(SEXP x) {
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_dataptr_or_null\n");
        Rprintf("           SEXP: %p\n", x);
    }

    if (is_materialized(x)) {
        SEXP data = get_materialized_data(x);
        return DATAPTR_RO(data);
    }

    SEXP data = copy_from_source(x);
    set_materialized_data(x, data);
    return DATAPTR_RO(data);
}

static int mosaic_integer_element(SEXP x, R_xlen_t i) {  // CONTINUE HERE
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_integer_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*INTSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == INTSXP, Rf_error, "type of data must be INTSXP");
        return INTEGER_ELT(data, i);
    }

    SEXP/*INTSXP*/ bitmap = get_bitmap(x);
    SEXP/*INTSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == INTSXP, Rf_error, "type of source must be INTSXP");
    R_xlen_t projected_index = bitmap_index_of_nth_set_bit(bitmap, i);
    return INTEGER_ELT(source, projected_index);
}

static double mosaic_numeric_element(SEXP x, R_xlen_t i) {
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_numeric_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*REALSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == REALSXP, Rf_error, "type of data must be REALSXP");
        return REAL_ELT(data, i);
    }

    SEXP/*INTSXP*/  bitmap = get_bitmap(x);
    SEXP/*REALSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == REALSXP, Rf_error, "type of source must be REALSXP");
    R_xlen_t projected_index = bitmap_index_of_nth_set_bit(bitmap, i);
    return REAL_ELT(source, projected_index);
}

static Rbyte mosaic_raw_element(SEXP x, R_xlen_t i) {
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_raw_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*RAWSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == RAWSXP, Rf_error, "type of data must be RAWSXP");
        return RAW_ELT(data, i);
    }

    SEXP/*INTSXP*/ bitmap = get_bitmap(x);
    SEXP/*RAWSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == RAWSXP, Rf_error, "type of source must be RAWSXP");
    R_xlen_t projected_index = bitmap_index_of_nth_set_bit(bitmap, i);
    return RAW_ELT(source, projected_index);
}

static Rcomplex mosaic_complex_element(SEXP x, R_xlen_t i) {
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_complex_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*CPLXSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == CPLXSXP, Rf_error, "type of data must be CPLXSXP");
        return COMPLEX_ELT(data, i);
    }

    SEXP/*INTSXP*/  bitmap = get_bitmap(x);
    SEXP/*CPLXSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == CPLXSXP, Rf_error, "type of source must be CPLXSXP");
    R_xlen_t projected_index = bitmap_index_of_nth_set_bit(bitmap, i);
    return COMPLEX_ELT(source, projected_index);
}

static int mosaic_logical_element(SEXP x, R_xlen_t i) {
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_logical_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*LGLSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == LGLSXP, Rf_error, "type of data must be LGLSXP");
        return LOGICAL_ELT(data, i);
    }

    SEXP/*INTSXP*/ bitmap = get_bitmap(x);
    SEXP/*LGLSXP*/ source = get_source(x);

    make_sure(TYPEOF(source) == LGLSXP, Rf_error, "type of source must be LGLSXP");
    R_xlen_t projected_index = bitmap_index_of_nth_set_bit(bitmap, i);
    return LOGICAL_ELT(source, projected_index);
}

static R_xlen_t mosaic_integer_get_region(SEXP x, R_xlen_t i, R_xlen_t n, int *buf) {
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_integer_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*INTSXP*/ data;
    if (!is_materialized(x)) {
        data = copy_from_source(x);
        set_materialized_data(x, data);
    } else {
        data = get_materialized_data(x);
    }

    make_sure(TYPEOF(data) == INTSXP, Rf_error, "type of data must be INTSXP");
    return INTEGER_GET_REGION(data, i, n, buf);
}

static R_xlen_t mosaic_numeric_get_region(SEXP x, R_xlen_t i, R_xlen_t n, double *buf) {
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_numeric_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*REALSXP*/ data;
    if (!is_materialized(x)) {
        data = copy_from_source(x);
        set_materialized_data(x, data);
    } else {
        data = get_materialized_data(x);
    }

    make_sure(TYPEOF(data) == REALSXP, Rf_error, "type of data must be REALSXP");
    return REAL_GET_REGION(data, i, n, buf);
}

static R_xlen_t mosaic_raw_get_region(SEXP x, R_xlen_t i, R_xlen_t n, Rbyte *buf) {
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_raw_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*RAWSXP*/ data;
    if (!is_materialized(x)) {
        data = copy_from_source(x);
        set_materialized_data(x, data);
    } else {
        data = get_materialized_data(x);
    }

    make_sure(TYPEOF(data) == RAWSXP, Rf_error, "type of data must be RAWSXP");
    return RAW_GET_REGION(data, i, n, buf);
}

static R_xlen_t mosaic_complex_get_region(SEXP x, R_xlen_t i, R_xlen_t n, Rcomplex *buf) {
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_complex_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*CPLXSXP*/ data;
    if (!is_materialized(x)) {
        data = copy_from_source(x);
        set_materialized_data(x, data);
    } else {
        data = get_materialized_data(x);
    }

    make_sure(TYPEOF(data) == CPLXSXP, Rf_error, "type of data must be CPLXSXP");
    return COMPLEX_GET_REGION(data, i, n, buf);
}

static R_xlen_t mosaic_logical_get_region(SEXP x, R_xlen_t i, R_xlen_t n, int *buf) {
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("mosaic_logical_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*LGLSXP*/ data;
    if (!is_materialized(x)) {
        data = copy_from_source(x);
        set_materialized_data(x, data);
    } else {
        data = get_materialized_data(x);
    }

    make_sure(TYPEOF(data) == LGLSXP, Rf_error, "type of data must be LGLSXP");
    return LOGICAL_GET_REGION(data, i, n, buf);
}

SEXP/*REALSXP*/ translate_indices_by_bitmap(SEXP/*REALSXP*/ screened_indices, SEXP/*bitmap*/ bitmap) {
    make_sure(TYPEOF(screened_indices) == REALSXP, Rf_error, "type of screened_indices must be REALSXP");
    make_sure(TYPEOF(bitmap) == INTSXP, Rf_error, "type of bitmap must be INTSXP");

	SEXP/*REALSXP*/ translated_indices = allocVector(REALSXP, XLENGTH(screened_indices));

	for (R_xlen_t i = 0; i < XLENGTH(screened_indices); i++) {
		double index_but_mistyped = REAL_ELT(screened_indices, i);

		if (ISNAN(index_but_mistyped)) {
			SET_REAL_ELT(translated_indices, i, index_but_mistyped);
			continue;
		}

		R_xlen_t index = (R_xlen_t) index_but_mistyped;
		R_xlen_t translated_index = bitmap_index_of_nth_set_bit(bitmap, index - 1) + 1;
		SET_REAL_ELT(translated_indices, i, translated_index);
	}

	return translated_indices;
}

SEXP/*bitmap*/ translate_bitmap(SEXP source, SEXP/*bitmap*/ bitmap, SEXP/*INTSXP|REALSXP*/ indices) {
	make_sure(TYPEOF(indices) == INTSXP ||TYPEOF(indices) == REALSXP, Rf_error,
			  "type of indices must be either INTSXP or REALSXP");

    R_xlen_t translated_bitmap_size = XTRUELENGTH(bitmap);
    SEXP translated_bitmap = bitmap_new(XLENGTH(source));

    R_xlen_t viewport_index = 0;
    R_xlen_t indices_index = 0;

    for (R_xlen_t i = 0; i < translated_bitmap_size; i++) {
    	if (bitmap_get(bitmap, i)) {
        	R_xlen_t index = (R_xlen_t) TYPEOF(indices) == REALSXP ? REAL_ELT(indices, indices_index)
        														   : INTEGER_ELT(indices, indices_index);

            if (viewport_index == index - 1) {
                bitmap_set(translated_bitmap, i);
                indices_index++;
            }
            viewport_index++;
        }
    }

    return translated_bitmap;
}

static SEXP mosaic_extract_subset(SEXP x, SEXP indices, SEXP call) {
    make_sure(x != NULL, Rf_error, "x cannot be null");
    make_sure(TYPEOF(indices) == INTSXP || TYPEOF(indices) == REALSXP, Rf_error,
    		  "type of indices must be either INTSXP or REALSXP");

    if (get_debug_mode()) {
        Rprintf("mosaic_extract_subset\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("        indices: %p\n", indices);
        Rprintf("           call: %p\n", call);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*INTSXP*/ bitmap = get_bitmap(x);
    SEXP           source = get_source(x);
    R_xlen_t       length = get_length(x);

    // No indices.
    R_xlen_t size = XLENGTH(indices);
    if (size == 0) {
        return allocVector(TYPEOF(source), 0);
    }

    SEXP/*REALSXP*/ screened_indices = screen_indices(indices, length);

    if (is_materialized(x)) {
        // TODO maybe instead just return a viewport into the materialized sexp?
        SEXP materialized_data = get_materialized_data(x);
        return copy_data_at_indices(materialized_data, screened_indices);
    }

    if (!are_indices_monotonic(screened_indices)) {
    	SEXP/*REALSXP*/ translated_indices = translate_indices_by_bitmap(screened_indices, bitmap);
        return copy_data_at_indices(source, translated_indices);
    }

    // Monotonic indices.
    SEXP/*bitmap*/ translated_bitmap = translate_bitmap(source, bitmap, indices);
    return mosaic_new(source, translated_bitmap, XLENGTH(indices));
}

// R_set_altstring_Set_elt_method
// static void string_set_elt(SEXP x, R_xlen_t i, SEXP v)

// R_set_altstring_Elt_method
// string_elt(SEXP x, R_xlen_t i)

void init_common_mosaic(R_altrep_class_t cls) {
    R_set_altrep_Duplicate_method(cls, mosaic_duplicate);
    R_set_altrep_Inspect_method(cls, mosaic_inspect);
    R_set_altrep_Length_method(cls, mosaic_length);

    R_set_altvec_Dataptr_method(cls, mosaic_dataptr);
    R_set_altvec_Dataptr_or_null_method(cls, mosaic_dataptr_or_null);
    R_set_altvec_Extract_subset_method(cls, mosaic_extract_subset);
}

void init_integer_mosaic(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altinteger_class("mosaic_integer_altrep", "viewports", dll);
    mosaic_integer_altrep = cls;

    init_common_mosaic(cls);

    R_set_altinteger_Elt_method(cls, mosaic_integer_element);
    R_set_altinteger_Get_region_method(cls, mosaic_integer_get_region);
}

void init_numeric_mosaic(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altreal_class("mosaic_numeric_altrep", "viewports", dll);
    mosaic_numeric_altrep = cls;

    init_common_mosaic(cls);

    R_set_altreal_Elt_method(cls, mosaic_numeric_element);
    R_set_altreal_Get_region_method(cls, mosaic_numeric_get_region);
}

void init_logical_mosaic(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altlogical_class("mosaic_logical_altrep", "viewports", dll);
    mosaic_logical_altrep = cls;

    init_common_mosaic(cls);

    R_set_altlogical_Elt_method(cls, mosaic_logical_element);
    R_set_altlogical_Get_region_method(cls, mosaic_logical_get_region);
}

void init_complex_mosaic(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altcomplex_class("mosaic_complex_altrep", "viewports", dll);
    mosaic_complex_altrep = cls;

    init_common_mosaic(cls);

    R_set_altcomplex_Elt_method(cls, mosaic_complex_element);
    R_set_altcomplex_Get_region_method(cls, mosaic_complex_get_region);
}

void init_raw_mosaic(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altraw_class("mosaic_raw_altrep", "viewports", dll);
    mosaic_raw_altrep = cls;

    init_common_mosaic(cls);

    R_set_altraw_Elt_method(cls, mosaic_raw_element);
    R_set_altraw_Get_region_method(cls, mosaic_raw_get_region);
}

// UFO Inits
void init_mosaic_altrep_class(DllInfo * dll) {
    init_integer_mosaic(dll);
    init_numeric_mosaic(dll);
    init_logical_mosaic(dll);
    init_complex_mosaic(dll);
    init_raw_mosaic(dll);
}

SEXP/*A*/ create_mosaic(SEXP/*A*/ source, SEXP/*INTSXP|REALSXP|LGLSXP*/ indices) {
	SEXPTYPE source_type = TYPEOF(source);
	SEXPTYPE indices_type = TYPEOF(indices);
	R_xlen_t source_length = XLENGTH(source);
	R_xlen_t indices_length = XLENGTH(indices);

    make_sure(source_type == INTSXP || source_type == REALSXP || source_type == CPLXSXP
           || source_type == LGLSXP || source_type == RAWSXP  || source_type == VECSXP
           || source_type == STRSXP, Rf_error,
		      "type of source must be one of INTSXP, REALSXP, RAWSXP, CPLXSXP, LGLSXP, VECSXP, or STRSXP");

    make_sure(indices_type == INTSXP  || indices_type == REALSXP
           || (indices_type == LGLSXP && (indices_length == source_length)), Rf_error,
        	  "type of indices must be either INTSXP or REALSXP or LGLSXP "
        	  "(if it is LGLSXP, the length of indices must be the same as the source)");

    if (indices_type == INTSXP || indices_type == REALSXP)
    	if (!are_indices_in_range(indices, 0, source_length))
    		Rf_error("Cannot use these indices with this source: out of range");

    if (get_debug_mode()) {
        Rprintf("mosaic_new\n");
        Rprintf("           SEXP: %p\n", source);
        Rprintf("        indices: %p\n", indices);
        Rprintf("   indices type: %li\n", indices_type);
        Rprintf(" indices length: %li\n", indices_length);
    }

    SEXP/*INTSXP*/ bitmap = bitmap_new(source_length);
    R_xlen_t how_many_set_bits = convert_indices_to_bitmap(indices, bitmap);
    return mosaic_new(source, bitmap, how_many_set_bits);
}
