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

#define MAKE_SURE
#include "make_sure.h"

static R_altrep_class_t prism_integer_altrep;
static R_altrep_class_t prism_numeric_altrep;
static R_altrep_class_t prism_logical_altrep;
static R_altrep_class_t prism_complex_altrep;
static R_altrep_class_t prism_raw_altrep;

static inline R_altrep_class_t class_from_sexp_type(SEXPTYPE type) { // @suppress("No return")
    switch (type) {
        case INTSXP:  return prism_integer_altrep;
        case REALSXP: return prism_numeric_altrep;
        case LGLSXP:  return prism_logical_altrep;
        case CPLXSXP: return prism_complex_altrep;
        case RAWSXP:  return prism_raw_altrep;
        default:      Rf_error("No ALTREP prism class for vector of type %s", type2str(type));
    }
}

SEXP prism_new(SEXP source, SEXP/*INTSXP|REALSXP*/ indices) {
    make_sure(sizeof(double) >= sizeof(R_xlen_t), Rf_error, "a vector of doubles cannot serve as a vector of R_xlen_t");
    make_sure(TYPEOF(indices) == REALSXP || TYPEOF(indices) == INTSXP, Rf_error, "type of indices should be either INTSXP or REALSXP");
    make_sure(source != NULL && source != R_NilValue, Rf_error, "source must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_new\n");
        Rprintf("           SEXP: %p\n", source);
        R_xlen_t indices_length = XLENGTH(indices);
        for (R_xlen_t i = 0; i < indices_length; i++) {
            R_xlen_t value = (R_xlen_t) (TYPEOF(indices) ? REAL_ELT(indices, i) : INTEGER_ELT(indices, i));
            Rprintf("           [%l2i]: %li\n", i, value);
        }
    }

    SEXP/*LISTSXP*/ data = allocSExp(LISTSXP);
    SETCAR (data, source);     // The original vector
    SET_TAG(data, R_NilValue); // Starts as R_NilValue, becomes a vector if it the prism is written to
    SETCDR (data, R_NilValue); // Nothing here

    return R_new_altrep(class_from_sexp_type(TYPEOF(source)), indices, data);
}

static inline SEXP/*INTSXP|REALSXP*/ get_indices(SEXP x) {
    return R_altrep_data1(x);
}

static inline R_xlen_t get_length(SEXP x) {
	return XLENGTH(get_indices(x));
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

SEXP prism_duplicate(SEXP x, Rboolean deep) {//TODO

    if (get_debug_mode()) {
        Rprintf("prism_duplicate\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("           deep: %i\n", deep);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*INTSXP|REALSXP*/ indices = get_indices(x);
    SEXP                   source  = get_source(x);

    if (deep) {
        SEXP prism = prism_new(source, duplicate(indices));
        if (is_materialized(x)) {
            SEXP data = get_materialized_data(x);
            set_materialized_data(x, duplicate(data));
        }
        return prism;
    } else {
        SEXP/*LISTSXP*/ data = allocSExp(LISTSXP);
        SETCAR (data, source);
        if (is_materialized(x)) {
            SET_TAG(data, get_materialized_data(x));
        } else {
            SET_TAG(data, R_NilValue);
        }
        SETCDR (data, R_NilValue);
        return R_new_altrep(class_from_sexp_type(TYPEOF(source)), indices, data);
    }
}

static Rboolean prism_inspect(SEXP x, int pre, int deep, int pvec, void (*inspect_subtree)(SEXP, int, int, int)) {

    Rprintf("prism_altrep %s\n", type2char(TYPEOF(x)));

    inspect_subtree(R_altrep_data1(x), pre, deep, pvec);
    inspect_subtree(R_altrep_data2(x), pre, deep, pvec);

    return FALSE;
}

static R_xlen_t prism_length(SEXP x) {
    if (get_debug_mode()) {
        Rprintf("prism_length\n");
        Rprintf("           SEXP: %p\n", x);
    }

    return get_length(x);
}

static void *prism_dataptr(SEXP x, Rboolean writeable) {
    make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_dataptr\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("      writeable: %i\n", writeable);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP data = get_materialized_data(x);
        return (writeable) ? DATAPTR(data) : ((void *) DATAPTR_RO(data));
    }

    SEXP/*INTSXP|REALSXP*/  indices = get_indices(x);
    SEXP                    source  = get_source(x);
    SEXP data = copy_data_at_indices(source, indices);
    set_materialized_data(x, data);
    return writeable ? DATAPTR(data) : (void *) DATAPTR_RO(data);
}

static const void *prism_dataptr_or_null(SEXP x) {
	make_sure(x != NULL, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_dataptr_or_null\n");
        Rprintf("           SEXP: %p\n", x);
    }

    if (is_materialized(x)) {
        SEXP data = get_materialized_data(x);
        return DATAPTR_RO(data);
    }

    SEXP/*INTSXP|REALSXP*/  indices = get_indices(x);
    SEXP                    source  = get_source(x);
    SEXP data = copy_data_at_indices(source, indices);
    set_materialized_data(x, data);
    return DATAPTR_RO(data);
}

static inline R_xlen_t translate_index(SEXP/*INTSXP|REALSXP*/ indices, R_xlen_t index) {
    SEXPTYPE type = TYPEOF(indices);
    make_sure(type == REALSXP || type == INTSXP, Rf_error, "type of indices should be either INTSXP or REALSXP");

    if (type == REALSXP) {
        return (R_xlen_t) REAL_ELT   (indices, index);
    } else {
        return (R_xlen_t) INTEGER_ELT(indices, index);
    }
}

static int prism_integer_element(SEXP x, R_xlen_t i) {
    make_sure(x != R_NilValue, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_integer_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*INTSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == INTSXP, Rf_error, "type of data should be INTSXP");
        return INTEGER_ELT(data, i);
    }

    SEXP/*INTSXP|REALSXP*/ indices = get_indices(x);
    SEXP/*INTSXP*/         source  = get_source(x);

    R_xlen_t projected_index = translate_index(indices, i) - 1;

    if (get_debug_mode()) {
        Rprintf("projected_index: %li\n", projected_index);
    }

    make_sure(TYPEOF(source) == INTSXP, Rf_error, "type of source should be INTSXP");
    return INTEGER_ELT(source, projected_index);
}

static double prism_numeric_element(SEXP x, R_xlen_t i) {
    make_sure(x != R_NilValue, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_numeric_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*INTSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == REALSXP, Rf_error, "type of data should be REALSXP");
        return INTEGER_ELT(data, i);
    }

    SEXP/*INTSXP|REALSXP*/ indices = get_indices(x);
    SEXP/*REALSXP*/        source  = get_source(x);

    R_xlen_t projected_index = translate_index(indices, i) - 1;

    if (get_debug_mode()) {
        Rprintf("projected_index: %li\n", projected_index);
    }

    make_sure(TYPEOF(source) == REALSXP, Rf_error, "type of source should be REALSXP");
    return REAL_ELT(source, projected_index);
}

static Rbyte prism_raw_element(SEXP x, R_xlen_t i) {
    make_sure(x != R_NilValue, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_raw_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*RAWSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == RAWSXP, Rf_error, "type of data should be RAWSXP");
        return RAW_ELT(data, i);
    }

    SEXP/*INTSXP|REALSXP*/ indices = get_indices(x);
    SEXP/*RAWSXP*/         source  = get_source(x);

    R_xlen_t projected_index = translate_index(indices, i) - 1;

    if (get_debug_mode()) {
        Rprintf("projected_index: %li\n", projected_index);
    }

    make_sure(TYPEOF(source) == RAWSXP, Rf_error, "type of source should be RAWSXP");
    return RAW_ELT(source, projected_index);
}

static Rcomplex prism_complex_element(SEXP x, R_xlen_t i) {
    make_sure(x != R_NilValue, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_complex_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*CPLXSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == CPLXSXP, Rf_error, "type of data should be CPLXSXP");
        return COMPLEX_ELT(data, i);
    }

    SEXP/*INTSXP|REALSXP*/ indices = get_indices(x);
    SEXP/*CPLXSXP*/        source  = get_source(x);

    R_xlen_t projected_index = translate_index(indices, i) - 1;

    if (get_debug_mode()) {
        Rprintf("projected_index: %li\n", projected_index);
    }

    make_sure(TYPEOF(source) == CPLXSXP, Rf_error, "type of source should be CPLXSXP");
    return COMPLEX_ELT(source, projected_index);
}

static int prism_logical_element(SEXP x, R_xlen_t i) {
    make_sure(x != R_NilValue, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_logical_element\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    if (is_materialized(x)) {
        SEXP/*LGLSXP*/ data = get_materialized_data(x);
        make_sure(TYPEOF(data) == LGLSXP, Rf_error, "type of data should be LGLSXP");
        return LOGICAL_ELT(data, i);
    }

    SEXP/*INTSXP|REALSXP*/ indices = get_indices(x);
    SEXP/*LGLSXP*/         source  = get_source(x);

    R_xlen_t projected_index = translate_index(indices, i) - 1;

    if (get_debug_mode()) {
        Rprintf("projected_index: %li\n", projected_index);
    }

    make_sure(TYPEOF(source) == LGLSXP, Rf_error, "type of source should be LGLSXP");
    return LOGICAL_ELT(source, projected_index);
}

static R_xlen_t prism_integer_get_region(SEXP x, R_xlen_t i, R_xlen_t n, int *buf) {
    make_sure(x != R_NilValue, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_integer_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*INTSXP*/ data;
    if (!is_materialized(x)) {
        SEXP/*INTSXP|REALSXP*/ indices = get_indices(x);
        SEXP                   source  = get_source(x);
        SEXP data = copy_data_at_indices(source, indices);
        set_materialized_data(x, data);
    } else {
        data = get_materialized_data(x);
    }

    make_sure(TYPEOF(data) == INTSXP, Rf_error, "type of data should be INTSXP");
    return INTEGER_GET_REGION(data, i, n, buf);
}

static R_xlen_t prism_numeric_get_region(SEXP x, R_xlen_t i, R_xlen_t n, double *buf) {
    make_sure(x != R_NilValue, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_numeric_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*REALSXP*/ data;
    if (!is_materialized(x)) {
        SEXP/*INTSXP|REALSXP*/ indices = get_indices(x);
        SEXP                   source  = get_source(x);
        SEXP data = copy_data_at_indices(source, indices);
        set_materialized_data(x, data);
    } else {
        data = get_materialized_data(x);
    }

    make_sure(TYPEOF(data) == REALSXP, Rf_error, "type of data should be REALSXP");
    return REAL_GET_REGION(data, i, n, buf);
}

static R_xlen_t prism_raw_get_region(SEXP x, R_xlen_t i, R_xlen_t n, Rbyte *buf) {
    make_sure(x != R_NilValue, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_raw_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*RAWSXP*/ data;
    if (!is_materialized(x)) {
        SEXP/*INTSXP|REALSXP*/ indices = get_indices(x);
        SEXP                   source  = get_source(x);
        SEXP data = copy_data_at_indices(source, indices);
        set_materialized_data(x, data);
    } else {
        data = get_materialized_data(x);
    }

    make_sure(TYPEOF(data) == RAWSXP, Rf_error, "type of data should be RAWSXP");
    return RAW_GET_REGION(data, i, n, buf);
}

static R_xlen_t prism_complex_get_region(SEXP x, R_xlen_t i, R_xlen_t n, Rcomplex *buf) {
    make_sure(x != R_NilValue, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_complex_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*CPLXSXP*/ data;
    if (!is_materialized(x)) {
        SEXP/*INTSXP|REALSXP*/ indices = get_indices(x);
        SEXP                   source  = get_source(x);
        SEXP data = copy_data_at_indices(source, indices);
        set_materialized_data(x, data);
    } else {
        data = get_materialized_data(x);
    }

    make_sure(TYPEOF(data) == CPLXSXP, Rf_error, "type of data should be CPLXSXP");
    return COMPLEX_GET_REGION(data, i, n, buf);
}

static R_xlen_t prism_logical_get_region(SEXP x, R_xlen_t i, R_xlen_t n, int *buf) {
    make_sure(x != R_NilValue, Rf_error, "x must not be null");

    if (get_debug_mode()) {
        Rprintf("prism_logical_get_region\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("          index: %li\n", i);
        Rprintf("           size: %li\n", n);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*LGLSXP*/ data;
    if (!is_materialized(x)) {
        SEXP/*INTSXP|REALSXP*/ indices = get_indices(x);
        SEXP                   source  = get_source(x);
        SEXP data = copy_data_at_indices(source, indices);
        set_materialized_data(x, data);
    } else {
        data = get_materialized_data(x);
    }

    make_sure(TYPEOF(data) == LGLSXP, Rf_error, "type of data should be LGLSXP");
    return LOGICAL_GET_REGION(data, i, n, buf);
}

SEXP/*REALSXP*/ map_indices_onto_source(SEXP/*INTSXP|REALSXP*/ indices, SEXP/*INTSXP|REALSXP*/ prism_indices, R_xlen_t size) {
	SEXP translated_indices = allocVector(REALSXP, size);

	switch (TYPEOF(indices)) {
	        case INTSXP:
	            for (R_xlen_t i = 0; i < size; i++) {
	                int index = INTEGER_ELT(indices, i);
	                if (index == NA_INTEGER) {
	                	SET_REAL_ELT(translated_indices, i, NA_REAL);
	                	continue;
	                }
	                R_xlen_t translated_index = translate_index(prism_indices, ((R_xlen_t) index) - 1);
	                SET_REAL_ELT(translated_indices, i, translated_index);
	            }
	            break;

	        case REALSXP:
	            for (R_xlen_t i = 0; i < size; i++) {
	                double index = REAL_ELT(indices, i);
	                if (ISNAN(index)) {
	                	SET_REAL_ELT(translated_indices, i, NA_REAL);
	                	continue;
	                }
	                R_xlen_t translated_index = translate_index(prism_indices, ((R_xlen_t) index) - 1);
	                SET_REAL_ELT(translated_indices, i, translated_index);
	            }
	            break;
	    }

	return translated_indices;
}

static SEXP prism_extract_subset(SEXP x, SEXP indices, SEXP call) {
	make_sure(x != NULL, Rf_error, "x must not be null");
    make_sure(TYPEOF(indices) == REALSXP || TYPEOF(indices) == INTSXP, Rf_error, "type of indices should be either INTSXP or REALSXP");

    if (get_debug_mode()) {
        Rprintf("prism_extract_subset\n");
        Rprintf("           SEXP: %p\n", x);
        Rprintf("        indices: %p\n", indices);
        Rprintf("           call: %p\n", call);
        Rprintf("is_materialized: %i\n", is_materialized(x));
    }

    SEXP/*INTSXP|REALSXP*/ prism_indices = get_indices(x);
    SEXP                   source        = get_source(x);
    R_xlen_t       		   length        = get_length(x);

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

    if (!do_indices_contain_NAs(screened_indices)) {
    	SEXP/*REALSXP*/ translated_indices = map_indices_onto_source(screened_indices, prism_indices, size);
    	return copy_data_at_indices(source, translated_indices);
    }

    // Non-NA indices.
    SEXP translated_indices = map_indices_onto_source(indices, prism_indices, size);
    return prism_new(source, translated_indices);
}

// R_set_altstring_Set_elt_method
// static void string_set_elt(SEXP x, R_xlen_t i, SEXP v)

// R_set_altstring_Elt_method
// string_elt(SEXP x, R_xlen_t i)

void init_common_prism(R_altrep_class_t cls) {
    R_set_altrep_Duplicate_method(cls, prism_duplicate);
    R_set_altrep_Inspect_method(cls, prism_inspect);
    R_set_altrep_Length_method(cls, prism_length);

    R_set_altvec_Dataptr_method(cls, prism_dataptr);
    R_set_altvec_Dataptr_or_null_method(cls, prism_dataptr_or_null);
    R_set_altvec_Extract_subset_method(cls, prism_extract_subset);
}

void init_integer_prism(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altinteger_class("prism_integer_altrep", "viewports", dll);
    prism_integer_altrep = cls;

    init_common_prism(cls);

    R_set_altinteger_Elt_method(cls, prism_integer_element);
    R_set_altinteger_Get_region_method(cls, prism_integer_get_region);
}

void init_numeric_prism(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altreal_class("prism_numeric_altrep", "viewports", dll);
    prism_numeric_altrep = cls;

    init_common_prism(cls);

    R_set_altreal_Elt_method   (cls, prism_numeric_element);
    R_set_altreal_Get_region_method   (cls, prism_numeric_get_region);
}

void init_logical_prism(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altlogical_class("prism_logical_altrep", "viewports", dll);
    prism_logical_altrep = cls;

    init_common_prism(cls);

    R_set_altlogical_Elt_method(cls, prism_logical_element);
    R_set_altlogical_Get_region_method(cls, prism_logical_get_region);
}

void init_complex_prism(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altcomplex_class("prism_complex_altrep", "viewports", dll);
    prism_complex_altrep = cls;

    init_common_prism(cls);

    R_set_altcomplex_Elt_method(cls, prism_complex_element);
    R_set_altcomplex_Get_region_method(cls, prism_complex_get_region);
}

void init_raw_prism(DllInfo * dll) {
    R_altrep_class_t cls = R_make_altraw_class("prism_raw_altrep", "viewports", dll);
    prism_raw_altrep = cls;

    init_common_prism(cls);

    R_set_altraw_Elt_method    (cls, prism_raw_element);
    R_set_altraw_Get_region_method    (cls, prism_raw_get_region);
}

void init_prism_altrep_class(DllInfo * dll) {
    init_integer_prism(dll);
    init_numeric_prism(dll);
    init_logical_prism(dll);
    init_complex_prism(dll);
    init_raw_prism(dll);
}

SEXP create_prism(SEXP source, SEXP/*INTSXP|REALSXP*/ indices) {
	SEXPTYPE source_type = TYPEOF(source);
	SEXPTYPE indices_type = TYPEOF(indices);
	R_xlen_t source_length = XLENGTH(source);
	//R_xlen_t indices_length = XLENGTH(indices);

	make_sure(source_type == INTSXP || source_type == REALSXP || source_type == CPLXSXP
	           || source_type == LGLSXP || source_type == RAWSXP  || source_type == VECSXP
	           || source_type == STRSXP, Rf_error,
			      "type of source must be one of INTSXP, REALSXP, RAWSXP, CPLXSXP, LGLSXP, VECSXP, or STRSXP");

	make_sure(indices_type == REALSXP || indices_type == INTSXP, Rf_error,
			  "type of indices should be either INTSXP or REALSXP");

  	if (!are_indices_in_range(indices, 0, source_length))
    		Rf_error("Cannot use these indices with this source: out of range");

    if (get_debug_mode()) {
        Rprintf("create prism\n");
        Rprintf("           SEXP: %p\n", source);
    }

    return prism_new(source, indices);
}
