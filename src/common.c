#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>

#include "common.h"

#define MAKE_SURE
#include "make_sure.h"

bool are_integer_indices_monotonic(SEXP/*INTSXP*/ indices) {
    make_sure(TYPEOF(indices) == INTSXP, Rf_error, "type of indices must be INTSXP");

    R_xlen_t size = XLENGTH(indices);
    int previous = NA_INTEGER;

    for (int i = 0; i < size; i++) {
        int current = INTEGER_ELT(indices, i);
        if (current == NA_INTEGER) {
            return false;
        }
        if (previous != NA_INTEGER && previous >= current) {
            return false;
        }
        previous = current;
    }
    return true;
}

bool are_numeric_indices_monotonic(SEXP/*REALSXP*/ indices) {
	make_sure(TYPEOF(indices) == REALSXP, Rf_error, "type of indices should be REALSXP");

    R_xlen_t size = XLENGTH(indices);
    double previous = NA_REAL;

    for (int i = 0; i < size; i++) {
        double current = REAL_ELT(indices, i);
        if (ISNAN(current)) {
            return false;
        }
        if (!ISNAN(previous) && previous >= current) {
            return false;
        }
        previous = current;
    }
    return true;
}

bool are_indices_monotonic (SEXP/*INTSXP | REALSXP*/ indices) {
    SEXPTYPE type = TYPEOF(indices);
    make_sure(type == INTSXP || type == REALSXP, Rf_error, "type of indices should be either INTSXP or REALSXP");

    switch (type) {
        case INTSXP:  return are_integer_indices_monotonic(indices);
        case REALSXP: return are_numeric_indices_monotonic(indices);
        default:      Rf_error("Slices can be indexed by integer or numeric vectors but found: %d\n", type);
    }
}

bool are_integer_indices_contiguous(SEXP/*INTSXP*/ indices) {
	make_sure(TYPEOF(indices) == INTSXP, Rf_error, "type of indices should be INTSXP");

    R_xlen_t size = XLENGTH(indices);
    int previous = NA_INTEGER;

    for (int i = 0; i < size; i++) {
        int current = INTEGER_ELT(indices, i);
        if (current == NA_INTEGER) {
            return false;
        }
        if (previous != NA_INTEGER && previous != current - 1) {
            return false;
        }
        previous = current;
    }

    return true;
}

bool are_numeric_indices_contiguous(SEXP/*REALSXP*/ indices) {
	make_sure(TYPEOF(indices) == REALSXP, Rf_error, "type of indices should be REALSXP");

    R_xlen_t size = XLENGTH(indices);
    double previous = NA_REAL;

    for (int i = 0; i < size; i++) {
        double current = REAL_ELT(indices, i);
        if (ISNAN(current)) {
            return false;
        }
        if (!ISNAN(previous) && previous != current + 1) {
            return false;
        }
        previous = current;
    }

    return true;
}

bool are_indices_contiguous(SEXP/*INTSXP | REALSXP*/ indices) {
    SEXPTYPE type = TYPEOF(indices);
    make_sure(type == INTSXP || type == REALSXP, Rf_error, "type of indices should be either INTSXP or REALSXP");

    switch (type) {
        case INTSXP:  return are_integer_indices_contiguous(indices);
        case REALSXP: return are_numeric_indices_contiguous(indices);
        default:      Rf_error("Slices can be indexed by integer or numeric vectors but found: %d\n", type);
    }
}

bool are_integer_indices_in_range(SEXP/*INTSXP*/ indices, R_xlen_t min, R_xlen_t max) {
    make_sure(TYPEOF(indices) == INTSXP, Rf_error, "type of indices should be INTSXP");
    make_sure(min <= max, Rf_error, "min should be less or equal to max");

    R_xlen_t size = XLENGTH(indices);
    for (int i = 0; i < size; i++) {
        int current = INTEGER_ELT(indices, i);
        if (current < min || current > max) {
            return false;
        }
    }

    return true;
}

bool are_numeric_indices_in_range(SEXP/*REALSXP*/ indices, R_xlen_t min, R_xlen_t max) {
	make_sure(TYPEOF(indices) == REALSXP, Rf_error, "type of indices should be REALSXP");
    make_sure(min <= max, Rf_error, "min must be less or equal to max");

    R_xlen_t size = XLENGTH(indices);
    for (int i = 0; i < size; i++) {
        double current = REAL_ELT(indices, i);
        if (current < min || current > max) {
            return false;
        }
    }

    return true;
}

bool are_indices_in_range(SEXP/*INTSXP | REALSXP*/ indices, R_xlen_t min, R_xlen_t max) {
    SEXPTYPE type = TYPEOF(indices);
    make_sure(type == INTSXP || type == REALSXP, Rf_error, "type of indices should be either INTSXP or REALSXP");
    make_sure(min <= max, Rf_error, "min must be less or equal to max");

    switch (type) {
        case INTSXP:  return are_integer_indices_in_range(indices, min, max);
        case REALSXP: return are_numeric_indices_in_range(indices, min, max);
        default:      Rf_error("Slices can be indexed by integer or numeric vectors but found: %d\n", type);
    }
}

bool do_integer_indices_contain_NAs(SEXP/*INTSXP*/ indices) {
	make_sure(TYPEOF(indices) == INTSXP, Rf_error, "type of indices should be INTSXP");
    R_xlen_t size = XLENGTH(indices);
    for (int i = 0; i < size; i++) {
        if (INTEGER_ELT(indices, i) == NA_INTEGER) {
            return false;
        }
    }
    return true;
}

bool do_numeric_indices_contain_NAs(SEXP/*REALSXP*/ indices) {
	make_sure(TYPEOF(indices) == REALSXP, Rf_error, "type of indices should be REALSXP");
    R_xlen_t size = XLENGTH(indices);
    for (int i = 0; i < size; i++) {
        if (ISNAN(REAL_ELT(indices, i))) {
            return false;
        }
    }
    return true;
}

bool do_indices_contain_NAs(SEXP/*INTSXP | REALSXP*/ indices) {
    SEXPTYPE type = TYPEOF(indices);
    make_sure(type == INTSXP || type == REALSXP, Rf_error, "type of indices should be either INTSXP or REALSXP");

    switch (type) {
        case INTSXP:  return do_integer_indices_contain_NAs(indices);
        case REALSXP: return do_numeric_indices_contain_NAs(indices);
        default:      Rf_error("Slices can be indexed by integer or numeric vectors but found: %d\n", type);
    }
}

R_xlen_t get_first_element_as_length(SEXP/*INTSXP | REALSXP*/ indices) {
    SEXPTYPE type = TYPEOF(indices);
    make_sure(type == INTSXP || type == REALSXP, Rf_error, "type of indices should be either INTSXP or REALSXP");
    make_sure(XLENGTH(indices) > 0, Rf_error, "indices cannot be empty");

    switch (type) {
        case INTSXP:  return (R_xlen_t) INTEGER_ELT(indices, 0);
        case REALSXP: return (R_xlen_t) REAL_ELT(indices, 0);
        default:      Rf_error("Slices can be indexed by integer or numeric vectors but found: %d\n", type);
    }
}

void copy_element(SEXP source, R_xlen_t source_index, SEXP target, R_xlen_t target_index) {
    make_sure(TYPEOF(source) == TYPEOF(target), Rf_error, "type of source and target must be the same");
    make_sure(TYPEOF(source) == INTSXP  || TYPEOF(source) == REALSXP || TYPEOF(source) == CPLXSXP
            || TYPEOF(source) == LGLSXP || TYPEOF(target) == RAWSXP  || TYPEOF(source) == VECSXP
            || TYPEOF(source) == STRSXP,
			Rf_error, "type of source must be one of INTSXP, REALSXP, CPLXSXP, LGLSXP, RAWSXP, VECSXP, or STRSXP");

    switch(TYPEOF(source)) {
        case INTSXP:  SET_INTEGER_ELT (target, target_index, INTEGER_ELT (source, source_index)); break;
        case REALSXP: SET_REAL_ELT    (target, target_index, REAL_ELT    (source, source_index)); break;
        case LGLSXP:  SET_LOGICAL_ELT (target, target_index, LOGICAL_ELT (source, source_index)); break;
        case CPLXSXP: SET_COMPLEX_ELT (target, target_index, COMPLEX_ELT (source, source_index)); break;
        case RAWSXP:  SET_RAW_ELT     (target, target_index, RAW_ELT     (source, source_index)); break;
        case STRSXP:  SET_VECTOR_ELT  (target, target_index, VECTOR_ELT  (source, source_index)); break;
        case VECSXP:  SET_VECTOR_ELT  (target, target_index, VECTOR_ELT  (source, source_index)); break;
        default:      Rf_error("Unsupported vector type: %d\n", TYPEOF(source));
    }
}

void set_element_to_NA(SEXP target, R_xlen_t target_index) {

	make_sure(TYPEOF(target) == INTSXP  || TYPEOF(target) == REALSXP || TYPEOF(target) == CPLXSXP
	       || TYPEOF(target) == LGLSXP  || TYPEOF(target) == STRSXP,
		   Rf_error, "type of source must be one of INTSXP, REALSXP, CPLXSXP, LGLSXP, or STRSXP");

    switch(TYPEOF(target)) {
        case INTSXP:  SET_INTEGER_ELT (target, target_index, NA_INTEGER); break;
        case REALSXP: SET_REAL_ELT    (target, target_index, NA_REAL);    break;
        case LGLSXP:  SET_LOGICAL_ELT (target, target_index, NA_LOGICAL); break;
        case CPLXSXP: {
                      Rcomplex NA_CPLX = { NA_REAL, NA_REAL };
                      SET_COMPLEX_ELT (target, target_index, NA_CPLX);    break;
        }
        case STRSXP:  SET_STRING_ELT  (target, target_index, NA_STRING);  break;
        default:      Rf_error("Unsupported vector type: %d\n", TYPEOF(target));
    }
}

SEXP copy_data_in_range(SEXP source, R_xlen_t start, R_xlen_t size) {

    SEXPTYPE type = TYPEOF(source);
    SEXP target = allocVector(type, size);

    for (R_xlen_t i = 0; i < size; i++) {
        copy_element(source, start + i, target, i);
    }

    return target;
}

SEXP copy_data_at_indices(SEXP source, SEXP/*INTSXP | REALSXP*/ indices) {

    SEXPTYPE type = TYPEOF(indices);
    make_sure(type == INTSXP || type == REALSXP, Rf_error, "type of indices should be either INTSXP or REALSXP");

    R_xlen_t size = XLENGTH(indices);
    SEXP target = allocVector(TYPEOF(source), size);

    switch (type) {
        case INTSXP:  {
            for (R_xlen_t i = 0; i < size; i++) {
                int index = INTEGER_ELT(indices, i);
                if (index == NA_INTEGER) {
                	set_element_to_NA(target, i);
                } else {
                	copy_element(source, ((R_xlen_t) index) - 1, target, i);
                }
            }
            break;
        }

        case REALSXP: {
            for (R_xlen_t i = 0; i < size; i++) {
                double index = REAL_ELT(indices, i);
                if (ISNAN(index)) {
                  	set_element_to_NA(target, i);
                } else {
                    copy_element(source, ((R_xlen_t) index) - 1, target, i);
                }
            }
            break;
        }

        default:
        	Rf_error("Unsupported vector type: %d\n", type);
    }

    return target;
}

SEXP copy_data_at_mask(SEXP source, SEXP/*LGLSXP*/ mask) {

    SEXPTYPE type = TYPEOF(mask);
    make_sure(type == LGLSXP, Rf_error, "type of mask must be LGLSXP");

    R_xlen_t mask_size = XLENGTH(mask);
    R_xlen_t target_size = mask_size;
    for (R_xlen_t i = 0; i < mask_size; i++) {
        if (LOGICAL_ELT(mask, i) == FALSE) {
            target_size--;
        }
    }

    SEXP target = allocVector(type, target_size);
    R_xlen_t copied_elements = 0;

    for (R_xlen_t index = 0; index < mask_size; index++) {
        Rboolean current = LOGICAL_ELT(mask, index);

        if (current == NA_LOGICAL) {
            set_element_to_NA(target, copied_elements);
            copied_elements++;
            continue;
        }

        if (current == TRUE) {
            copy_element(source, index, target, copied_elements);
            copied_elements++;
        }
    }

    make_sure(XLENGTH(target) == copied_elements, Rf_error,
    		  "the number of copied elements is different than the size of the output vector");
    return target;
}

viewport_type_t recommend_vieport_type_for_indices(SEXP/*INTSXP | REALSXP*/ indices) {
    make_sure(false, Rf_error, "not implemented"); //FIXME
    return VIEWPORT_NONE;
}
