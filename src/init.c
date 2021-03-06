#include "debug.h"
#include "slices.h"
#include "mosaics.h"
#include "prisms.h"

#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

// List of functions provided by the package.
static const R_CallMethodDef CallEntries[] = {

    {"slice",  (DL_FUNC) &create_slice, 3},
    {"mosaic",  (DL_FUNC) &create_mosaic, 2},
    {"prism",  (DL_FUNC) &create_prism, 2},

    // Turn on debug mode.
    {"viewport_set_debug_mode",  (DL_FUNC) &set_debug_mode,  1},

    // Terminates the function list. Necessary.
    {NULL, NULL, 0} 
};

void attribute_visible R_init_viewports(DllInfo *dll) {
    init_slice_altrep_class(dll);
    init_mosaic_altrep_class(dll);
    init_prism_altrep_class(dll);
}

void attribute_visible R_unload_viewports(DllInfo *dll) {
}
