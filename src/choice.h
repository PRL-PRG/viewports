#pragma once

typedef enum {
    VIEWPORT_NONE,
    VIEWPORT_SLICE,
    VIEWPORT_MOSAIC,
    VIEWPORT_PRISM,
} viewport_type_t;

viewport_type_t select_best_viewport_type    (SEXP source, SEXP indices);
SEXP            create_best_viewport_or_clone(SEXP source, SEXP indices);