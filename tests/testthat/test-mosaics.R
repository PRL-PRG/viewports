context("Mosaic viewports")

test_that("create a mosaic from integer indices", {
    source   <- 1:10000
    viewport <- mosaic(source, as.integer(5:10))

    expect_equal(length(viewport), 6)
    expect_vector(viewport)
    expect_type(viewport, "integer")

    for (i in 1:5) { expect_equal(viewport[i], 5+i-1) }
    expect_equal(viewport, 5:10)
})

test_that("create a mosaic from numeric indices", {
    source   <- 1:10000
    viewport <- mosaic(source, as.numeric(5:10))

    expect_type(viewport, "integer")
    expect_equal(length(viewport), 6)
    expect_vector(viewport)

    for (i in 1:5) { expect_equal(viewport[i], 5+i-1) }
    expect_equal(viewport, 5:10)
})

test_that("create a mosaic from logical mask", {
    source   <- 1:10000
    mask     <- mapply(function(x) x >= 5 && x <= 10, 1:10000)
    viewport <- mosaic(source, mask)

    expect_type(viewport, "integer")
    expect_vector(viewport)

    for (i in 1:5) { expect_equal(viewport[i], 5+i-1) }
    expect_equal(viewport, 5:10)
})

test_that("write to source", {
    source <- as.integer(1:10000)
    viewport <- mosaic(source, 100:200)

    expect_equal(viewport[1], 100)
    source[100] <- 42
    expect_equal(viewport[1], 100)  # Because source was copied on write and now is not the same vector as when viewport
                                    # was created. Boo.
})

test_that("write to mosaic", {
    source <- 1:10000
    viewport <- mosaic(source, 100:200)

    expect_equal(viewport[1], 100)
    viewport[1] <- 42
    expect_equal(viewport[1], 42)
})

test_that("check range", {
     source <- 1:100
     expect_error(mosaic(source, 200:400))
})

test_that("dereference out of range", {
    source <- 1:100
    viewport <- mosaic(source, 10:20)

    expect_equal(viewport[100], as.integer(NA))
})

test_that("region partially out of range", {
    source <- 1:100
    viewport <- mosaic(source, 10:19)
    definition <- source[10:19]

    expect_equal(viewport[5:15], definition[5:15])
})

test_that("integer test", {
    source <- as.integer(c(1,2,3,4,5,6,7,8,9,10))
    viewport <- mosaic(source, 4:7)

    expect_vector(viewport)
    expect_type(viewport, "integer")

    expect_equal(viewport[1], 4)
    expect_equal(viewport[2], 5)
    expect_equal(viewport[3], 6)
    expect_equal(viewport[4], 7)
})

test_that("numeric test", {
    source <- as.numeric(c(1,2,3,4,5,6,7,8,9,10))
    viewport <- mosaic(source, 4:7)

    expect_vector(viewport)
    expect_type(viewport, "double")

    expect_equal(viewport[1], 4)
    expect_equal(viewport[2], 5)
    expect_equal(viewport[3], 6)
    expect_equal(viewport[4], 7)
})

test_that("logical test", {
    source <- as.logical(c(T, F, T, F, T, F, T, F, T, F, T, F))
    viewport <- mosaic(source, 4:7)

    expect_vector(viewport)
    expect_type(viewport, "logical")

    expect_equal(viewport[1], F)
    expect_equal(viewport[2], T)
    expect_equal(viewport[3], F)
    expect_equal(viewport[4], T)
})

test_that("complex test", {
    source <- as.complex(c(1,2,3,4,5,6,7,8,9,10))
    viewport <- mosaic(source, 4:7)

    expect_vector(viewport)
    expect_type(viewport, "complex")

    expect_equal(viewport[1], 4+0i)
    expect_equal(viewport[2], 5+0i)
    expect_equal(viewport[3], 6+0i)
    expect_equal(viewport[4], 7+0i)
})

test_that("raw test", {
    source <- as.raw(c(1,2,3,4,5,6,7,8,9,10))
    viewport <- mosaic(source, 4:7)

    expect_vector(viewport)
    expect_type(viewport, "raw")

    expect_equal(viewport[1], as.raw(04))
    expect_equal(viewport[2], as.raw(05))
    expect_equal(viewport[3], as.raw(06))
    expect_equal(viewport[4], as.raw(07))
})

test_that("sum integer test", {
	source <- as.integer(1:100)
    viewport <- mosaic(source, 10:19)
    
    expect_equal(sum(viewport), sum(source[10:19]))
})

test_that("sum numeric test", {
	source <- as.numeric(1:100)
    viewport <- mosaic(source, 10:19)
    
    expect_equal(sum(viewport), sum(source[10:19]))
})

test_that("sum boolean test", {
	source <- rep(c(T,F), 100)
    viewport <- mosaic(source, 10:19)
    
    expect_equal(sum(viewport), sum(source[10:19]))
})

# test_that("sum raw test", {
# 	source <- as.raw(1:100)
#   viewport <- mosaic(source, 10:19)
#     
#   expect_equal(sum(viewport), sum(source[10:19]))
# })

test_that("sum complex test", {
	source <- as.complex(1:100)
    viewport <- mosaic(source, 10:19)
    
    expect_equal(sum(viewport), sum(source[10:19]))
})

test_that("subset contiguous monotonic", {
    source <- 1:100
    viewport <- mosaic(source, 10:19)
    expect_equal(viewport[2:4], source[10:19][2:4])
})

test_that("subset non-contiguous monotonic", {
    source <- 1:100
    viewport <- mosaic(source, 10:19)
    expect_equal(viewport[c(1,3,5)], source[10:19][c(1,3,5)])
})

test_that("subset non-contiguous non-monotonic", {
    source <- 1:100
    viewport <- mosaic(source, 10:19)
    expect_equal(viewport[c(1,5,3)], source[10:19][c(1,5,3)])
})
