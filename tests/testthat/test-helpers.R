# Tests for the helpers functions :O

test_that("The largest object is kept", {
  testmat = matrix(0, nrow = 100, ncol = 10)
  testmat[,6:10 ] <- 1
  testmat2 = keep_largest(testmat)
  expect_equal(testmat, testmat2)
})


test_that("Don't fill non-existent holes", {
  testmat = matrix(0, nrow = 100, ncol = 10)
  testmat[2:49,  ] <- 1
  testmat2 = fill_holes(testmat)
  expect_equal(testmat, testmat2)
  testmat = (!testmat)*1L
  testmat2 = fill_holes(testmat)
  expect_equal(testmat, testmat2)

})

test_that("Get the right cardinal direction", {
  mat = matrix(0, nrow = 100, ncol = 100)
  expect_equal(card_select(mat, 1), matrix(c(1,1), nrow = 1))
  expect_equal(card_select(mat, 5), matrix(c(100,1), nrow = 1))
  expect_equal(card_select(mat, 9), matrix(c(100,100), nrow = 1))
  expect_equal(card_select(mat, 13), matrix(c(1,100), nrow = 1))
})
