

test_that("prior string editing works correctly", {
  expect_equal(minus.append("moi", 0.0), "moi")
  expect_equal(minus.append("moi", 3), "moi - 3")
  expect_equal(minus.prepend("joo", 0), "joo")
  expect_equal(minus.prepend("joo", 1), " - (joo)")
})
