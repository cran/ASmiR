test_that("RNA sequence validation works", {
  d <- data.frame(sequence_id = c("a", "b", "c"), sequence = c("AUGC", "augu", "AUNT"))
  v <- asmir_validate_sequences(d)
  expect_true(v$valid[1])
  expect_true(v$valid[2])
  expect_false(v$valid[3])
})

test_that("bundled model files are discoverable", {
  p <- asmir_model_files("all")
  expect_length(p, 8)
  expect_true(all(file.exists(p)))
})
