# CLAIM (README, Dataset Groups): the per-file summary gives one row per file (per
# channel facet, per metric) tagged with the file's group, so between-group
# statistics can use the cell as the unit of replication.
test_that("per_file_summary resolves facets to files and attaches groups", {
  nn <- data.frame(distance = c(1:10, 11:15, 21:30),
                   dataset = rep(c("a.npy", "b.npy", "c.npy \u00b7 Ch 1"), c(10, 5, 10)))
  out <- per_file_summary(list(batch_nn = list(df = nn, metrics = "distance", unit = "NN"),
                               rog = list(df = NULL, metrics = "rog")),
                          group_map = c(a.npy = "WT", c.npy = "KO"),
                          files = c("a.npy", "b.npy", "c.npy"))
  expect_equal(nrow(out), 3L)
  expect_equal(out$file[out$dataset == "c.npy \u00b7 Ch 1"], "c.npy")
  expect_equal(out$group[out$file == "c.npy"], "KO")
  expect_true(is.na(out$group[out$file == "b.npy"]))           # unassigned file keeps a row
  expect_equal(out$n[out$file == "a.npy"], 10L)
  expect_equal(out$median[out$file == "a.npy"], 5.5)
})

test_that("per_file_summary returns an empty frame when nothing has been run", {
  expect_equal(nrow(per_file_summary(list())), 0L)
  expect_equal(nrow(per_file_summary(list(x = list(df = NULL, metrics = "d")))), 0L)
})
