library(testthat)
library(genBART)

context("columnname warnings")
test_that("warnings are thrown when columnname is NULL or misspecified", {
  expect_warning(desInfo(y = tb.expr, design = tb.design), "Must enter columnname to match design and expression")
  expect_warning(desInfo(y = tb.expr, design = tb.design, columnname = "sample_id"), 
                 "columnname values in design do not match any columnnames in y")
})

context("design and expression matching")
test_that("design matches expression even when original inputs do not match", {
  test.design <- tb.design
  test.design <- test.design[sample(1:nrow(test.design), nrow(test.design)), ]
  des.info <- desInfo(y = tb.expr, design = test.design, columnname = "columnname")
  expect_equal(as.character(des.info$design$columnname), colnames(des.info$y))
  des.info <- desInfo(y = tb.expr, design = test.design[-c(sample(1:nrow(test.design), 5)), ],
                      columnname = "columnname")
  expect_message(desInfo(y = tb.expr, design = test.design[-c(sample(1:nrow(test.design), 5)), ],
                         columnname = "columnname"), "More samples in y than design. Throwing out excess samples.")
  expect_equal(as.character(des.info$design$columnname), colnames(des.info$y))
  expect_message(desInfo(y = tb.expr[, -c(sample(1:ncol(tb.expr), 5))], design = test.design,
                         columnname = "columnname"), "More samples in design than y. Throwing out excess samples.")
  des.info <- desInfo(y = tb.expr[, -c(sample(1:ncol(tb.expr), 5))], design = test.design,
                      columnname = "columnname")
  expect_equal(as.character(des.info$design$columnname), colnames(des.info$y))
})

context("hc correctly returned")
test_that("hc is returned correctly", {
  hc <- desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                control_var = "clinical_status", control_val = "Latent")$hc
  expect_true(hc)
  hc <- suppressWarnings(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                control_var = "clinical_status")$hc)
  expect_false(hc)
  hc <- suppressWarnings(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                control_val = "Latent")$hc)
  expect_false(hc)
  expect_warning(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                         control_var = "clinical_status"), "control_var specified without control_val")
  expect_warning(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                         control_val = "Latent"), "control_val specified without control_var")
  
})

expect_all_warnings <- function(warnings, test.warnings) {
  expect_equal(warnings, test.warnings)
}

context("long design warnings/messages")
test_that("longitudinal designs are specified correctly", {
  expect_all_warnings(capture_warnings(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                         long = TRUE)), c("long is TRUE but time_var not specified", 
                                          "long is TRUE but baseline_var not specified", 
                                          "long is TRUE but baseline_val not specified"))
  expect_all_warnings(capture_warnings(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, time_var = "timepoint")), c("long is TRUE but baseline_var not specified", 
                                                                                        "long is TRUE but baseline_val not specified"))
  expect_all_warnings(capture_warnings(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, time_var = "timepoint", baseline_var = "timepoint")), 
                      c("long is TRUE but baseline_val not specified"))
  expect_all_warnings(capture_warnings(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, time_var = "timepoint", baseline_val = 0)), 
                      c("long is TRUE but baseline_var not specified"))
  expect_all_warnings(capture_warnings(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, baseline_var = "timepoint", baseline_val = 0)), 
                      c("long is TRUE but time_var not specified"))
  expect_all_warnings(capture_warnings(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, baseline_var = "timepoint")), 
                      c("long is TRUE but time_var not specified", "long is TRUE but baseline_val not specified"))
  expect_all_warnings(capture_warnings(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, baseline_val = 0)), 
                      c("long is TRUE but time_var not specified", "long is TRUE but baseline_var not specified"))
  expect_all_warnings(capture_warnings(desInfo(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, time_var = "timepoint", baseline_var = "timepoint", baseline_val = 0)), 
                      capture_warnings(5==5))
  expect_message(desInfo(y = tb.expr, design = tb.design, columnname = "columnname", time_var = "timepoint", 
                         baseline_var = "timepoint", baseline_val = 0), "time_var specified. Setting long to TRUE.")
  expect_message(desInfo(y = tb.expr, design = tb.design, columnname = "columnname", time_var = "timepoint"),
                 "time_var specified without baseline_var and baseline_val")
  expect_warning(desInfo(y = tb.expr, design = tb.design, columnname = "columnname", time_var = "timepoint", 
                         baseline_var = "timepoint"), "baseline_var specified without baseline_val")
  expect_warning(desInfo(y = tb.expr, design = tb.design, columnname = "columnname", time_var = "timepoint", 
                         baseline_val = 0), "baseline_val specified without baseline_var")
})





