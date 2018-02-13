library(testthat)
library(genBart)

context("columnname warnings")
test_that("warnings are thrown when columnname is NULL or misspecified", {
  expect_warning(metaData(y = tb.expr, design = tb.design), "Must enter columnname to match design and expression")
  expect_warning(metaData(y = tb.expr, design = tb.design, columnname = "sample_id"), 
                 "columnname values in design do not match any columnnames in y")
  expect_warning(metaData(y = tb.expr, design = tb.design, columnname = "Columnname"),
                 "columnname parameter given is not a column name in design. 
Please check spelling.")
}) 

context("design and expression matching")
test_that("design matches expression even when original inputs do not match", {
  test.design <- tb.design[sample(1:nrow(tb.design), nrow(tb.design)), ]
  meta <- metaData(y = tb.expr, design = test.design, columnname = "columnname")
  expect_equal(as.character(meta$design$columnname), colnames(meta$y)) 
  meta <- metaData(y = tb.expr, design = test.design[-c(sample(1:nrow(test.design), 5)), ], columnname = "columnname")
  expect_message(metaData(y = tb.expr, design = test.design[-c(sample(1:nrow(test.design), 5)), ], 
                          columnname = "columnname"), "More samples in y than design. Throwing out excess samples.")
  expect_equal(as.character(meta$design$columnname), colnames(meta$y))
  expect_message(metaData(y = tb.expr[, -c(sample(1:ncol(tb.expr), 5))], design = test.design,
                          columnname = "columnname"), "More samples in design than y. Throwing out excess samples.")
  meta <- metaData(y = tb.expr[, -c(sample(1:ncol(tb.expr), 5))], design = test.design, columnname = "columnname")
  expect_equal(as.character(meta$design$columnname), colnames(meta$y))
  test.design <- tb.design[-c(sample(1:nrow(tb.design), 5)), ]
  test.design$columnname <- as.character(test.design$columnname)
  test.design$columnname[1:3] <- c("newSample1", "newSample2", "newSample3")
  expect_message(metaData(y = tb.expr, design = test.design, columnname = "columnname"), 
                 "More samples in y than design, but some samples in design are not 
in y. Throwing out excess and unmatched samples from design and y.")
  meta <- metaData(y = tb.expr, design = test.design, columnname = "columnname")
  expect_equal(as.character(meta$design$columnname), colnames(meta$y))
  test.expr <- tb.expr[, -c(sample(1:ncol(tb.expr), 5))]
  colnames(test.expr)[1:3] <- c("newSample1", "newSample2", "newSample3")
  expect_message(metaData(y = test.expr, design = tb.design, columnname = "columnname"), 
                 "More samples in design than y, but some samples in y are not in 
design. Throwing out excess and unmatched samples from design and y.")
  meta <- metaData(y = tb.expr, design = test.design, columnname = "columnname")
  expect_equal(as.character(meta$design$columnname), colnames(meta$y))
  test.expr <- tb.expr
  colnames(test.expr)[1:3] <- c("newSample1", "newSample2", "newSample3")
  test.design <- tb.design
  test.design$columnname <- as.character(test.design$columnname)
  test.design$columnname[1:3] <- c("newSample1", "newSample2", "newSample3")
  expect_message(metaData(y = test.expr, design = tb.design, columnname = "columnname"),
                 "There are an equal number of samples in design and y, but some 
sample names do not match. Throwing out unmatched samples.")
  meta <- metaData(y = test.expr, design = tb.design, columnname = "columnname")
  expect_equal(as.character(meta$design$columnname), colnames(meta$y))
  expect_message(metaData(y = tb.expr, design = test.design, columnname = "columnname"),
                 "There are an equal number of samples in design and y, but some 
sample names do not match. Throwing out unmatched samples.")
  meta <- metaData(y = tb.expr, design = test.design, columnname = "columnname")
  expect_equal(as.character(meta$design$columnname), colnames(meta$y))
})

context("design with controls")
test_that("designs with controls are specified correctly", {
  expect_warning(metaData(y = tb.expr, design = tb.design, columnname = "columnname", control.var = "clinical_status"), 
                 "control.var specified but not control.val. Please specify 
control.val")
  expect_warning(metaData(y = tb.expr, design = tb.design, columnname = "columnname", control.val = "Latent"), 
                 "control.val specified but not control.var. Please specify 
control.var")
  expect_warning(metaData(y = tb.expr, design = tb.design, columnname = "columnname", control.var = "Clinical_status"),
                 "control.var parameter given is not a column name in design. 
Please check spelling.")
})

context("sample.id and subject.id correctly specified")
test_that("designs with controls are specified correctly", {
  expect_warning(metaData(y = tb.expr, design = tb.design, columnname = "columnname", control.var = "clinical_status",
                          control.val = "Latent", sample.id = "sample.id"), "sample.id parameter given is not a column name in design. 
Please check spelling."
                )
  expect_warning(metaData(y = tb.expr, design = tb.design, columnname = "columnname", control.var = "clinical_status",
                          control.val = "Latent", sample.id = "sample_id", subject.id = "monkey.id"), "subject.id parameter given is not a column name in design. 
Please check spelling."
  )
})

expect_all_warnings <- function(warnings, test.warnings) {
  expect_equal(warnings, test.warnings)
}

context("longitudinal design warnings/messages")
test_that("longitudinal designs are specified correctly", {
  expect_all_warnings(capture_warnings(metaData(y = tb.expr, design = tb.design, columnname = "columnname",
                         long = TRUE)), c("long is TRUE but time.var not specified", 
                                          "long is TRUE but baseline.var not specified", 
                                          "long is TRUE but baseline.val not specified"))
  expect_all_warnings(capture_warnings(metaData(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, time.var = "timepoint")), c("long is TRUE but baseline.var not specified", 
                                                                                        "long is TRUE but baseline.val not specified"))
  expect_warning(metaData(y = tb.expr, design = tb.design, columnname = "columnname", long = TRUE, time.var = "Timepoint"),
                 "time.var parameter given is not a column name in design. 
Please check spelling.")
  expect_all_warnings(capture_warnings(metaData(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, time.var = "timepoint", baseline.var = "timepoint")), 
                      c("long is TRUE but baseline.val not specified"))
  expect_warning(metaData(y = tb.expr, design = tb.design, columnname = "columnname", long = TRUE, time.var = "timepoint",
                          baseline.var = "Timepoint"), "baseline.var parameter given is not a column name in design. 
Please check spelling.")
  expect_all_warnings(capture_warnings(metaData(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, time.var = "timepoint", baseline.val = 0)), 
                      c("long is TRUE but baseline.var not specified"))
  expect_all_warnings(capture_warnings(metaData(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, baseline.var = "timepoint", baseline.val = 0)), 
                      c("long is TRUE but time.var not specified"))
  expect_all_warnings(capture_warnings(metaData(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, baseline.var = "timepoint")), 
                      c("long is TRUE but time.var not specified", "long is TRUE but baseline.val not specified"))
  expect_all_warnings(capture_warnings(metaData(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, baseline.val = 0)), 
                      c("long is TRUE but time.var not specified", "long is TRUE but baseline.var not specified"))
  expect_all_warnings(capture_warnings(metaData(y = tb.expr, design = tb.design, columnname = "columnname",
                                               long = TRUE, time.var = "timepoint", baseline.var = "timepoint", baseline.val = 0)), 
                      capture_warnings(5==5))
  expect_message(metaData(y = tb.expr, design = tb.design, columnname = "columnname", time.var = "timepoint", 
                         baseline.var = "timepoint", baseline.val = 0), "time.var specified. Setting long to TRUE.")
  expect_message(metaData(y = tb.expr, design = tb.design, columnname = "columnname", time.var = "timepoint"),
                 "time.var specified but not baseline.var and baseline.val")
  expect_warning(metaData(y = tb.expr, design = tb.design, columnname = "columnname", time.var = "timepoint", 
                         baseline.var = "timepoint"), "baseline.var specified but not baseline.val. Please specify 
baseline.val.")
  expect_warning(metaData(y = tb.expr, design = tb.design, columnname = "columnname", time.var = "timepoint", 
                         baseline.val = 0), "baseline.val specified but not baseline.var. Please specify 
baseline.var.")
})





