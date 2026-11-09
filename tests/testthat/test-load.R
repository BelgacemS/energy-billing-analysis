test_that("parse_kwh removes thousands separators", {
  expect_equal(parse_kwh(c("2 500", "687", "590 123", "1 234")),
               c(2500, 687, 590123, 1234))
  expect_error(parse_kwh(c("12", "abc")), "Unparsable")
})

test_that("days_in_month handles leap years and December", {
  d <- as.Date(c("2019-02-10", "2020-02-01", "2021-09-15", "2021-12-31"))
  expect_equal(days_in_month(d), c(28L, 29L, 30L, 31L))
})

test_that("load_bills normalises a small hand-made file", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("NAME,BEGIN,END,KWH_PERIOD",
               "0002KO,01/02/2020,28/02/2020,1 160",
               "0001BR,15/09/2021,30/09/2021,300",
               "0001BR,01/08/2021,31/08/2021,3 100"), path)
  b <- load_bills(path)
  expect_equal(b$site, c("0001BR", "0001BR", "0002KO"))
  expect_equal(b$suffix, c("BR", "BR", "KO"))
  expect_equal(b$month, as.Date(c("2021-08-01", "2021-09-01", "2020-02-01")))
  expect_equal(b$billed_days, c(31L, 16L, 28L))
  expect_equal(b$month_days, c(31L, 30L, 29L))
  expect_equal(b$kwh_per_day, c(100, 10, 40))

  writeLines(c("NAME,BEGIN,END,KWH_PERIOD", "0001BR,15/01/2021,14/02/2021,900"), path)
  expect_error(load_bills(path), "two calendar months")
})

test_that("log_kwh_matrix leaves NA for a missing month and rejects duplicates", {
  b <- data.frame(site = c("A", "A", "B"),
                  month = as.Date(c("2020-01-01", "2020-03-01", "2020-01-01")),
                  kwh = c(31, 31, 62), kwh_per_day = c(1, 1, 2))
  m <- log_kwh_matrix(b)
  expect_equal(dim(m), c(2L, 3L))
  expect_equal(colnames(m), c("2020-01", "2020-02", "2020-03"))
  expect_true(is.na(m["A", "2020-02"]))
  expect_equal(m["B", "2020-01"], log(2))
  expect_error(log_kwh_matrix(rbind(b, b[1, ])), "same site and month")
  b$kwh[1] <- 0
  expect_error(log_kwh_matrix(b), "Non-positive")
})

test_that("real data: known values and complete 92 x 35 panel", {
  b <- load_bills(root_path(config$data_path))
  expect_equal(nrow(b), 3220)
  expect_false(anyNA(b))

  jan <- b[b$site == "0007BR" & b$month == as.Date("2019-01-01"), ]
  expect_equal(jan$kwh, 2500)
  expect_equal(jan$month_days, 31L)
  expect_equal(jan$billed_days, 31L)

  sep <- b[b$site == "0628BR" & b$month == as.Date("2021-09-01"), ]
  expect_equal(sep$billed_days, 16L)
  expect_equal(sep$month_days, 30L)

  m <- log_kwh_matrix(b)
  expect_equal(dim(m), c(92L, 35L))
  expect_equal(colnames(m)[c(1, 35)], c("2019-01", "2021-11"))
  expect_false(anyNA(m))
  expect_equal(m["0007BR", "2019-01"], log(2500 / 31))
})
