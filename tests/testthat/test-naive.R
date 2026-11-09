test_that("naive_per_day_forecast scales the last daily rate to the target month", {
  b <- data.frame(site = c("A", "A", "B", "B"), suffix = c("BR", "BR", "KO", "KO"),
                  month = as.Date(c("2021-10-01", "2021-11-01", "2021-10-01", "2021-11-01")),
                  kwh_per_day = c(5, 10, 7, 2))
  f <- naive_per_day_forecast(b, as.Date("2021-12-01"))
  expect_equal(f$site, c("A", "B"))
  expect_equal(f$kwh, c(310, 62))
  f2 <- naive_per_day_forecast(b[b$month < as.Date("2021-11-01"), ], as.Date("2021-11-01"))
  expect_equal(f2$kwh, c(150, 210))
  expect_error(naive_per_day_forecast(b, as.Date("2022-01-01")), "month before")
})

test_that("real data: 92 December forecasts without missing values", {
  b <- load_bills(root_path(config$data_path))
  f <- naive_per_day_forecast(b, config$forecast_month)
  expect_equal(nrow(f), 92)
  expect_false(anyNA(f$kwh))
  nov <- b[b$site == "0007BR" & b$month == as.Date("2021-11-01"), ]
  expect_equal(f$kwh[f$site == "0007BR"], nov$kwh / 30 * 31)
})
