
coef_specs <- list(
  '1' = list(
    a1 = list(type = "norm", mean = -5, sd = 1),
    a2 = list(type = "norm", mean = 7, sd = 0.5)
  ),
  '2' = list(
    b1 = list(type = "unif", min = 4, max = 6),
    b2 = list(type = "norm", mean = 0, sd = 0.5)
  ),
  '3' = list(
    c1 = list(type = "norm", mean = -3, sd = 0.6),
    c2 = list(type = "norm", mean = 2, sd = 0.25),
    c3 = list(type = "norm", mean = -2, sd = 0.5),
    c4 = list(type = "norm", mean = 2, sd = 0.75)
  ),
  '4' = list(
    d1 = list(type = "norm", mean = -2, sd = 0.5),
    d2 = list(type = "norm", mean = 3, sd = sqrt(0.75^2))
  ),
  '5' = list(
    e1 = list(type = "unif", min = 4, max = 5),
    e2 = list(type = "norm", mean = 2, sd = sqrt(0.2^2))
  ),
  '6' = list(
    f1 = list(type = "norm", mean = 4, sd = sqrt(1^2)),
    f2 = list(type = "norm", mean = -3, sd = sqrt(0.25^2)),
    f3 = list(type = "norm", mean = 1, sd = 0.5)
  )
)

