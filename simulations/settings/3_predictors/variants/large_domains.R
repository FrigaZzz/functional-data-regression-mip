time_domains = list(
  list(-10, 10),
  list(-10, 10),
  list(-10, 10)
)
mu_funcs <- list(
  function(t) sin(2 * pi * t),
  function(t) cos(3 * pi * t),
  function(t) cos(5 * pi * t)
)
