set.seed(76137)

# Set dimensions
n_list = c(10^2, 10^3)
err_list = c("gamma", "lognormal")
link_list = c("linear", "cubic", "trigonometric")
settings = expand.grid(n = n_list, err = err_list, link = link_list, stringsAsFactors = F)
n_settings = nrow(settings)
n_rep = 100
data = vector(mode = "list", length = n_settings)
names(data) = do.call(paste, c(settings, sep = "-"))

# Utils
draw_err = function(n, distr){
  switch(distr,
         gamma = rgamma(n, 2, 2),
         lognormal = exp(rnorm(n, 0, 1)))
}

f = function(x, link){
  switch(link,
         linear = x,
         cubic = x^3,
         trigonometric = sin(2*(4*x-2) + 2*exp(-16^2*(x-0.5)^2))
         )
}

# Log-normal error
# For each setting draw covariate and response
# I fix covariate distribution as U(-2, 2) as in Semiparametric Mode Regression section 3.3
for (i in 1:nrow(settings)){
  n = settings[i, "n"]
  err = settings[i, "err"]
  link = settings[i, "link"]
  
  tmp = array(NA, dim = c(n_rep, n, 2), 
              dimnames = list("rep"  = 1:n_rep,
                              "unit" = 1:n,
                              "var"  = c("y", "x")))
  for (b in 1:n_rep){
    tmp_x = runif(n = n, min = -2, max = 2)
    tmp_err = draw_err(n, err)
    tmp_y = f(tmp_x, link) + tmp_err
    tmp[b, , "y"] = tmp_y
    tmp[b, , "x"] = tmp_x
  }
  
  data[[i]] = tmp
}

