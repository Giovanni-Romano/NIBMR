mean_IG = function(a, b){
  if (a > 1){
    b / (a-1)
  } else {
    Inf
  }
}

var_IG = function(a, b){
  if (a > 2){
    b / (a-2) / (a-1)
  } else {
    Inf
  }
}

p_tail <- function(q, a, b) {
  pgamma(1 / q, shape = a, rate = b)
}

find_b = function(a, q, alpha){
  qgamma(alpha, shape = a, rate = 1) * q
}

q_val = 5^2
thr_val = 0.001
a_val = seq(1.1, 5, by = 0.1)
b_val = sapply(a_val, find_b, q = q_val, alpha = thr_val)
plot(a_val, b_val, type= "b", main = "b as function of a")

par_val = cbind(a_val, b_val)
mean_val = apply(par_val, 1, function(p) mean_IG(p[1], p[2]))
var_val = apply(par_val, 1, function(p) var_IG(p[1], p[2]))
tail_val = apply(par_val, 1, function(p) p_tail(q_val, p[1], p[2]))

cbind(a = a_val, b = b_val, mean = mean_val, var = var_val, tail = tail_val)
