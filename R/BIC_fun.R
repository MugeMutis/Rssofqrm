BIC_fun = function(y, yfit, tau, h){
  n = length(y)
  
  res = y - yfit
  
  BIC_val = sum(check_loss(res, tau)) + (h * log(n))
  
  return(BIC_val)
}