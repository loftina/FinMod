#' @title Barndorff‐Nielsen and Shephard Realized Variance
#' @description
#' @usage BNS(ticker, start_date, end_date, prediction_period, days_to_drop, lambda1, delta1, gamma1, phi1, testing)
#' @param ticker Stock ticker
#' @param start_date Starting date of data to use
#' @param end_date Ending date of data to use
#' @param prediction_period How many days to predict using the model
#' @param days_to_drop Number of days of realized variances to drop when fitting the model
#' @param lambda The mean reverting variable
#' @param delta Influences the tail thickness of the disribution
#' @param gamma Influences the scale of the distribution
#' @param phi Volatility of the volatility.
#' @param testing Whether to use prediction_period days prior to end_date as the prediction period. If false, the prediction period takes place after end_date
#' @return Realized variance according to the Barndorff‐Nielsen and Shephard model when the correct parameter values are used
BNS <- function(ticker, start_date, end_date, prediction_period, days_to_drop, lambda, delta, gamma, phi, testing) {
  stock <- getSymbols(Symbols = ticker, from = start_date, to = end_date, auto.assign = F)
  close_price <- stock[,4]
  log_returns <- diff(log(close_price))[-1]
  time <- 1:length(log_returns)

  realized_variance <- cumsum(log_returns^2)/time
  if (days_to_drop > 0)
  {
    realized_variance <- realized_variance[seq(-1, -days_to_drop)]
  }
  time <- time[1:(length(time)-days_to_drop)]

  if(testing) {
    train_realized_variance <- realized_variance[1:(length(realized_variance)-prediction_period)]
    train_time <- time[1:(length(time)-prediction_period)]
    full_time <- time
  }
  else {
    train_realized_variance <- realized_variance
    train_time <- time
    full_time <- 1:(length(time)+prediction_period)
  }

  nlm.rvar <- nlsLM(formula = train_realized_variance~BNSRealizedVariance(rep(train_realized_variance[1],length(train_time)), train_time, l, d, g, phi), start = c(l=lambda, d=delta, g=gamma), control = nls.lm.control(maxiter = 1000))

  predicted_realized_variance <- BNSRealizedVariance(rep(train_realized_variance[1], length(full_time)), full_time, coef(nlm.rvar)[1], coef(nlm.rvar)[2], coef(nlm.rvar)[3], phi)

  CIvar <- sqrt(mean((train_realized_variance-predict(nlm.rvar))^2))

  UCvar <- predicted_realized_variance + 1.96*CIvar
  LCvar <- predicted_realized_variance - 1.96*CIvar

  plot(time, realized_variance,
       ylim=c(min(LCvar),max(realized_variance)),
       xlim=c(0,length(full_time)),
       main=ticker,
       ylab="Realized Variance",
       xlab="Date")
  lines(full_time, predicted_realized_variance, type = 'l', col='green')
  lines(full_time, UCvar, col='red')
  lines(full_time, LCvar, col='red')
  abline(v=length(train_time), lty=2, col='blue')

  return_data <- list(
    'realized variance' = data.frame(dates=1:prediction_period,
                                     realized_variance=predicted_realized_variance[(length(predicted_realized_variance)-prediction_period+1):length(predicted_realized_variance)],
                                     UC95=UCvar[(length(UCvar)-prediction_period+1):length(UCvar)],
                                     LC95=LCvar[(length(LCvar)-prediction_period+1):length(LCvar)])
  )

  return(return_data)
}

#' @title Barndorff‐Nielsen and Shephard Realized Variance
#' @description The Barndorff‐Nielsen and Shephard model's approximation for realized variance.
#' @usage BNSRealizedVariance(realizedVar_0, time, lambda, delta, gamma, phi)
#' @param realizedVar_0 Initial realized variance of a stock repeated T times
#' @param time Number of realized variance observations
#' @param lambda The mean reverting variable
#' @param delta Influences the tail thickness of the disribution
#' @param gamma Influences the scale of the distribution
#' @param phi Volatility of the volatility.
#' @return Realized variance according to the Barndorff‐Nielsen and Shephard model when the correct kappa value is used
BNSRealizedVariance <- function(realizedVar_0, time, lambda, delta, gamma, phi) {
  return ( (1/time)*(lambda^(-1)*(1-exp(-lambda*time))*realizedVar_0 + ExpectedZ1(delta, gamma)*(time-lambda^(-1)*(1-exp(-lambda*time)))) + phi^2*lambda*VarianceZ1(delta, gamma) )
}

#' @title
#' @description
#' @usage ExpectedZ1(delta, gamma)
#' @param delta Influences the tail thickness of the disribution
#' @param gamma Influences the scale of the distribution
#' @return Expected value of Z
ExpectedZ1 <- function (delta, gamma) {
  return ( delta*gamma^(-1) )
}

#' @title
#' @description
#' @usage VarianceZ1(delta, gamma)
#' @param delta Influences the tail thickness of the disribution
#' @param gamma Influences the scale of the distribution
#' @return Variance of Z
VarianceZ1 <- function (delta, gamma) {
  return ( 2*delta*gamma^(-3) )
}
