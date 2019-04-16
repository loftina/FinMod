#' @title Heston Realized Variance and Volitility
#' @description
#' @usage Heston(ticker, start_date, end_date, prediction_period, days_to_drop, kappa, theta, gamma, variance, testing)
#' @param ticker Stock ticker
#' @param start_date Starting date of data to use
#' @param end_date Ending date of data to use
#' @param prediction_period How many days to predict using the model
#' @param days_to_drop Number of days of realized variances to drop when fitting the model
#' @param kappa Mean reverting parameter.
#' @param theta The long term average.
#' @param gamma The volatiltiy of the volatility.
#' @param variance Boolean where true predicts variance and false predicts volitility
#' @param testing Whether to use prediction_period days prior to end_date as the prediction period. If false, the prediction period takes place after end_date
#' @return Realized variance or volitility according to the Heston model when the correct parameter values are used
Heston <- function(ticker, start_date, end_date, prediction_period, days_to_drop, kappa, theta, gamma, variance, testing) {
  stock <- getSymbols(Symbols = ticker, from = start_date, to = end_date, auto.assign = F)
  close_price <- stock[,4]
  log_returns <- diff(log(close_price))[-1]
  time <- 1:length(log_returns)

  realized_variance <- cumsum(log_returns^2)/time
  if (days_to_drop > 0)
  {
    realized_variance <- realized_variance[seq(-1,-days_to_drop)]
  }
  realized_volitility <- sqrt(realized_variance)

  time <- time[1:(length(time)-days_to_drop)]

  if(testing) {
    train_realized_variance <- realized_variance[1:(length(realized_variance)-prediction_period)]
    train_realized_volitility <- realized_volitility[1:(length(realized_volitility)-prediction_period)]
    train_time <- time[1:(length(time)-prediction_period)]
    full_time <- time
  }
  else {
    train_realized_variance <- realized_variance
    train_realized_volitility <- realized_volitility
    train_time <- time
    full_time <- 1:(length(time)+prediction_period)
  }

  nlm.rvar <- nlsLM(formula = train_realized_variance~HestonRealizedVariance(rep(train_realized_variance[1],length(train_time)), train_time, k, t), start = c(k=kappa, t=theta), control = nls.lm.control(maxiter = 1000))

  dates <- rep(time(stock)[nrow(stock)], prediction_period)
  dates <- dates + seq(1,prediction_period)
  dates

  predicted_realized_variance <- HestonRealizedVariance(rep(train_realized_variance[1], length(full_time)), full_time, coef(nlm.rvar)[1], coef(nlm.rvar)[2])

  CIvar <- sqrt(mean((train_realized_variance-predict(nlm.rvar))^2))

  UCvar <- predicted_realized_variance + 1.96*CIvar
  LCvar <- predicted_realized_variance - 1.96*CIvar

  nlm.rvol <- nlsLM(formula = train_realized_volitility~HestonRealizedVolitility(rep(train_realized_variance[1],length(train_time)), train_time, k, t, g), start=c(k=kappa, t=theta, g=gamma))

  predicted_realized_volitility <- HestonRealizedVolitility(rep(train_realized_variance[1], length(full_time)), full_time, coef(nlm.rvol)[1], coef(nlm.rvol)[2], coef(nlm.rvol)[3])

  CIvol <- sqrt(mean((train_realized_volitility-predict(nlm.rvol))^2))

  UCvol <- predicted_realized_volitility + 1.96*CIvol
  LCvol <- predicted_realized_volitility - 1.96*CIvol

  if(variance) {
    plot(time, realized_variance,
         ylim=c(min(LCvar), max(realized_variance)),
         xlim=c(0,length(full_time)),
         main=ticker,
         ylab="Realized Variance",
         xlab="Date")
    lines(full_time, predicted_realized_variance, type='l', col='green')
    lines(full_time, UCvar, col='red')
    lines(full_time, LCvar, col='red')
    abline(v=length(train_time), lty=2, col='blue')
  }
  else {
    plot(time, realized_volitility,
         ylim=c(min(LCvol), max(realized_volitility)),
         xlim=c(0,length(full_time)),
         main=ticker,
         ylab="Realized Volitility",
         xlab="Date")
    lines(full_time, predicted_realized_volitility, type='l', col='green')
    lines(full_time, UCvol, col='red')
    lines(full_time, LCvol, col='red')
    abline(v=length(train_time), lty=2, col='blue')
  }

  if (variance){
    return_data <- data.frame(dates=1:prediction_period,
                              realized_variance=predicted_realized_variance[(length(predicted_realized_variance)-prediction_period+1):length(predicted_realized_variance)],
                              UC95=UCvar[(length(UCvar)-prediction_period+1):length(UCvar)],
                              LC95=LCvar[(length(LCvar)-prediction_period+1):length(LCvar)])
  }
  else {
    return_data <- data.frame(dates=1:prediction_period,
                              realized_volitility=predicted_realized_volitility[(length(predicted_realized_volitility)-prediction_period+1):length(predicted_realized_volitility)],
                              UC95=UCvol[(length(UCvol)-prediction_period+1):length(UCvol)],
                              LC95=LCvol[(length(LCvol)-prediction_period+1):length(LCvol)])
  }

  return(return_data)
}

#' @title Heston Realized Variance
#' @description The Heston model's approximation for realized variance.
#' @usage HestonRealizedVariance(realizedVar_0, time, kappa, theta)
#' @param realizedVar_0 Initial realized variance of a stock repeated T times
#' @param time Number of realized variance observations
#' @param kappa Mean reverting parameter.
#' @param theta The long term average.
#' @return Realized variance according to the Heston model when the correct kappa and theta values are used
HestonRealizedVariance <- function(realizedVar_0, time, kappa, theta) {
  return ( ((1-exp(-kappa*time))/(kappa*time))*(realizedVar_0-theta^2)+theta^2 )
}


#' @title Heston Variance of Realized Variance
#' @description The variance of the realized variance of a stock according to
#' the Heston model.
#' @usage HestonVariance(realizedVar_0, time, kappa, theta, gamma)
#' @param realizedVar_0 Initial realized variance of a stock repeated T times
#' @param time Number of realized variance observations
#' @param kappa Mean reverting parameter.
#' @param theta The long term average.
#' @param gamma The volatiltiy of the volatility.
#' @return The variance of the realized variance according to the Heston model when correct kappa, theta, and gamma values are used
HestonVariance <- function(realizedVar_0, time, kappa, theta, gamma) {
  return ( (gamma^2*exp(-2*kappa*time))/(2*kappa^3*time^2)*((2*exp(2*kappa*time)-4*exp(kappa*time)*kappa*time-2)*(realizedVar_0^2-theta^2)+(2*exp(2*kappa*time)*kappa*time-3*exp(2*kappa*time)+4*exp(kappa*time)-1)*theta^2) )
}


#' @title Heston Realized Volitility
#' @description The Heston model's approximation for realized volitility.
#' @usage HestonRealizedVolitility(realizedVar_0, time, kappa, theta, gamma)
#' @param realizedVar_0 Initial realized variance of a stock repeated T times
#' @param time Number of realized variance observations
#' @param kappa Mean reverting parameter.
#' @param theta The long term average.
#' @param gamma The volatiltiy of the volatility.
#' @return Realized volitility according to the Heston model when the correct kappa, theta, and gamma values are used
HestonRealizedVolitility <- function (realizedVar_0, time, kappa, theta, gamma) {
  return ( sqrt(HestonRealizedVariance(realizedVar_0, time, kappa, theta))-HestonVariance(realizedVar_0, time, kappa, theta, gamma)/(8*((realizedVar_0/(kappa*time))*(exp(kappa*time)-1))^(3/2)) )
}
