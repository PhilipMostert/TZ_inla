beta <- seq(-20, 20, len = 500)
prec <- 0.01
density <- dnorm(beta, 0, 1/prec)

fixed_covariate <- c(scale(seq(-50, 50, len = 100)))
lin_prob <- outer(beta, fixed_covariate, '*')
prob <- inla.link.cloglog(lin_prob, inverse = TRUE)

plot(fixed_covariate, prob[1,], type = "l", col = "blue", ylim = c(0,1),
     xlab = 'Scaled covariate value', ylab = 'Probability')
lines(fixed_covariate, prob[125,], col = "green", ylim = c(0,1))
lines(fixed_covariate, prob[200,], col = "darkgreen", ylim = c(0,1))
lines(fixed_covariate, prob[250,], col = "black")
lines(fixed_covariate, prob[300,], col = "darkorange")
lines(fixed_covariate, prob[375,], col = "orange")
lines(fixed_covariate, prob[500,], col = "red")
legend("right", lty = 1, 
       col = c('blue',  'green', 'darkgreen', 'black', 'darkorange','orange', 'red'),
       legend = c('beta= -20', 
                  'beta= -10',
                  'beta = -4',
                  'beta = 0',
                  'beta = 4',
                  'beta = 10',
                  'beta = 20'))

plot(beta, dnorm(beta, 0, 1/0.001), type= 'l', ylim = c(0, 0.5), col = 'black')
lines(beta, dnorm(beta, 0, 1/0.01),  col = 'orange')

lines(beta, dnorm(beta, 0, 1/0.1),  col = 'blue')
lines(beta, dnorm(beta, 0, 1/1), col = 'green')
legend('topright', lty = 1, col = c('black', 'orange', 'blue', 'green'),
       legend = c('prec = 0.001', 'prec = 0.01',
                  'prec = 0.1', 'prec = 1'))
