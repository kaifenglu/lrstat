library(lrstat)

# input parameters
w = c(0.2, 0.8, 0, 0, 0, 0)
G = matrix(c(0, 1, 0, 0, 0, 0,  0, 0, 0.2, 0.8, 0, 0, 
             0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0.2, 0.8,
             0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 1, 0), 
           nrow=6, byrow=TRUE)

alpha = 0.025
asf = rep("sfOF", 6)
incid = matrix(c(1, 1, 0, rep(1, 15)), nrow=6, ncol=3, byrow=TRUE)
maxinfo = c(259, 183, 321, 195, 392, 230)

info = matrix(c(
  186, 259, 299, 
  96, 148, 183, 
  225, 296, 324, 
  100, 161, 195, 
  245, 341, 381, 
  122, 177, 210), 
  nrow=6, ncol=3, byrow=TRUE)

p = matrix(c(
  1.813147e-06, 1.559218e-07, 5.532412e-09,
  4.940081e-03, 7.710547e-05, 6.456224e-07,
  2.269443e-02, 6.872631e-03, 6.002544e-04,
  5.849374e-02, 2.413134e-03, 4.860418e-05,
  9.993609e-04, 5.871850e-04, 1.870486e-03,
  1.129593e-01, 1.282247e-01, 1.384501e-01), 
  nrow=6, ncol=3, byrow=TRUE)


(reject <- fseqbon(w=w, G=G, alpha=alpha, asf=asf, incid=incid, 
                   maxinfo=maxinfo, info=info, p=p))


# look 1, Itau = {1, 2, 3, 4, 5, 6} 
I = seq(1,6)

# nominal significance levels for H1
u = getBound(
  kMax = 2, informationRates = c(info[1,1]/maxinfo[1], 1), 
  alpha = w[1]*alpha, typeAlphaSpending = asf[1])

(astar = 1 - pnorm(u[1]))

# since p[1,1] < astar, reject H1, update w and G
(newgraph = updateGraph(w, G, I, 1))
w = newgraph$w
G = newgraph$G
I = newgraph$I

# nominal significance levels for H2
u = getBound(
  kMax = 2, informationRates = c(info[2,1]/maxinfo[2], 1), 
  alpha = w[2]*alpha, typeAlphaSpending = asf[2])
(astar = 1 - pnorm(u[1]))

# since p[2,1] > astar, cannot reject H2
# the remaining hypotheses have zero weights and cannot be rejected in look 1
# go to look 2, Itau = {2, 3, 4, 5, 6}

# nominal significance level for H2
u = getBound(
  kMax = 3, informationRates = c(info[2,1:2]/maxinfo[2], 1), 
  alpha = w[2]*alpha, typeAlphaSpending = asf[2])
(astar = 1 - pnorm(u[2]))

# since p[2,2] < astar, reject H2, update w and G
(newgraph = updateGraph(w, G, I, 2))
w = newgraph$w
G = newgraph$G
I = newgraph$I

# nominal significance level for H3
u = getBound(
  kMax = 3, informationRates = c(info[3,1:2]/maxinfo[3], 1), 
  alpha = w[3]*alpha, typeAlphaSpending = asf[3])
(astar = 1 - pnorm(u[2]))

# since p[3,2] > astar, cannot reject H3

# nominal significance level for H4
u = getBound(
  kMax = 3, informationRates = c(info[4,1:2]/maxinfo[4], 1), 
  alpha = w[4]*alpha, typeAlphaSpending = asf[4])
(astar = 1 - pnorm(u[2]))

# since p[4,2] < astar, reject H4, update w and G
(newgraph = updateGraph(w, G, I, 4))
w = newgraph$w
G = newgraph$G
I = newgraph$I

# since no change in significance level for H3, H3 cannot be rejected

# nominal significance level for H5
u = getBound(
  kMax = 3, informationRates = c(info[5,1:2]/maxinfo[5], 1), 
  alpha = w[5]*alpha, typeAlphaSpending = asf[5])
(astar = 1 - pnorm(u[2]))

# since p[5,2] < astar, reject H5, update w and G
(newgraph = updateGraph(w, G, I, 5))
w = newgraph$w
G = newgraph$G
I = newgraph$I

# since no change in significance level for H3, H3 cannot be rejected

# nominal significance level for H6
u = getBound(
  kMax = 3, informationRates = c(info[6,1:2]/maxinfo[6], 1), 
  alpha = w[6]*alpha, typeAlphaSpending = asf[6])
(astar = 1 - pnorm(u[2]))

# since p[6,2] > astar, H6 cannot be rejected, go to look 3, Itau = {3, 6}

# nominal significance level for H3, note the observed maximum information
# is different from the planned maximum information, need to adjust
# the critical value at the last look
u1 = getBound(
  kMax = 3, informationRates = c(info[3,1:2]/maxinfo[3], 1), 
  alpha = w[3]*alpha, typeAlphaSpending = asf[3])

probs = exitprob(b = u1, a = c(rep(-6,2), u1[3]), theta = 0, I = info[3,])

u = getBound(
  kMax = 3, informationRates = c(info[3,1:3]/info[3,3]),
  alpha = w[3]*alpha, typeAlphaSpending = "user",
  userAlphaSpending = c(cumsum(probs$exitProbUpper)[1:2], w[3]*alpha))
(astar = 1 - pnorm(u[3]))

# since p[3,3] < astar, H3 is rejected, update w and G
(newgraph = updateGraph(w, G, I, 3))
w = newgraph$w
G = newgraph$G
I = newgraph$I

# nominal significance level for H6, note the observed maximum information
# is different from the planned maximum information, need to adjust
# the critical value at the last look
u1 = getBound(
  kMax = 3, informationRates = c(info[6,1:2]/maxinfo[6], 1), 
  alpha = w[6]*alpha, typeAlphaSpending = asf[6])

probs = exitprob(b = u1, a = c(rep(-6,2), u1[3]), theta = 0, I = info[6,])

u = getBound(
  kMax = 3, informationRates = c(info[6,1:3]/info[6,3]),
  alpha = w[6]*alpha, typeAlphaSpending = "user",
  userAlphaSpending = c(cumsum(probs$exitProbUpper)[1:2], w[6]*alpha))
(astar = 1 - pnorm(u[3]))

# since p[6,3] > astar, H6 cannot be rejected, the procedure stops
# the final decision is to reject H1 at look 1, reject H2, H4, and H5 
# at look 2, reject H3 at look 3, cannot reject H6
