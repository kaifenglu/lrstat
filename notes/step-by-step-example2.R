library(lrstat)

eps = 1.0e-8

I = seq(1,6) # six hypotheses

# input parameters
w = c(0.1, 0.4, 0.1, 0.4, 0, 0) 

G = matrix(c(0, 1, 0, 0, 0, 0,  
             0, 0, 0.2, 0.8, 0, 0, 
             0, 0, 0, 1, 0, 0,  
             0.2*(1-eps), 0.8*(1-eps), 0, 0, 0.2*eps, 0.8*eps,
             0, 0, 0, 0, 0, 1,  
             0, 0, 0, 0, 1, 0), 
           nrow=6, byrow=TRUE) 

alpha = 0.025

asf = rep("sfHSD", 6)

asfpar = rep(-4, 6)

incid = matrix(c(1, 1, 0, 
                 1, 1, 1,
                 1, 1, 0,
                 1, 1, 1,
                 1, 1, 1,
                 1, 1, 1), 
               nrow=6, ncol=3, byrow=TRUE)

maxinfo = c(259, 183, 274, 195, 392, 230)
info = matrix(c(186, 259, 295,
                103, 156, 183,
                190, 289, 330,
                101, 163, 197,
                222, 340, 393,
                114, 189, 232), 
              nrow=6, ncol=3, byrow=TRUE)

p = matrix(c(0.0007349585, 2.929598e-05, 5.107926e-06,
             0.0046922022, 1.708761e-03, 1.673969e-03,
             0.0039662778, 2.113065e-02, 2.329271e-02,
             0.0031230745, 7.175618e-03, 1.950093e-02,
             0.2743761599, 8.778801e-03, 9.652760e-04,
             0.5640949238, 3.029222e-01, 1.429600e-01),
           nrow=6, ncol=3, byrow=TRUE)


(reject <- fseqbon(w=w, G=G, alpha=alpha, asf=asf, asfpar=asfpar, 
                   incid=incid, maxinfo=maxinfo, info=info, p=p))


# look 1, Itau = {1, 2, 3, 4, 5, 6} 

# nominal significance levels for H1
u = getBound(
  kMax = 2, informationRates = c(info[1,1]/maxinfo[1], 1), 
  alpha = w[1]*alpha, typeAlphaSpending = asf[1], 
  parameterAlphaSpending = asfpar[1])

(astar = 1 - pnorm(u[1]))

# since p[1,1] < astar, reject H1, update w and G
(newgraph = updateGraph(w, G, I, 1))
w = newgraph$w
G = newgraph$G
I = newgraph$I

# nominal significance levels for H2
u = getBound(
  kMax = 2, informationRates = c(info[2,1]/maxinfo[2], 1), 
  alpha = w[2]*alpha, typeAlphaSpending = asf[2], 
  parameterAlphaSpending = asfpar[2])
(astar = 1 - pnorm(u[1]))

# since p[2,1] > astar, cannot reject H2

# nominal significance levels for H3
u = getBound(
  kMax = 2, informationRates = c(info[3,1]/maxinfo[3], 1), 
  alpha = w[3]*alpha, typeAlphaSpending = asf[3], 
  parameterAlphaSpending = asfpar[3])
(astar = 1 - pnorm(u[1]))

# since p[3,1] > astar, cannot reject H3

# nominal significance levels for H4
u = getBound(
  kMax = 2, informationRates = c(info[4,1]/maxinfo[4], 1), 
  alpha = w[4]*alpha, typeAlphaSpending = asf[4], 
  parameterAlphaSpending = asfpar[4])
(astar = 1 - pnorm(u[1]))

# since p[4,1] > astar, cannot reject H4

# the remaining hypotheses have zero weights and cannot be rejected in look 1
# go to look 2, Itau = {2, 3, 4, 5, 6}


# nominal significance level for H2
u = getBound(
  kMax = 3, informationRates = c(info[2,1:2]/maxinfo[2], 1), 
  alpha = w[2]*alpha, typeAlphaSpending = asf[2], 
  parameterAlphaSpending = asfpar[2])
(astar = 1 - pnorm(u[2]))

# since p[2,2] < astar, reject H2, update w and G
(newgraph = updateGraph(w, G, I, 2))
w = newgraph$w
G = newgraph$G
I = newgraph$I

# nominal significance level for H3
u = getBound(
  kMax = 2, informationRates = c(info[3,1]/maxinfo[3], 1), 
  alpha = w[3]*alpha, typeAlphaSpending = asf[3], 
  parameterAlphaSpending = asfpar[3])

# update bound at look 2 (last look for H3)
u = getBound(
  kMax = 2, informationRates = info[3,1:2]/info[3,2], 
  alpha = w[3]*alpha, typeAlphaSpending = "user",
  userAlphaSpending = c(1-pnorm(u[1]), w[3]*alpha))
(astar = 1 - pnorm(u[2]))

# since p[3,2] > astar, cannot reject H3

# nominal significance level for H4
u = getBound(
  kMax = 3, informationRates = c(info[4,1:2]/maxinfo[4], 1), 
  alpha = w[4]*alpha, typeAlphaSpending = asf[4], 
  parameterAlphaSpending = asfpar[4])
(astar = 1 - pnorm(u[2]))

# since p[4,2] < astar, reject H4, update w and G
(newgraph = updateGraph(w, G, I, 4))
w = newgraph$w
G = newgraph$G
I = newgraph$I

# nominal significance level for H3
u = getBound(
  kMax = 2, informationRates = c(info[3,1]/maxinfo[3], 1), 
  alpha = w[3]*alpha, typeAlphaSpending = asf[3], 
  parameterAlphaSpending = asfpar[3])

# update bound at look 2 (last look for H3)
u = getBound(
  kMax = 2, informationRates = info[3,1:2]/info[3,2], 
  alpha = w[3]*alpha, typeAlphaSpending = "user",
  userAlphaSpending = c(1-pnorm(u[1]), w[3]*alpha))
(astar = 1 - pnorm(u[2]))

# since p[3,2] < astar, reject H3, update w and G
(newgraph = updateGraph(w, G, I, 3))
w = newgraph$w
G = newgraph$G
I = newgraph$I



# nominal significance level for H5
u = getBound(
  kMax = 3, informationRates = c(info[5,1:2]/maxinfo[5], 1), 
  alpha = w[5]*alpha, typeAlphaSpending = asf[5],
  parameterAlphaSpending = asfpar[5])
(astar = 1 - pnorm(u[2]))

# since p[5,2] > astar, H5 cannot be rejected

# nominal significance level for H6
u = getBound(
  kMax = 3, informationRates = c(info[6,1:2]/maxinfo[6], 1), 
  alpha = w[6]*alpha, typeAlphaSpending = asf[6],
  parameterAlphaSpending = asfpar[6])
(astar = 1 - pnorm(u[2]))

# since p[6,2] > astar, H6 cannot be rejected, go to look 3

# nominal significance level for H5, note the observed maximum information
# is slightly different from the planned maximum information, need to adjust
# the critical value at the last look
u = getBound(
  kMax = 3, informationRates = c(info[5,1:2]/maxinfo[5], 1), 
  alpha = w[5]*alpha, typeAlphaSpending = asf[5], 
  parameterAlphaSpending = asfpar[5])

probs = exitprob(b = u, a = c(rep(-6,2), u[3]), theta = 0, I = info[5,])

u = getBound(
  kMax = 3, informationRates = c(info[5,1:3]/info[5,3]),
  alpha = w[5]*alpha, typeAlphaSpending = "user",
  userAlphaSpending = c(cumsum(probs$exitProbUpper)[1:2], w[5]*alpha))
(astar = 1 - pnorm(u[3]))

# since p[5,3] < astar, H5 is rejected, update w and G
(newgraph = updateGraph(w, G, I, 5))
w = newgraph$w
G = newgraph$G
I = newgraph$I

# nominal significance level for H6, note the observed maximum information
# is different from the planned maximum information, need to adjust
# the critical value at the last look
u = getBound(
  kMax = 3, informationRates = c(info[6,1:2]/maxinfo[6], 1), 
  alpha = w[6]*alpha, typeAlphaSpending = asf[6], 
  parameterAlphaSpending = asfpar[6])

probs = exitprob(b = u, a = c(rep(-6,2), u[3]), theta = 0, I = info[6,])

u = getBound(
  kMax = 3, informationRates = c(info[6,1:3]/info[6,3]),
  alpha = w[6]*alpha, typeAlphaSpending = "user",
  userAlphaSpending = c(cumsum(probs$exitProbUpper)[1:2], w[6]*alpha))
(astar = 1 - pnorm(u[3]))

# since p[6,3] > astar, H6 cannot be rejected, the procedure stops
# the final decision is to reject H1 at look 1, reject H2, H4, and H3 
# at look 2, reject H5 at look 3, cannot reject H6
