library(deSolve)
library(phaseR)
rm(list = ls()) #réinitialiser tt les var
#x1=Sp
#x2=A
#x3=R
#params=pars <- c(alpha = 0.5, b = 20, c = 0.1, d = 0.01, e = 0.001)

# 
# syst<-function(t,x,parameters){
#   with(as.list(c(parameters,x)),{
#     dSp<- (-gama_*Sp-beta_*(1-ksi)*Sp*A-beta_*ksi*c1/(1+c1)*Sp*Sp+delta*(1-A-Sp)+mu*(1-A-Sp)+mu_etoile*A)
#     dA <-(gama_*Sp+sigma*(1-A-Sp)+beta_*(1-ksi)*Sp*A+beta_*ksi*c1/(1+c1)*Sp*Sp+nu*(1-A-Sp)*A-(zeta+mu_etoile)*A) 
#     #dR <-(zeta*A-nu*R*A-(delta+sigma+mu)*R)
#     #R=(1-A-Sp)
#     list(c(dSp,dA))
#   })
# }

syst<-function(t,x,params){
  x[3]=(1-x[2]-x[1])
  dx1<- (-params[2]*x[1]-params[1]*(1-params[3])*x[1]*x[2]-params[1]*params[3]*params[10]/(1+params[10])*x[1]*x[1]+params[8]*(1-x[2]-x[1])+params[4]*(1-x[2]-x[1])+params[5]*x[2])
  dx2 <-(params[2]*x[1]+params[9]*(1-x[2]-x[1])+params[1]*(1-params[3])*x[1]*x[2]+params[1]*params[3]*params[10]/(1+params[10])*x[1]*x[1]+params[7]*(1-x[2]-x[1])*x[2]-(params[6]+params[5])*x[2]) 
  dx3 <-(params[6]*x[2]-params[7]*x[3]*x[2]-(params[8]+params[9]+params[4])*x[3])
  list(c(dx1,dx2))
  
}

alpha=0.15
epsilon=1
beta=0.0036
ksi=.74
gama=.00744
zeta=1
delta=0.1
nu=.2
sigma=.7
mu=.007288
mu_etoile=0.01155

c1=alpha/(epsilon+gama+mu)
beta_=beta/(1+c1)
gama_=gama/(1+c1)
ksi_=ksi/(1+c1)
c2=c1/(1+c1)
pars<-c(beta_,gama_,ksi_,mu,mu_etoile,zeta,nu,delta,sigma,c1)
times <- seq (0,900, by= 0.01)#seq(0, 20000, length.out = 1000)#"lenght.out" c'est pour le nombre de points, c'est cette 'resolution en sortie' qui compte car si on varie l'amplitude de l'intervalle de temps sur lequel
#on calcule la solution, la taille des pas de temps s'adapte
# set initial conditions
y0 <- c(0.6,0.3)
# solve the ODE system
sol2 <- ode(y = y0, times = times, parms = pars, func = syst,method="ode45")


plot(sol2)
solSp=sol2[,2]
solA=sol2[,3]
#solR=sol2[,4]
tot=solSp+solA
plot(times,tot)
par(mfrow=c(1,1))
plot(solSp~times,col="red",ylim=c(0,1))#,solR,solSp, col=c("red","green","blue"))
#lines(solR,col="green")
lines(solA,col="blue")
lines(tot)


Eq <- findEquilibrium(syst, y0 = c(0.6,0.3),
                     parameters = pars, summary = TRUE)
numericalSolution(syst, y0, parameters=pars,ylim= c(0,1),tlim=c(0,100))
