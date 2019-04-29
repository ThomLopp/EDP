%
% dS = -alpha*sigma-beta*(1-ksi)*S*A-beta*ksi*S*P+epsilon*P +delta*R+mu*(P+R)+mu_*A;
% dP = alpha *S -(epsilon+gama+mu)*P;
% dA = gama*P+ sigma*R+beta*(1-ksi)*S*A+beta*ksi*S*P+nu*R*A-(zeta+mu_)*A;
% dR = zeta*A-nu*R*A-(delta+sigma+mu)*R;
tspan = [0 3000]; 
S0 = 100000;    
P0 = 0;
A0 = 0;
R0 = 0;
[T,Y] = ode15s(@osc,tspan,[S0,P0,A0,R0]); 
plot(T,Y(:,1),'o')
