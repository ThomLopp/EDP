function dydt = osc(t,y)
    %constantes
    alpha=0.15;
    epsilon=1;
    beta=0.0036;
    ksi=.74;
    gama=.00744;
    zeta=1;
    delta=0.1;
    nu=.2;
    sigma=.7;
    mu=.007288;
    mu_=0.01155;
    

%y est le vecteur de solutions y = S,P,A,R
    dydt=zeros(4,1);
        
    dydt(1)=-alpha*sigma-beta*(1-ksi)*y(1)*y(3)-beta*ksi*y(1)*y(2)+epsilon*y(2) +delta*y(4)+mu*(y(2)+y(4))+mu_*y(3);
    dydt(2)=alpha *y(1) -(epsilon+gama+mu)*y(2);
    dydt(3)=gama*y(1)+ sigma*y(4)+beta*(1-ksi)*y(1)*y(3)+beta*ksi*y(1)*y(2)+nu*y(4)*y(3)-(zeta+mu_)*y(3);
    dydt(4)=zeta*y(3)-nu*y(4)*y(3)-(delta+sigma+mu)*y(4);
end
    