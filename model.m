%% Turing 2D: FitzHugh-Nagumo 

% parametre des equations
Dv = 0.001;%DS
Dw = 0.0010; %DA
Dr = 0.00;%DR

% parametres de simulation, espace
% parametre des equations
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
    mu=0.007288;
    mu_etoile=0.01155;
       
    c1=alpha/(epsilon+gama+mu);
    beta_=beta/(1+c1);
    gama_=gama/(1+c1);
    ksi_=ksi/(1+c1);
    c2=c1/(1+c1);
    
s = 0.6; % taille domaine (carre sxs)
h = 0.2;
x0 = 0;
x1 = s;
y0 = 0;
y1 = s;
x = x0:h:x1;
y = y0:h:y1;
[X,Y] = meshgrid(x,y);
J = length(x);
J2 = J*J;

% variable dynamiques

Sp = zeros(J2,1); % stocke seulement l'etat au temps t
Sp_ =reshape(Sp,J,J);
A = zeros(J2,1);
A_ = reshape(A,J,J);
R = zeros(J2,1);
newSp = zeros(J2,1);
newA = zeros(J2,1);
newR = zeros(J2,1);

% Condition periodiques
L = sparse(1:J2,1:J2,-4); % matrice creuse, compacte en memoire
coinhautgauche = 1;
coinbasgauche = J;
coinhautdroit = J*(J-1)+1;
coinbasdroit = J2;
bordgauche = 2:J-1;
bordhaut = J+1:J:J*(J-2)+1;
bordbas = 2*J:J:J*(J-1);
borddroit = J*(J-1)+2:J2-1;
bord = [coinhautgauche, coinhautdroit, coinbasgauche, coinbasdroit, ...
    bordgauche, bordhaut, bordbas, borddroit];
 
interieur = setdiff(1:J2, bord);


%conditions au bord
C=zeros(J2,1);
C=C + sparse(1:J,1,1,J2,1);
C=C + sparse((J2-J):J2,1,1,J2,1);
C(1:J:end)=1;
C((J-1):J:end)=1;
C((J2-J):1:end)=1;
%C = sparse(1:J2,1:J2,1);


figure(1); clf;
surf(X,Y,reshape(C,J,J),'EdgeColor','none');
view(2)
drawnow;

pause;

% interieur
L = L + sparse(1:(J2-1),(1:(J2-1))+1,1,J2,J2);%diagonal sup
L = L + sparse(2:(J2),(2:(J2))-1,1,J2,J2);%diagonal inf
L = L + sparse(1:(J2-J),(1:(J2-J))+J,1,J2,J2);%diagonal sup +J
L = L + sparse((J+1):J2,((J+1):J2)-J,1,J2,J2);%diagonal sup -J
 
%on veut que les bords ne soient pas prise en compte par le Laplacien
L = L + sparse(bordhaut,bordhaut+1,0,J2,J2);
L = L + sparse(bordhaut,bordhaut+J-1,0,J2,J2);
L = L + sparse(bordhaut,bordhaut+J,0,J2,J2);
L = L + sparse(bordhaut,bordhaut-J,0,J2,J2);


L = L + sparse(bordgauche,bordgauche+1,0,J2,J2);
L = L + sparse(bordgauche,bordgauche-1,0,J2,J2);
L = L + sparse(bordgauche,bordgauche+J,0,J2,J2);
L = L + sparse(bordgauche,bordgauche+J*(J-1),0,J2,J2);

L = L + sparse(bordbas,bordbas-(J-1),0,J2,J2);
L = L + sparse(bordbas,bordbas-1,0,J2,J2);
L = L + sparse(bordbas,bordbas+J,0,J2,J2);
L = L + sparse(bordbas,bordbas-J,0,J2,J2);

L = L + sparse(borddroit,borddroit+1,0,J2,J2);
L = L + sparse(borddroit,borddroit-1,0,J2,J2);
L = L + sparse(borddroit,borddroit-J*(J-1),0,J2,J2);
L = L + sparse(borddroit,borddroit-J,0,J2,J2);

% coins
L(coinhautgauche,coinhautgauche+1) = 1;
L(coinhautgauche,coinhautgauche+J-1) = 1;
L(coinhautgauche,coinhautgauche+J) = 1;
L(coinhautgauche,coinhautgauche+J*(J-1)) = 1;


L(coinbasgauche,coinbasgauche-(J-1)) = 1;
L(coinbasgauche,coinbasgauche-1) = 1;
L(coinbasgauche,coinbasgauche+J) = 1;
L(coinbasgauche,coinbasgauche+J*(J-1)) = 1;

L(coinhautdroit,coinhautdroit+1) = 1;
L(coinhautdroit,coinhautdroit+J-1) = 1;
L(coinhautdroit,coinhautdroit-J*(J-1)) = 1;
L(coinhautdroit,coinhautdroit-J) = 1;

L(coinbasdroit,coinbasdroit-(J-1)) = 1;
L(coinbasdroit,coinbasdroit-1) = 1;
L(coinbasdroit,coinbasdroit-J*(J-1)) = 1;
L(coinbasdroit,coinbasdroit-J) = 1;

% condition initiales
Sp_(1,1)= 0.5;

A_(1,1)=0.5;
% Sp_(141:251,1:100)=0.1;
% Sp_(1:100,141:251)=0.1;
% Sp_(141:251,141:251)=0.1;
Sp=reshape(Sp_,J2,1);
A=reshape(A_,J2,1);
%Sp =   0.1*(rand(J^2,1)); (1:100,1:100)
%A = 0 + 0.1*(-0 + 0*rand(J^2,1));
%R =  0.0 + 0.1*(-0 + 0*rand(J^2,1));

% parametres de simulation, temps
t0 = 0;
tfinal = 20000000; 
t = t0;
% k doit etre < a h^2/2/d/D ou dim d = 2
k = min(1,0.2*( h^2/4/max(Dv,Dw) ));
 

figure(1); clf;
surf(X,Y,reshape(A,J,J),'EdgeColor','none');
view(2)
drawnow;
tk = 0;
s=sum(Sp);
s

pause

% BOUCLE PRINCIPALE
while t < tfinal
    pause;
    newSp = C.*Sp + Sp + k*(-gama_.*Sp-beta_.*(1-ksi).*Sp.*A-beta_.*ksi_.*c1.*Sp.*Sp+delta.*R+mu.*(R)+mu_etoile.*A) + ...
        Dv*k/h^2*L*Sp;
    newA = A + k.*(gama_.*Sp+sigma.*R+beta_.*(1-ksi).*Sp.*A+beta_.*ksi_.*c1.*Sp.*Sp+nu.*R.*A-(zeta+mu_etoile).*A) +  ...
        Dw*k/h^2*L*A;
    newR = R + k*(zeta.*A-nu.*R.*A-(delta+sigma+mu).*R) +...
        Dr*k/h^2*L*R;
    Sp = newSp;
    A = newA;
    R = newR;
    if tk > 10 
        surf(X,Y,reshape(A,J,J),'EdgeColor','none');
        view(2  )
        drawnow;
        sp=sum(Sp);
        sp
        a=sum(Sp+R+A);
        a
        t
        
        
    end
    t = t + k;
    tk = tk + k;
end
