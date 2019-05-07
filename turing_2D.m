%% Turing 2D: FitzHugh-Nagumo 

% parametre des equations
epsilon = 0.08;
b = 0.7;
c = 0.8;
Dv = 0.01;
Dw = 2.5; 
I = 0.32;

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
    mu=.007288;
    mu_etoile=0.01155;
    R=0;
    
    c1=alpha/(epsilon+gama+mu);
    beta_=beta/(1+c1);
    gama_=gama/(1+c1);
    ksi_=ksi/(1+c1);
    c2=c1/(1+c1);
    
s = 20; % taille domaine (carre sxs)
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
v = zeros(J2,1); % stocke seulement l'etat au temps t
w = zeros(J2,1);
newv = zeros(J2,1);
neww = zeros(J2,1);


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

% interieur
L = L + sparse(interieur,interieur+1,1,J2,J2);
L = L + sparse(interieur,interieur-1,1,J2,J2);
L = L + sparse(interieur,interieur+J,1,J2,J2);
L = L + sparse(interieur,interieur-J,1,J2,J2);

% bords
L = L + sparse(bordhaut,bordhaut+1,1,J2,J2);
L = L + sparse(bordhaut,bordhaut+J-1,1,J2,J2);
L = L + sparse(bordhaut,bordhaut+J,1,J2,J2);
L = L + sparse(bordhaut,bordhaut-J,1,J2,J2);


L = L + sparse(bordgauche,bordgauche+1,1,J2,J2);
L = L + sparse(bordgauche,bordgauche-1,1,J2,J2);
L = L + sparse(bordgauche,bordgauche+J,1,J2,J2);
L = L + sparse(bordgauche,bordgauche+J*(J-1),1,J2,J2);

L = L + sparse(bordbas,bordbas-(J-1),1,J2,J2);
L = L + sparse(bordbas,bordbas-1,1,J2,J2);
L = L + sparse(bordbas,bordbas+J,1,J2,J2);
L = L + sparse(bordbas,bordbas-J,1,J2,J2);

L = L + sparse(borddroit,borddroit+1,1,J2,J2);
L = L + sparse(borddroit,borddroit-1,1,J2,J2);
L = L + sparse(borddroit,borddroit-J*(J-1),1,J2,J2);
L = L + sparse(borddroit,borddroit-J,1,J2,J2);

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
v = -1 + 0.1*(-1 + 2*rand(J^2,1)); 
w = -0.3 + 0.1*(-1 + 2*rand(J^2,1));


% parametres de simulation, temps
t0 = 0;
tfinal = 200; 
t = t0;
% k doit etre < a h^2/2/d/D ou dim d = 2
k = min(1,0.2*( h^2/4/max(Dv,Dw) ));
 

figure(1); clf;
surf(X,Y,reshape(v,J,J),'EdgeColor','none');
view(2)
drawnow;
tk = 0;
pause

% BOUCLE PRINCIPALE
while t < tfinal
    drawnow;
    newv = v + k*(-gama_.*v-beta_.*(1-ksi).*v.*w-beta_.*ksi_.*c1.*v.*v+delta.*R+mu.*(v.*c2+R)+mu_etoile.*w) + ...
        Dv*k/h^2*L*v;
    neww = w + k.*(gama_.*v+sigma.*R+beta_.*(1-ksi).*v.*w+beta_.*ksi_.*c1.*v.*v-nu.*R.*w-(zeta+mu_etoile).*w) +  ...
        Dw*k/h^2*L*w;
    v = newv;
    w = neww;
    if tk > 1 
        t
        surf(X,Y,reshape(v,J,J),'EdgeColor','none');
        view(2)
        drawnow;
        tk = 0;
    end
    t = t + k;
    tk = tk + k;
end
