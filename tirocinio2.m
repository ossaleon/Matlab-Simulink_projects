clear
clc

%input
gamma_1 = -2;
gamma_2 = 1;
A = [0 1; gamma_1 gamma_2];
B = [0; 1];
%realizzo un sistema lineare date le matrici A, B, C, D
%C e D sono ininfluenti per la risposta nello stato, quindi le poniamo
%nulle
%l'ingresso è arbitrario
C = [0, 0; 0 0];
D = [0; 0];

sys = ss(A,B,C,D);

t_segn = 3;
T = [1:0.1:t_segn+3];

%% ingresso rampa a partire da 0
t = (0:0.1:t_segn)';
unitstep = t>=0;
u0 = t.*unitstep;

%% rampa quadratica a partire da 0
t = (0:0.1:t_segn)';
unitstep = t>=0;
ramp = t.*unitstep;
u0 = t.^2.*unitstep;

%% ingresso costante
t = (0:0.1:t_segn)';
k = 1;
u0 = k.*ones(t_segn/0.1+1, 1);

%%

%calcolo la risposta nello stato (che, partendo dallo stato nullo,
%coincide con la risposta forzata nello stato)
%lsim calcola la storia della risposta nello stato
[~,~,x] = lsim(sys, u0, t);

x_segn = (x(t_segn/0.1 + 1, :))';

%vettore J_Ti(u_star)
J = [];
for Ti = T
   u_star = controllo_ottimo(A, B, x_segn, Ti);
   temp = indice_costo_minimo(u_star);
   J = [J temp];
end 
J_t_segn = trapz(t, u0.^2);

%% grafico di J_Ti(u*) al variare di Ti
figure(1)
plot(T, J, 'b')
hold on
plot(T, J_t_segn*ones(size(T)), 'r')
title('Andamento di J al variare del tempo finale T_i')
xlabel('T_i')
ylabel('J')
legend('J_T(u*)', 'J_tsegn(u_0)')
hold off

%% grafico di u* nell'intervallo [0, Ti] e delle componenti della corrispondente risposta nello stato x*
T_fin = 5;
u_star = controllo_ottimo(A, B, x_segn, T_fin);
figure(2)
plot([0:0.1:T_fin], u_star', 'g');
xlabel('t')
ylabel('u*_Tfin')
title("Andamento di u* nell'intervallo [0, T_i]")

[~,~,x_star] = lsim(sys, u_star, [0:0.1:T_fin]');
figure(3)
plot([0:0.1:T_fin], (x_star(:, 1))', 'r')
hold on
plot([0:0.1:T_fin], (x_star(:, 2))', 'b')
xlabel('t')
ylabel('x*_1 - x*_2')
legend('x*_1', 'x*_2')
title("Andamento di x*_1 e x*_2 nell'intervallo [0, T_i]")
hold off

%% FUNZIONI

%calcolo della matrice gramiana G al tempo finale T (funzione di T)
function G = calcolo_G(A, B)
    syms tau T_fin
    exp_matrix = expm(A*tau);
    product = exp_matrix*B*B'*exp_matrix';
    G = int(product, tau, 0, T_fin);
end   

%restituisce lo storico della u_star (in colonna) dato uno stato da
%raggiungere x e un tempo finale T
function u_star = controllo_ottimo(A, B, x, T)
    temp = [0: 0.1: T];
    syms G_integrale(T_fin)
    G_integrale(T_fin) = calcolo_G(A, B);
    G = G_integrale(T);
    beta = G\x;
    u_star = [];
    [~, cols] = size(temp);
    for i = 1:cols
        u = B'*expm(A'*(T-temp(1, i)))*beta;
        u= double(u);
        u_star = [u_star; u];
    end    
end
 
%restituisce l'indice di costo minimo di un dato controllo ottimo u_star
function J = indice_costo_minimo(u_star)
   [rows, ~] = size(u_star);
   t = [0:0.1:(rows-1)*0.1]';
   J = trapz(t, u_star.^2);
end
