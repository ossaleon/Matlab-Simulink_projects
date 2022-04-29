%% parametri caso autovalori immaginari puri

clear
clc


%parametri
m1 = 2;                                                                                                   
m2 = 3;
k = 3;
c = 0;

%parametri misurati
sigma = 0;
m1_m = m1 + sigma*randn;
m2_m = m2 + sigma*randn;
k_m = k + sigma*randn;
c_m = c + sigma*randn;

%% parametri caso autovalori complessi coniugati

clear
clc

%parametri
m1 = 2;                                                                                                   
m2 = 3;
k = 3;
c = 0.5;

%parametri misurati
sigma = 0;
m1_m = m1 + sigma*randn;
m2_m = m2 + sigma*randn;                                                                     
k_m = k + sigma*randn;
c_m = c + sigma*randn;

%%

%riferimento r(t) = M1*t + M0
M1 = 1;
M0 = 2;

%parametri del riferimento sinusoidale
w = 1;
M = 2;
phi = pi/2;

%matrici del sistema P nominali
A_n = [0 0 1 0; 0 0 0 1; -k/m1 k/m1 -c/m1 c/m1; k/m2 -k/m2 c/m2 -c/m2];
B_n = [0; 0; 1/m1; 0];
C_n = [1 0 0 0];
D_n = 0;
%x0 = [1;10;3;2];
x0 = [0; 2; 0; 0];

%matrici del sistema P con errori di misura nei parametri
A = [0 0 1 0; 0 0 0 1; -k_m/m1_m k_m/m1_m -c_m/m1_m c_m/m1_m; k_m/m2_m -k_m/m2_m c_m/m2_m -c_m/m2_m];
B = [0; 0; 1/m1_m; 0];
C = [1 0 0 0];
D = 0;


%il sistema (A, B, C, D) contiene già due poli in 0
%perciò Cm = 1;
A_segn = A;
B_segn = B;
C_segn = C;
D_segn = D;


%Ackermann per il calcolo di F
%lambda1 = -2 - 1i;
%lambda2 = -1;
%lambda3 = -1.5;
%lambda4 = -2 + 1i;

lambda1 = -0.75 - 1i;
lambda2 = -0.25;
lambda3 = -0.75;
lambda4 = -0.75 + 1i;




P_segn = [B_segn A_segn*B_segn (A_segn^2)*B_segn (A_segn^3)*B_segn];
P_segn_inv = inv(P_segn);
tau_n = P_segn_inv(4, :);
poles = (A_segn - lambda1*eye(4))*(A_segn - lambda2*eye(4))*(A_segn - lambda3*eye(4))*(A_segn - lambda4*eye(4));
F = -tau_n*poles;
F = real(F);
%eig(A_segn + B_segn*F)


%Ackermann per il calcolo di V
A_segn_d = A_segn';
B_segn_d = C_segn';
P_segn_d = [B_segn_d A_segn_d*B_segn_d (A_segn_d^2)*B_segn_d (A_segn_d^3)*B_segn_d];
P_segn_d_inv = inv(P_segn_d);
tau_n = P_segn_d_inv(4, :);

%lambda1 = -2;
%lambda2 = -2 + 1i;
%lambda3 = -2 - 1i;
%lambda4 = -2.5;

lambda1 = -1;
lambda2 = -1 + 1i;
lambda3 = -1 - 1i;
lambda4 = -1.5;

poles = (A_segn_d - lambda1*eye(4))*(A_segn_d - lambda2*eye(4))*(A_segn_d - lambda3*eye(4))*(A_segn_d - lambda4*eye(4));
F_d = -tau_n*poles;
V = -F_d';

%matrici del compensatore
Ac= A_segn - V*C_segn + B_segn*F - V*D_segn*F;
Bc = V;
Cc = -F;
Dc = 0;
