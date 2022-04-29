%% parametri caso autovalori immaginari puri

clear
clc

%parametri
m1 = 2;                                                                                                   
m2 = 4;
k = 2;
c = 0;

%parametri misurati
sigma = 0.1;
m1_m = m1 + sigma*randn;
m2_m = m2 + sigma*randn;
k_m = k + sigma*randn;
c_m = c + sigma*randn;

%% parametri caso autovalori complessi coniugati

clear
clc

%parametri
m1 = 2;                                                                                                   
m2 = 4;
k = 3;
c = 0.5;

%parametri misurati
sigma = 0.1;
m1_m = m1 + sigma*randn;
m2_m = m2 + sigma*randn;                                                                     
k_m = k + sigma*randn;
c_m = c + sigma*randn;

%%

%riferimento r(t) = M1*t + M0
M1 = 1;
M0 = 1;

%parametri del riferimento sinusoidale
w = 1;
M = 1;
phi = 2;

%matrici del sistema P
A_n = [0 0 1 0; 0 0 0 1; -k/m1 k/m1 -c/m1 c/m1; k/m2 -k/m2 c/m2 -c/m2];
B_n = [0; 0; 1/m1; 0];
C_n = [1 0 0 0];
D_n = 0;
%x0 = [1;10;3;2];
x0 = [0;0;0;0];

%matrici del sistema P con errori di misura nei parametri
A = [0 0 1 0; 0 0 0 1; -k_m/m1_m k_m/m1_m -c_m/m1_m c_m/m1_m; k_m/m2_m -k_m/m2_m c_m/m2_m -c_m/m2_m];
B = [0; 0; 1/m1_m; 0];
C = [1 0 0 0];
D = 0;


%matrici del sistema C_M
%C_M = 1/(s^2 + w^2);
Am = [0 -w^2; 1 0];
Bm = [1; 0];
Cm = [0 1];
Dm = 0;

%matrici della cascata C_M, P
A_segn = [Am zeros(2, 4); B*Cm A];
B_segn = [Bm; B*Dm];
C_segn = [D*Cm C];
D_segn = [D*Dm];

%Ackermann per il calcolo di F
lambda1 = -2 - 1i;
lambda2 = -1;
lambda3 = -1.5;
lambda4 = -2 + 1i;
lambda5 = -1 + 1i;
lambda6 = -1 - 1i;



P_segn = [B_segn A_segn*B_segn (A_segn^2)*B_segn (A_segn^3)*B_segn (A_segn^4)*B_segn (A_segn^5)*B_segn];
P_segn_inv = inv(P_segn);
tau_n = P_segn_inv(6, :);
poles = (A_segn - lambda1*eye(6))*(A_segn - lambda2*eye(6))*(A_segn - lambda3*eye(6))*(A_segn - lambda4*eye(6))*(A_segn - lambda5*eye(6))*(A_segn - lambda6*eye(6));
F = -tau_n*poles;
F = real(F);
%eig(A_segn + B_segn*F)


%Ackermann per il calcolo di V
A_segn_d = A_segn';
B_segn_d = C_segn';
P_segn_d = [B_segn_d A_segn_d*B_segn_d (A_segn_d^2)*B_segn_d (A_segn_d^3)*B_segn_d (A_segn_d^4)*B_segn_d (A_segn_d^5)*B_segn_d];
P_segn_d_inv = inv(P_segn_d);
tau_n = P_segn_d_inv(6, :);

lambda1 = -2;
lambda2 = -4 + 2i;
lambda3 = -4 - 2i;
lambda4 = -3;
lambda5 = -2 +2i;
lambda6 = -2 -2i;

poles = (A_segn_d - lambda1*eye(6))*(A_segn_d - lambda2*eye(6))*(A_segn_d - lambda3*eye(6))*(A_segn_d - lambda4*eye(6))*(A_segn_d - lambda5*eye(6))*(A_segn_d - lambda6*eye(6));
F_d = -tau_n*poles;
V = -F_d';
%eig(A_segn - V*C_segn)

%matrici dello stabilizzatore
As = A_segn - V*C_segn + B_segn*F - V*D_segn*F;
Bs = V;
Cs = -F;
Ds = 0;

%matrici del compensatore
Ac = [As zeros(6, 2); Bm*Cs Am];
Bc = [Bs; Bm*Ds];
Cc = [Dm*Cs Cm];
Dc = Dm*Ds;


