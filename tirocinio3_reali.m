% Autovalori reali e distinti
clear
clc

%parametri 
m1 = 3;                                                                                                   
m2 = 4;
k = 3;
c = 5;

%parametri misurati
sigma = 0.2;
m1_m = m1  + sigma*randn;                                                                                                   
m2_m = m2  + sigma*randn;
k_m = k  + sigma*randn;
c_m = c + sigma*randn;

%matrici del sistema primale
A_n = [0 0 1 0; 0 0 0 1; -k/m1 k/m1 -c/m1 c/m1; k/m2 -k/m2 c/m2 -c/m2];
B_n = [0; 0; 1/m1; 0];
C_n = [1 0 0 0];
M_n = [0; 0; 0; 1/m2];
D_n = 0;

%matrici del sistema primale con errori di misura nei parametri
A = [0 0 1 0; 0 0 0 1; -k_m/m1_m k_m/m1_m -c_m/m1_m c_m/m1_m; k_m/m2_m -k_m/m2_m c_m/m2_m -c_m/m2_m];
B = [0; 0; 1/m1_m; 0];
C = [1 0 0 0];
M = [0; 0; 0; 1/m2_m];
D = 0;

%matrici del sistema duale
Ad = A';
Bd = C';

Pd = [Bd Ad*Bd Ad*Ad*Bd Ad*Ad*Ad*Bd];
Pd_inv = inv(Pd);
tau_n = Pd_inv(4, :);

%ampiezza del disturbo di x
d = 0;
%potenza del disturbo di y
np = 0;

%% Autovalori reali e distinti, convergenza lenta

lambda1 = -2 - 1i;
lambda2 = -1;
lambda3 = -1.5;
lambda4 = -2 + 1i;

poles = (Ad - lambda1*eye(4))*(Ad - lambda2*eye(4))*(Ad - lambda3*eye(4))*(Ad - lambda4*eye(4));
Fd = -tau_n*poles;
Fd = real(Fd);
%eig(Ad + Bd*Fd) %per la verifica degli autovalori del sistema a ciclo
%chiuso

V = -Fd';
H = A-V*C;
R = B-V*D;
%x0 = [1;10;3;2];
%x0 = [0;2;0;0];
x0 = [0;0;0;0];
%csi0 = [0.9;10;3.1;1.9];
%csi0 = [0;0;0;0]; %se non si conosce la stima di x
csi0 = [1;2;1;0];
e0 = x0 - csi0;

%% Autovalori reali e distinti, convergenza veloce

%lambda1 = -3;
%lambda2 = -5 + 2i;
%lambda3 = -5 - 2i;
%lambda4 = -4;

lambda1 = -2;
lambda2 = -4 + 2i;
lambda3 = -4 - 2i;
lambda4 = -3;


poles = (Ad - lambda1*eye(4))*(Ad - lambda2*eye(4))*(Ad - lambda3*eye(4))*(Ad - lambda4*eye(4));
Fd = -tau_n*poles;
Fd = real(Fd);
%eig(Ad + Bd*Fd) %per la verifica degli autovalori del sistema a ciclo
%chiuso

V = -Fd';
H = A-V*C;
R = B-V*D;
x0 = [1;10;3;2];
%x0 = [0;2;0;0];
%csi0 = [0.9;10;3.1;1.9];
csi0 = [0;0;0;0]; %se non si conosce la stima di x
e0 = x0 - csi0;