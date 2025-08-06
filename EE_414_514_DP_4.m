%----------------------------------------------------------------------
%% EE-414-514-sn620094-Design Project 4
%% Chebyshev
%----------------------------------------------------------------------

clearvars
clc
close all

%----------------------------------------------------------------------
% Units 
%----------------------------------------------------------------------

G = 10^9;
Meg = 10^6;
k = 10^+3;
c = 10^-2;
m = 10^-3;
u = 10^-6;
n = 10^-9;
%----------------------------------------------------------------------
% 
%----------------------------------------------------------------------
IFigure = 0;
NF = 32;
dfreq = 1;
df = 1*Meg;
j = 1*j;

%----------------------------------------------------------------------
Zg = 50;
ZL = 100;
Z0 = Zg;
Znp = ZL;
Gamma_L = (1/2)*log(ZL/Z0);
Pbmax_dB = -30; % dB given
pbmax = 10^(Pbmax_dB/20);
%----------------------------------------------------------------------
f0 = 8*G;
BWf = 6.4*G;
fL = f0 - (1/2)*BWf;
fH = f0 + (1/2)*BWf;

%----------------------------------------------------------------------
theta_0 = 90; % Degrees
theta_L = (fL/f0)*theta_0;
theta_H = (fH/f0)*theta_0;
uL = + 1;
u0 = uL/cosd(theta_L);

N_P = acosh(abs(Gamma_L)/pbmax ) / acosh(u0);
NP= ceil(N_P);
NT = NP;
%----------------------------------------------------------------------
uo =  cosh((1/NP) * acosh(Gamma_L / pbmax));
new_theta_L = acosd(uL/uo);
new_theta_H = 2*theta_0 - new_theta_L;
BW_theta = new_theta_H - new_theta_L;
BW_theta_percent = (BW_theta/theta_0)*100;
%----------------------------------------------------------------------
new_fL = (new_theta_L/theta_0)*f0;
new_fH = (new_theta_H/theta_0)*f0;
Print_Real_Unit('new_fL',new_fL,'Hz')
Print_Real_Unit('new_fH',new_fH,'Hz')
%----------------------------------------------------------------------

% u = uo*cosd(theta_0);
% uz = 0;
% j = 1*j;
% A = sign(Gamma_L) * pbmax;
% 
% w = exp(j * (2 * theta_0));
% theta_zi = acos(uz/uo);
% wzi = exp(j * (2*theta_zi)); % stage 2
% Print_Polar('wzi',wzi)

%[T_N, uzi, uz_poly, X_Sign] = EE322_Chebyshev_Poly(N_P);
[T_NP, uz, uz_poly] = EE414_Chebyshev_Poly(NP);

Print_Real('uz',uz)
Print_Real('T_N',T_NP)
Print_Real('Uz_poly',uz_poly )
%Print_Real('X_sign',X_Sign)

theta_zi = acosd(uz/uo);
wzi = Polar_2_Rect(1,-2*theta_zi);
In = real(poly(wzi));
Aw = Gamma_L / sum(In);
Gamma_k = Aw * In;
% Z1 = Z0^(3/4) * ZL^(1/4);
% Z2 = Z0^(1/4) * ZL^(3/4);
Z0 = 50;

Z1 = Z0*((1+ Gamma_k(1))/(1 - Gamma_k(1)));
Z2 = Z1*((1+ Gamma_k(2))/(1 - Gamma_k(2)));
Z3 = Z2*((1+ Gamma_k(3))/(1 - Gamma_k(3)));
Z4 = Z3*((1+ Gamma_k(4))/(1 - Gamma_k(4)));
Print_Break

Print_Polar('wzi',wzi)
Print_Real('Ik',In)

Print_Real('Gamma_k',Gamma_k)
Print_Real_Unit('Z1',Z1,'Ohms')
Print_Real_Unit('Z2',Z2,'Ohms')
Print_Real_Unit('Z3',Z3,'Ohms')
Print_Real_Unit('Z4',Z4,'Ohms')

RLmin = 30;

S11_dB_RL = -RLmin;

S11_max_RL = 10^(((S11_dB_RL)/20));

S21min_RL = sqrt(1-abs(S11_max_RL)^2);
Print_Real('S21_min (RL)',S21min_RL,'W/W'); %Good
S21_dB_RL = 20*log10(S21min_RL); % IL
Print_Real('S21_min (RL)',S21_dB_RL,'dB'); %Good

Ap_dB = abs(S21_dB_RL);
%Ap_dB = 0.02; % Round down
Print_Real('Ap (RIPPLE)',Ap_dB,'dB'); %Good

delta = ((fH - fL)/f0);

wp = 2*pi*f0;
Ap = 10^(Ap_dB/10); % W/W
Xp = Ap - 1;
epsilon = sqrt(Xp);
% ALS = 26.0206;
% ALS_W = 10^(ALS/10);
% XLS = ALS_W - 1;
% Omega_LS = (1/delta)*((fls/f0)-(f0/fls));
% NLS = acosh(sqrt(XLS)/epsilon ) / acosh(abs(Omega_LS));
% N_LS = ceil(NLS);
% NLS_Elements = ceil(N_LS/2);
% Print_Real('NLS',NLS)
% %----------------------------------------------------------------------
% % FHS
% %----------------------------------------------------------------------
% AHS = 26.0206;
% AHS_W = 10^(AHS/10);
% XHS = AHS_W - 1;
% Omega_HS = (1/delta)*((fhs/f0)-(f0/fhs));
% NHS = acosh(sqrt(XHS)/epsilon ) / acosh(abs(Omega_HS));
% Print_Real('NHS',NHS)



f_min = 0*G;
f_max = 16*G;
freq = f_min : df : f_max;
freq = sort(freq);
freq = freq';
% I_f0_Amp = freq_Amp == f0;
%Z0 = 50;
N_Freq = length(freq);
S_Filter = zeros(N_Freq,2,2);
for kk = 1 : N_Freq
fk = freq(kk);
theta_P = (fk / f0 )* theta_0;
%T0 = eye(2); % ?
T1 = EE414_ABCD_TRL(Z0, theta_P);
T2 = EE414_ABCD_TRL(Z1, theta_P);
T3 = EE414_ABCD_TRL(Z2,theta_P);
T4 = EE414_ABCD_TRL(Z3,theta_P);
T5 = EE414_ABCD_TRL(ZL, theta_P);
T = T1*T2*T3*T4*T5;
S_Filter(kk, :, :) = ABCD_to_S(T, [Z0, ZL]);
end

%----------------------------------------------------------------------
%
%----------------------------------------------------------------------

S11_Filter = S_Filter(:,1,1);
S11_Filter_Mag = abs(S11_Filter);
S11_Filter_dB = 20*log10(S11_Filter_Mag);
S21_Filter = S_Filter(:,2,1);
S21_Filter_Mag = abs(S21_Filter);
S21_Filter_dB = 20*log10(S21_Filter_Mag);


IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S11_Filter_dB, 'r', 'linewidth', 6)
hold on
plot(new_fL/G,interp1(freq/G,S11_Filter_dB,new_fL/G),'ro','linewidth',10);
plot(new_fH/G,interp1(freq/G,S11_Filter_dB,new_fH/G),'ro','linewidth',10);
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{11} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
%title('MATLAB Ideal Lumped Elements')
legend('|{\itS}_{11}|')
axis([f_min/G, f_max/G, -40, 0])
set(gca, 'xtick', f_min/G : 1 : f_max/G);
set(gca, 'linewidth', 2.0)

% 
IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S21_Filter_dB, 'r', 'linewidth', 6)
hold on
plot(new_fL/G,interp1(freq/G,S21_Filter_dB,new_fL/G),'ro','linewidth',10);
plot(f0/G,interp1(freq/G,S21_Filter_dB,f0/G),'ro','linewidth',10);
plot(new_fH/G,interp1(freq/G,S21_Filter_dB,new_fH/G),'ro','linewidth',10);
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{21} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
%title('MATLAB Ideal Lumped Elements')
legend('|{\itS}_{21}|')
axis([f_min/G, f_max/G, -0.25, 0])
set(gca, 'xtick', f_min/G : 1 : f_max/G);
set(gca, 'linewidth', 2.0)
%-----------------------------------------------


%% Butterworth
N_P1 = log(Gamma_L/pbmax)/log(u0);
NP1= ceil(N_P1);
NT1 = NP1;
%----------------------------------------
uo1 =  exp((1/NP1)*log(Gamma_L/pbmax));
new_theta_L1 = acosd(uL/uo1);
new_theta_H1 = 2*theta_0 - new_theta_L1;
BW_theta1 = new_theta_H1 - new_theta_L1;
BW_theta_percent1 = (BW_theta1/theta_0)*100;
%----------------------------------------------------------------------
new_fL1 = (new_theta_L1/theta_0)*f0;
new_fH1 = (new_theta_H1/theta_0)*f0;
Print_Real_Unit('new_fL1',new_fL1,'Hz')
Print_Real_Unit('new_fH1',new_fH1,'Hz')
%---------------------------------------------------------------------
u1 = uo1*cosd(theta_0);
 uz1 = 0;
 j = 1*j;
 A = sign(Gamma_L) * pbmax;
 
 w1 = exp(j * (2 * theta_0));
 theta_zi1 = acos(uz/uo);
 wzi1 = exp(j * (2*theta_zi1)); % stage 2
 Print_Polar('wzi1',wzi1)
 Print_Polar('w1',w1)
%IK = real(poly(wzi));
% % Butterworth
 Aw1 = 2^(-NP1)*Gamma_L;
 C_k_N = zeros(NP1,1);
 for kk = 0 : NP1
     C_k_N(kk+1) = nchoosek(NP1,kk); % 
 end
 Ik = C_k_N;
 Gamma_k1 = Aw1 * Ik;
 %---------------------------------------------------------------------
Z11 = Z0*exp(2*Gamma_k1(1));
Z21 = Z11*exp(2*Gamma_k1(2));
Z31 = Z21*exp(2*Gamma_k1(3));
Z41 = Z31*exp(2*Gamma_k1(4));
Z51 = Z41*exp(2*Gamma_k1(5));
Z61 = Z51*exp(2*Gamma_k1(6));
Print_Break
Print_Real('Gamma_k1',Gamma_k1)
Print_Real_Unit('Z11',Z11,'Ohms')
Print_Real_Unit('Z21',Z21,'Ohms')
Print_Real_Unit('Z31',Z31,'Ohms')
Print_Real_Unit('Z41',Z41,'Ohms')
Print_Real_Unit('Z51',Z51,'Ohms')
Print_Real_Unit('Z61',Z61,'Ohms')

RLmin1 = 30;

S11_dB_RL1 = -RLmin1;

S11_max_RL1 = 10^(((S11_dB_RL1)/20));

S21min_RL1 = sqrt(1-abs(S11_max_RL1)^2);
%Print_Real('S21_min (RL)',S21min_RL1,'W/W'); %Good
S21_dB_RL1 = 20*log10(S21min_RL1); % IL
%Print_Real('S21_min (RL)',S21_dB_RL,'dB'); %Good

Ap_dB1 = abs(S21_dB_RL1);
%Ap_dB = 0.02; % Round down
%Print_Real('Ap (RIPPLE)',Ap_dB,'dB'); %Good

delta1 = ((fH - fL)/f0);

wp1 = 2*pi*f0;
Ap1 = 10^(Ap_dB1/10); % W/W
Xp1 = Ap1 - 1;
epsilon1 = sqrt(Xp1);
f_min1 = 0*G;
f_max1 = 16*G;
freq1 = f_min1 : df : f_max1;
freq1 = sort(freq1);
freq = freq1';
% I_f0_Amp = freq_Amp == f0;
%Z0 = 50;
N_Freq = length(freq);
S_Filter = zeros(N_Freq,2,2);
for kk = 1 : N_Freq
fk = freq(kk);
theta_P = (fk / f0 )* theta_0;
%T0 = eye(2); % ?
T11 = EE414_ABCD_TRL(Z0, theta_P);
T21 = EE414_ABCD_TRL(Z11, theta_P);
T31 = EE414_ABCD_TRL(Z21, theta_P);
T41 = EE414_ABCD_TRL(Z31, theta_P);
T51 = EE414_ABCD_TRL(Z41, theta_P);
T61 = EE414_ABCD_TRL(Z51, theta_P);
%T7 = EE414_ABCD_TRL(Z6, theta_P);
T81 = EE414_ABCD_TRL(ZL, theta_P);
T1 = T11*T21*T31*T41*T51*T61*T81;
S_Filter(kk, :, :) = ABCD_to_S(T1, [Z0, ZL]);
end

S11_Filter1 = S_Filter(:,1,1);
S11_Filter_Mag1 = abs(S11_Filter1);
S11_Filter_dB1 = 20*log10(S11_Filter_Mag1);
S21_Filter1 = S_Filter(:,2,1);
S21_Filter_Mag1 = abs(S21_Filter1);
S21_Filter_dB1 = 20*log10(S21_Filter_Mag1);


IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S11_Filter_dB1, 'r', 'linewidth', 6)
hold on
plot(new_fL/G,interp1(freq/G,S11_Filter_dB1,new_fL/G),'ro','linewidth',10);
plot(new_fH/G,interp1(freq/G,S11_Filter_dB1,new_fH/G),'ro','linewidth',10);
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{11} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
%title('MATLAB Ideal Lumped Elements')
legend('|{\itS}_{11}|')
axis([f_min1/G, f_max1/G, -40, 0])
set(gca, 'xtick', f_min1/G : 2 : f_max1/G);
set(gca, 'linewidth', 2.0)

% 
IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S21_Filter_dB1, 'r', 'linewidth', 6)
hold on
plot(new_fL/G,interp1(freq/G,S21_Filter_dB1,new_fL/G),'ro','linewidth',10);
plot(f0/G,interp1(freq/G,S21_Filter_dB1,f0/G),'ro','linewidth',10);
plot(new_fH/G,interp1(freq/G,S21_Filter_dB1,new_fH/G),'ro','linewidth',10);
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{21} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
%title('MATLAB Ideal Lumped Elements')
legend('|{\itS}_{21}|')
axis([f_min1/G, f_max1/G, -0.25, 0])
set(gca, 'xtick', f_min1/G : 2 : f_max1/G);
set(gca, 'linewidth', 2.0)
%-----------------------------------------------

IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S21_Filter_dB1, 'r', 'linewidth', 6)
hold on
plot(freq/G, S21_Filter_dB, 'g', 'linewidth', 6)
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{11} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
%title('MATLAB Ideal Lumped Elements')
legend('Butterworth','Chebyshev','HFSS w/ Loss')
axis([f_min/G, f_max/G, -0.25, 0])
set(gca, 'xtick', f_min/G : 1 : f_max/G);
set(gca, 'linewidth', 2.0)




