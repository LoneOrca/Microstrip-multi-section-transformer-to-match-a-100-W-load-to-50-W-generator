%----------------------------------------------------------------------
%% EE-414-514-sn620094-Design Project 4
%% Butterworth
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
%----------------------------------------------------------------------
%-----------------------------------------------
N_P = log(Gamma_L/pbmax)/log(u0);
NP= ceil(N_P);
NT = NP;
%----------------------------------------------------------------------
uo =  exp((1/NP)*log(Gamma_L/pbmax));
new_theta_L = acosd(uL/uo);
new_theta_H = 2*theta_0 - new_theta_L;
BW_theta = new_theta_H - new_theta_L;
BW_theta_percent = (BW_theta/theta_0)*100;
%----------------------------------------------------------------------
new_fL = (new_theta_L/theta_0)*f0;
new_fH = (new_theta_H/theta_0)*f0;
Print_Real_Unit('new_fL',new_fL,'Hz')
Print_Real_Unit('new_fH',new_fH,'Hz')
%---------------------------------------------------------------------
u = uo*cosd(theta_0);
 uz = 0;
 j = 1*j;
 A = sign(Gamma_L) * pbmax;
 
 w = exp(j * (2 * theta_0));
 theta_zi = acos(uz/uo);
 wzi = exp(j * (2*theta_zi)); % stage 2
 Print_Polar('wzi',wzi)
 Print_Polar('w',w)
%IK = real(poly(wzi));
% % Butterworth
 Aw = 2^(-NP)*Gamma_L;
 C_k_N = zeros(NP,1);
 for kk = 0 : NP
     C_k_N(kk+1) = nchoosek(NP,kk); % 
 end
 Ik = C_k_N;
 Gamma_k = Aw * Ik;
 %---------------------------------------------------------------------
 Z0 = 50;

Z1 = Z0*exp(2*Gamma_k(1));
Z2 = Z1*exp(2*Gamma_k(2));
Z3 = Z2*exp(2*Gamma_k(3));
Z4 = Z3*exp(2*Gamma_k(4));
Z5 = Z4*exp(2*Gamma_k(5));
Z6 = Z5*exp(2*Gamma_k(6));
Print_Break

%Print_Polar('wzi',wzi)
%Print_Real('Ik',Ik)

Print_Real('Gamma_k',Gamma_k)
Print_Real_Unit('Z1',Z1,'Ohms')
Print_Real_Unit('Z2',Z2,'Ohms')
Print_Real_Unit('Z3',Z3,'Ohms')
Print_Real_Unit('Z4',Z4,'Ohms')
Print_Real_Unit('Z5',Z5,'Ohms')
Print_Real_Unit('Z6',Z6,'Ohms')

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
T3 = EE414_ABCD_TRL(Z2, theta_P);
T4 = EE414_ABCD_TRL(Z3, theta_P);
T5 = EE414_ABCD_TRL(Z4, theta_P);
T6 = EE414_ABCD_TRL(Z5, theta_P);
%T7 = EE414_ABCD_TRL(Z6, theta_P);
T8 = EE414_ABCD_TRL(ZL, theta_P);
T = T1*T2*T3*T4*T5*T6*T8;
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
set(gca, 'xtick', f_min/G : 2 : f_max/G);
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
set(gca, 'xtick', f_min/G : 2 : f_max/G);
set(gca, 'linewidth', 2.0)
%-----------------------------------------------
I_f0 = find(freq==f0);
S11_f0 = S11_Filter_dB(I_f0);
Pint_Polar('S11',S11_f0)

