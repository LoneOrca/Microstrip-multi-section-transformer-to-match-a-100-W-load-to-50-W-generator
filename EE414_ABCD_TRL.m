%------------------------------
function T = EE414_ABCD_TRL(Z0x, theta_x)
%------------------------------
% ABCD matrix Transmission line
%------------------------------
j=1j;
K=cosd(theta_x);
L=sind(theta_x);
Y0x=1/Z0x;
A=K;
B=j*Z0x*L;
C=j*Y0x*L;
D=K;
T = [A, B; C, D];
%------------------------------