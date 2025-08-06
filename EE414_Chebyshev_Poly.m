%--------------------------------------------------------------------

function [T_NP, uz, uz_poly] = ...
    EE414_Chebyshev_Poly(NP)

%--------------------------------------------------------------------
% [T_NP] = EE4XY_Chebyshev_Poly(NP);
% T_NP = Chebyshev Polynomial of Nth order
% Roots of T_NP(u)
%  T_NP(uz) = 0
%--------------------------------------------------------------------
% 12/06/2013 By John J. Burke, PhD, PE.
% 11/13/2014 By John J. Burke, PhD, PE.
% 11/05/2015 By John J. Burke, PhD, PE.
% 05/21/2020 By John J. Burke, PhD, PE.
% 06/03/2020 By John J. Burke, PhD, PE.
% 03/13/2021 By John J. Burke, PhD, PE.
% 02/06/2024 By John J. Burke, PhD, PE.
%--------------------------------------------------------------------

u = [1, 0];
T_NP_m1 = 1;
T_NP = [1, 0];

%--------------------------------------------------------------------

for ii = 2:NP
    T_Old = T_NP;
    TA =  2*conv(T_NP, u);
    TB =  [0, 0, T_NP_m1];
    T_NP = TA - TB;
    T_NP_m1 = T_Old;
end

%--------------------------------------------------------------------

uz_poly = roots(T_NP);
uz_poly = sort(uz_poly);

%--------------------------------------------------------------------

kk = 1 : 2 : (2*NP-1);
uz = cosd((kk/NP)*90);
uz = sort(uz);
uz = transpose(uz);

%--------------------------------------------------------------------
