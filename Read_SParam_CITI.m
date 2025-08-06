%--------------------------------------------------------------------

function [freq, S, Mult, S_Type] = ...
    Read_SParam_CITI(FileName)

%--------------------------------------------------------------------
% Read_SParam_CITI
%--------------------------------------------------------------------
%   Written by John J. Burke 08/01/2017
%   Modified by John J. Burke 08/01/2017
%   Modified by John J. Burke 03/23/2018
%   Modified by John J. Burke 05/05/2018
%   Modified by John J. Burke 03/30/2020
%   Modified by John J. Burke 05/10/2022
%   Modified by John J. Burke 05/10/2024
%--------------------------------------------------------------------

Mult = 1;
fid = fopen(FileName);
NPorts = 0;

%--------------------------------------------------------------------

while ~feof(fid) % loop over the following until the end of the file is reached.
    line = fgets(fid); % read in one line
    if contains(line, 'VAR freq')
        [S_Type] = Sparam_Type(line);
        Data = ...
            str2num(regexprep(line, '^(\D+)(\d+)(.*)', '$2'));
        N_freq = Data;
    end
    if contains(line, 'DATA PortName')
        NPorts = NPorts + 1;
    end
    if contains(line, 'BEGIN')
        break
    end
end
    
%--------------------------------------------------------------------

[freq, Nfreq] = ...
    Freq_Extract(fid, N_freq);
S = zeros(Nfreq, NPorts, NPorts);

%--------------------------------------------------------------------

for jj = 1 : NPorts
    for kk = 1 : NPorts
        S(:, jj, kk) = ...
            SParam_Read(fid, N_freq, S_Type);
    end
end

%--------------------------------------------------------------------

fclose(fid);

%--------------------------------------------------------------------











%--------------------------------------------------------------------

function S_ij = ...
    SParam_Read(fid, Nfreq, S_Type)

%--------------------------------------------------------------------

S_ij = zeros(Nfreq, 1);
fgets(fid); % read in one line

%--------------------------------------------------------------------

for kk = 1 : Nfreq
        line = fgets(fid); % read in one line
        Data = str2num(line);
        C1 = Data(1);
        C2 = Data(2);
        S_ij(kk, 1) = ...
            Sparam_Extract(C1, C2, S_Type);
end

%--------------------------------------------------------------------

fgets(fid); % read in one line

%--------------------------------------------------------------------










%--------------------------------------------------------------------

function [S_Type] = ...
    Sparam_Type(line)

%--------------------------------------------------------------------

line = lower(line);
if contains(line, 'ma')
    S_Type = 'MA';
elseif contains(line, 'db')
    S_Type = 'dB';
elseif contains(line, 'ri')
    S_Type = 'RI';
else
    S_Type = 'MA';
end

%--------------------------------------------------------------------









%--------------------------------------------------------------------

function [S] = ...
    Sparam_Extract(C1, C2, S_Type)

%--------------------------------------------------------------------

switch S_Type
    case 'MA'
        S = Polar_2_Rect(C1, C2);
    case 'RI'
        S = C1 + 1j*C2;
    case 'dB'
        C1 = 10.^(C1/20);
        S = Polar_2_Rect(C1, C2);
end

%--------------------------------------------------------------------









%--------------------------------------------------------------------

function [freq, N_freq] = ...
    Freq_Extract(fid, N_freq)

%--------------------------------------------------------------------

freq = zeros(N_freq, 1);
I_Flag = 0;

%--------------------------------------------------------------------

for kk = 1 : N_freq
    line = fgets(fid); % read in one line
    Data = str2num(line);
    freq(kk) = Data;
    if contains(line, 'END')
        I_Flag = 1;
        break
    end
end
    
%--------------------------------------------------------------------

if (I_Flag == 0)
        fgets(fid); % read in one line
else
    if (kk < N_freq), N_freq = kk; end
end
%--------------------------------------------------------------------
