% 0. Initialize Parameters

n = 1251;                   % Number of locations to evaluate bridge failure
L = 1250;                   % Length of bridge
x = linspace(0, L, n);      % Define x coordinate
total_loads = zeros(1, n);  % Initialize total_loads for SFD(x)
y_zero = zeros(1,n);        % Initialize y-axis

%% 1. Point Loading Analysis (SFD, BMD)
P = -185;
[total_loads] = ApplyPL(550, P, x, total_loads); % Construct load vector
[total_loads] = ApplyPL(L, P, x, total_loads); % Construct load vector
[SFD_PL, BMD_PL] = MakeSFDBMD(total_loads); % Construct SFD, BMD
plot(x, SFD_PL)
plot(x, BMD_PL)
set(gca,'YDir','reverse')

%% 2. Define cross-sections
% There are many (more elegant ways) to construct cross-section objects
xc = [0 550 L]; % Location, x, of cross-section change
bft = [100 100 100]; % Top Flange Width
tft = [2.54 2.54 2.54]; % Top Flange Thickness
hw = [100 120 100]; % Web Height
tw = [1.27 1.27 1.27]; % Web Thickness (Assuming 2 separate webs)
bfb = [80 80 80]; % Bottom Flange Width
tfb = [1.27 1.27 1.27]; % Bottom Flange Thickness
a = [400 400 400]; % Diaphragm Spacing
% Optional but you need to ensure that your geometric inputs are correctly implemented
% VisualizeBridge( {CrossSectionInputs} ); 

%% 3. Define Material Properties
SigT = 30;
SigC = 6;
E = 4000;
TauU = 4;
TauG = 2;
mu = 0.2;

function [ total_loads ] = ApplyPL( xP, P, x, total_loads )
% Constructs load vector from application of one point load. Assumes fixed location of supports
% Input: location and magnitude of point load. The previous load vector can be entered as input to 
% construct the load vector of multiple point loads.
% Output: total_loads, 1-D array of length n
% xP = distance from left (A) support (x = 0)
xP_B = 1060;
cur_P_B = -P * xP / xP_B; % P_B due to current P
P_A = total_loads(1) -P - cur_P_B; % Overall P_A
P_B = cur_P_B + total_loads(xP_B+1) - total_loads(xP_B); % Overall P_B

total_loads(1) = P_A;
total_loads(xP) = P;
total_loads(xP_B+1) = P_B;
end

function [ SFD_PL, BMD_PL ] = MakeSFDBMD(total_loads)
% Constructs SFD, BMD from all point loads applied.
% Input: total_loads vector containing the current loading of the beam.
% Output: SFD_PL, BMD_PL, both 1-D arrays of length n.
    n = 1251;
    SFD_PL = zeros(1, n);
    SFD_PL(1) = total_loads(1);
    for i = 2:length(total_loads)
        SFD_PL(i) = total_loads(i) + SFD_PL(i-1);
    end
    BMD_PL = cumsum(SFD_PL);
end