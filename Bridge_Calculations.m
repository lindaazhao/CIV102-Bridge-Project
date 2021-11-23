% 0. Initialize Parameters

n = 1251;               % Number of locations to evaluate bridge failure
L = 1250;               % Length of bridge
x = linspace(0, L, n);  % Define x coordinate
% SFD_PL = zeros(1, n)    % Initialize SFD(x)
total_loads = zeros(1, n); % Initialize total_loads for SFD(x)
y_zero = zeros(1,n);    % Initialize y-axis
BMD_PL = zeros(1,n);    % Initialize BMD(x)

%% 1. Point Loading Analysis (SFD, BMD)
P = -185;
[total_loads, BMD_PL] = ApplyPL(550, P, x, total_loads); % Construct SFD, BMD
[total_loads, BMD_PL] = ApplyPL(L, P, x, total_loads); % Construct SFD, BMD
SFD_PL = MakeSFD(total_loads)
BMD_PL = MakeBMD(SFD_PL)
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
%%
function [ total_loads, BMD_PL ] = ApplyPL( xP, P, x, total_loads )
% Constructs SFD and BMD from application of 1 Point Load. Assumes fixed location of supports
% Input: location and magnitude of point load. The previous SFD can be entered as input to construct SFD of multiple point loads
% Output: SFD, BMD both 1-D arrays of length n
% xP = distance from left (A) support (x = 0)
xP_B = 1060;
cur_P_B = -P * xP / xP_B; % P_B due to current P
P_A = total_loads(1) -P - cur_P_B; % Overall P_A
P_B = cur_P_B + total_loads(xP_B+1) - total_loads(xP_B); % Overall P_B

total_loads(1) = P_A;
total_loads(xP) = P;
total_loads(xP_B+1) = P_B;
BMD_PL = 0;
end

function [SFD_PL] = MakeSFD(total_loads)
    n = 1251;
    SFD_PL = zeros(1, n);
    SFD_PL(1) = total_loads(1);
    for i = 2:length(total_loads)
        SFD_PL(i) = total_loads(i) + SFD_PL(i-1);
    end
end

function [BMD_PL] = MakeBMD(SFD_PL)
    BMD_PL = cumsum(SFD_PL);
end