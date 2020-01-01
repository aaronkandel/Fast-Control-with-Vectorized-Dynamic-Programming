% Vectorized Dynamic Programming for Nonlinear Control
% Aaron Kandel

%{
In this script, I implement a vectorized dynamic programming procedure to solve
the optimal fast-charging problem for a lithium-ion battery. I use a simple
equivalent circuit model of a lithium-ion battery with a single RC pair and
nonlinear output voltage equation.  The principal constraint is on the
output voltage remaining within a safe operating range.
%}


clc
clear
close all

%% Initialize Battery Parameters:

% Load and Fit OCV Data:
VOC_data = csvread('Voc.dat',1,0);
soc = VOC_data(:,1);
voc = VOC_data(:,2);
% Fit cubic polynomial to OCV:
Reg = [ones(length(soc),1),soc, soc.^2, soc.^3];
theta = inv(Reg'*Reg)*Reg'*voc;

% Charging Parameters:
tmax = 250; % [s] max time
dt = 0.5; % [s] timestep
tint = 0:dt:tmax; % [s] time vector
I_max = 36; % [A] maximum control input
I_min = 0; % [A] minimum current
C_Batt = 2.3*3600;% [A-hr]
R_0 = 0.01; % [Ohm]
R_1 = 0.01; % [Ohm]
C_1 = 2500; % [F]
dt = 1; % [s]

% State Transition Matrices:
A = [1, 0;0, (1-(dt/(R_1(1)*C_1(1))))]; 
B = [dt/C_Batt(1);dt/C_1(1)];
C = [theta',1];
D = R_0(1);

% Charging Parameters:
SOC0 = 0.2; % [-] Initial State-of-Charge
SOC_targ = 0.5; % [-] Target State-of-Charge
VsimLim = 3.6; % [V] Voltage Constraint Boundary

% Define State and Control Variable Meshes:
x1mesh = (SOC0:0.0025:0.55)'; % [-] SOC
x2mesh = (0:0.005:0.2)'; % [V] V_{RC}
umesh = I_min:0.5:I_max; % [A] Input Current 

% Terminal State Index:
xf_index = find(x1mesh==SOC_targ);

% Define Big Arrays
cost2Go = 1e5*ones(length(x1mesh)*length(x2mesh),length(umesh)); %);
csoc = zeros(length(x1mesh)*length(x2mesh),length(umesh));
nsoc = zeros(length(x1mesh)*length(x2mesh),length(umesh));
tCost = zeros(length(x1mesh)*length(x2mesh),length(umesh));

% Collect Discrete States:
x1a = repmat(x1mesh, length(x2mesh), 1);
x1a = sort(x1a);
x2a = repmat(x2mesh, length(x1mesh), 1);
x = [x1a, x2a];



%% Obtain Transition Costs:

xvec = x';
dx1 = (x1mesh(2) - x1mesh(1));
dx2 = (x2mesh(3) - x2mesh(2));

for u = 1:length(umesh) % Iterate over control inputs
    
    % Compute state transitions, nearest state indices:
    nxvec = (A*xvec + B*umesh(u))'; 
    xyvec = [ones(length(x), 1), nxvec(:, 1), nxvec(:, 1).^2, nxvec(:, 1).^3, nxvec(:,2)];
    nyv = C*xyvec' + D*umesh(u);
    nxv1 = round(nxvec(:,1)./dx1) .* dx1;
    nxv2 = round(nxvec(:,2)./dx2) .* dx2;
    nxInd1 = (nxv1 - min(x1mesh))./dx1 + 1;
    nxInd2 = (nxv2 - min(x2mesh))./dx2 + 1;
    nxInd = [nxInd1, nxInd2];
    nxInd = round(nxInd);
    
    % Compute next state index:
    nxCol(:,u) = ((nxInd(:,1)-1)*length(x2mesh) + (nxInd(:,2))); 
    
    % Compute State Transition Costs:
    tCost(:,u) = (nxvec(:,1) - SOC_targ).^2 + 100*(nyv' >= VsimLim);
    % Transition cost depends on SOC reference tracking and an indicator
    % function on the voltage constraint.
end % END FOR

% Set transition cost at terminal point to zero:
indx = find(xvec(1,:) >= SOC_targ);
cost2Go(indx, 1) = 0;



%% Initialize Vectorization for DP Code:

% Maintain feasible state indices:
N = (tmax-0)/dt;
newXind = nxCol;
rch = find(newXind > length(x1mesh)*length(x2mesh));
newXind(rch) = length(x1mesh)*length(x2mesh);
[cost2GoNow, minctg] = min(cost2Go,[],2); 


% Initialize cost-to-go, next-cost, and best input at final timestep:
[rowC2G,colC2G] = find(cost2Go == 0);
nextCostNow = cost2GoNow;
nextCost = cost2Go;
newXind = round(newXind);
nextCost(nxInd) = nextCostNow(newXind(nxInd));

ctg = 1e5*ones(length(x1mesh)*length(x2mesh),length((N:-1:1))+1);
ctg(:,end) = cost2GoNow;

bestU = zeros(length(x1mesh)*length(x2mesh),length((N:-1:1))+1);
bu = zeros(length(x1mesh)*length(x2mesh),1);
bu(rowC2G) = umesh(colC2G);
bestU(:,end) = bu;

%% Backwards Recursion to Compute Control:

nx = length(x1mesh)*length(x2mesh);
nu = length(umesh);

% Add algorithm function folder:
cd ..
addpath Algorithm

bestU = DP_Loop(N, nx, nu, newXind, nextCost, tCost, cost2GoNow);

%% Simulate Results:

tsim = 250;
SOC = 0.2;
Vrc = 0;
VOC = C*[1, SOC, SOC^2, SOC^3, 0]';
Vsim = VOC;
xsim = [SOC;Vrc];
ysim = Vsim;
usim = zeros(1, tsim);

for i = 1:tsim
    % Find Nearest State:
    xtest = xsim(:,i);
    xnorm = xtest./[dx1;dx2];
    xnorm = round(xnorm).*[dx1;dx2];
    nx = (xnorm - [min(x1mesh);min(x2mesh)])./[dx1;dx2] + 1;
    xind = ((nx(1)-1)*length(x2mesh) + (nx(2))); % *length(x1mesh) + (nxInd(:,3))
    xind  = round(xind);
    usim(i) = umesh(bestU(xind, i));
    
    xsim(:, i+1) = A*xsim(:,i) + B*usim(i);
    cvec = [1, xsim(1,i), xsim(1,i)^2, xsim(1,i)^3, xsim(2,i)];
    ysim(i) = C*cvec' + D*usim(i);
    
end % END FOR

tvec = 0:dt:tsim;
figure
subplot(231)
plot(tvec, xsim(1,:))
grid on
xlabel('Time [s]')
ylabel('State of Charge [-]')
subplot(232)
plot(tvec, xsim(2,:))
grid on
xlabel('Time [s]')
ylabel('Capacitor Voltage [V]')
subplot(233)
plot(tvec(1:end-1), usim)
grid on
xlabel('Time [s]')
ylabel('Input Current [A]')
subplot(2,3,[4,5,6])
hold on
plot(tvec(1:end-1), ysim)
plot([0, tsim], [VsimLim, VsimLim], '--r', 'Linewidth', 2)
legend('Voltage', 'Constraint')
grid on
xlabel('Time [s]')
ylabel('Voltage [V]')








