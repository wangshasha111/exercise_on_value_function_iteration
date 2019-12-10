%% Problem Set 2
% Estimation of the Real Business Cycle model with Epstein-Zin
% Preferences, Productivity Shocks and Capital Share Shocks
% (7) Spectral Method using Chebyshev Polynomials
%
% Felipe Ruiz Mazin
% December 10, 2018

%% 0. Housekeeping

clear all
close all
clc

tic

%%  1. Calibration

inputs.paramNames = 'beta';
inputs.paramNames = char(inputs.paramNames, 'delta');
inputs.paramNames = char(inputs.paramNames, 'psi');
inputs.paramNames = char(inputs.paramNames, 'rho');
inputs.paramNames = char(inputs.paramNames, 'eta');

bbeta  = 0.99;  % Discount factor
ddelta = 0.1;  % Depreciation rate
ppsi   = -9;     % Risk aversion factor
rrho   = 0.5;    % Elasticity of intertemporal substitution factor: 1/(1-rrho)

inputs.params     = NaN(5,1);
inputs.params(1)  = bbeta;   % Discount factor
inputs.params(2)  = ddelta;  % Depreciation rate
inputs.params(3)  = ppsi;    % Risk aversion factor
inputs.params(4)  = rrho;    % Elasticity of intertemporal substitution factor: 1/(1-rrho)


% Transition matrices
transitionProd  = [0.9727 0.0273 0 0 0; ...
                   0.0041 0.9806 0.0153 0 0; ...
                   0 0.0082 0.9836 0.0082 0; ...
                   0 0 0.0153 0.9806 0.0041; ...
                   0 0 0 0.0273 0.9727];

transitionAlpha = [0.9 0.07 0.03; ...
                   0.05 0.9 0.05; ...
                   0.03 0.07 0.9];
               
% Collapsing both states into only one with size 15
% Transition matrix:

transitionMatrix        = kron(transitionProd,transitionAlpha);
inputs.transitionMatrix = transitionMatrix;

%% 2. Steady State

laborSteadyState = 50;
y = fsolve(@(x) sstate(x,laborSteadyState), [4 5]);

capitalSteadyState     = y(1);
eeta                   = y(2);
inputs.params(5)       = eeta; % eta
outputSteadyState      = capitalSteadyState^0.3 * laborSteadyState^0.7;
consumptionSteadyState = capitalSteadyState^0.3 * laborSteadyState^0.7 - (1 - (1 - inputs.params(2)))*capitalSteadyState;
utilitySteadyState     = log(capitalSteadyState^0.3 * laborSteadyState^0.7 - inputs.params(2)*capitalSteadyState) - y(2)*laborSteadyState^2/2;
inputs.steadyState     = capitalSteadyState;

fprintf(' Output = %2.6f, Capital = %2.6f, Consumption = %2.6f\n', outputSteadyState, capitalSteadyState, consumptionSteadyState); 
fprintf('\n')

inputs.simulParams = NaN(5,1);

gridProd = [-0.0673; -0.0336; 0; 0.0336; 0.0673];
sizeProd = length(gridProd);
inputs.gridProd       = gridProd;
inputs.simulParams(2) = sizeProd;

gridAlpha = [0.25; 0.3; 0.35];
sizeAlpha = length(gridAlpha);
inputs.gridAlpha      = gridAlpha;
inputs.simulParams(3) = sizeAlpha;

sizeShocks = sizeProd*sizeAlpha;

% Euler errors
eulerSize  = 3000;                  % # of grid points for  capital (to compute euler errors)
inputs.eulerSize = eulerSize;

% Simulation
T         = 10000;                 % Number of periods for the simulation of the economy
dropT     = 1000;                  % Burn-in
inputs.simulParams(4) = T;
inputs.simulParams(5) = dropT;

%% 3. Spectral Method using Chebychev Polynomials

% Define boundaries for capital
coverGrid  = 0.5;
capitalMin = capitalSteadyState*(1-coverGrid);
capitalMax = capitalSteadyState*(1+coverGrid);
interval   = capitalSteadyState*2*coverGrid;
inputs.simulParams(6) = capitalMin;
inputs.simulParams(7) = capitalMax;
inputs.interval   = interval;

tic

nodeNum = 8;
inputs.simulParams(1) = nodeNum;
M       = nodeNum*sizeProd*sizeAlpha;

% Find Zeros of the Chebyshev Polynomial of order 8 
zerosCheby = -cos((2*(1:nodeNum)'-1)*pi/(2*nodeNum));

% Define Chebyshev polynomials
chebyPoly      = ones(nodeNum,nodeNum);
chebyPoly(:,2) = zerosCheby;

for i1 = 3:nodeNum
    chebyPoly(:,i1) = 2*zerosCheby.*chebyPoly(:,i1-1)-chebyPoly(:,i1-2);
end
inputs.chebyPoly = chebyPoly;

% Project collocation points in the K space
gridCapital           = ((zerosCheby+1)*(capitalMax-capitalMin))/2+capitalMin;
inputs.gridCapital    = gridCapital;

% Initial Guess for Chebyshev coefficients
thetaGuess = zeros(M,2);

for indexShock = 1:sizeProd*sizeAlpha
    thetaGuess((indexShock-1)*nodeNum+1,1)   = utilitySteadyState;
    thetaGuess((indexShock-1)*nodeNum+1,2) = laborSteadyState;
end

% Solve for Chebyshev coefficients
theta = residual_function(thetaGuess,inputs);

toc

%----------------------------------------------------------------
% 3. Compute Euler Errors and Decision rules
%----------------------------------------------------------------
%%
gridCapitalComplete = zeros(eulerSize,1);

for i = 1:eulerSize
    gridCapitalComplete(i) = capitalMin + (i-1)*interval/(eulerSize-1);
end

[policyFunction,consumptionFunction,laborFunction ,valueFunction,eulerError,maxError]= ...
         eulerr_grid(theta,inputs);
     
[vSeries,kSeries,cSeries,lSeries,ySeries,eeSeries,rbSeries,rkSeries,rkCondSeries]= ...
         simulation(theta,inputs);

mean_error    = sum(eeSeries)/(T-dropT);
max_error_sim = max(eeSeries);

disp(' ')
disp('Integral of Euler Equation Error:')
disp(mean_error)
disp('Max Euler Equation Error Simulation:')
disp(max_error_sim)

%----------------------------------------------------------------
% 4. Figures
%----------------------------------------------------------------

%% Decision Rules
% (VALUE FUNCTIONS)

figure
j1 = [1 4 7 10 13];
for ii = 1:5
plot(gridCapitalComplete,valueFunction(:,j1(ii)))
hold on
end
xlim([capitalMin, capitalMax])
ylim([3.46, 3.62])
xlabel('Capital')
ylabel('Value')
tit=title('Value Function (\alpha = 0.25)');
leg=legend(strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','southeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'VFunc_025_Cheby.png')
figure
j2 = [2 5 8 11 14];
for ii = 1:5
plot(gridCapitalComplete,valueFunction(:,j2(ii)))
hold on
end
xlim([capitalMin, capitalMax])
ylim([3.46, 3.62])
xlabel('Capital')
ylabel('Value')
tit=title('Value Function (\alpha = 0.30)');
leg=legend(strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','southeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'VFunc_030_Cheby.png')
figure
j3 = [3 6 9 12 15];
for ii = 1:5
plot(gridCapitalComplete,valueFunction(:,j3(ii)))
hold on
end
xlim([capitalMin, capitalMax])
ylim([3.46, 3.62])
xlabel('Capital')
ylabel('Value')
tit=title('Value Function (\alpha = 0.35)');
leg=legend(strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','southeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'VFunc_035_Cheby.png')

%% (POLICY FUNCTIONS FOR CAPITAL)

figure
j1 = [1 4 7 10 13];
plot(gridCapitalComplete,gridCapitalComplete,'--k','LineWidth',0.7), hold on
for ii = 1:5
    plot(gridCapitalComplete,policyFunction(:,j1(ii)))
    hold on
end
xlim([capitalMin, capitalMax])
ylim([capitalMin, capitalMax])
xlabel('Current Capital')
ylabel('Capital Next Period')
tit=title('Policy Function for Capital (\alpha = 0.25)');
leg=legend('45° line',strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','southeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'KFunc_025_Cheby.png')
figure
j2 = [2 5 8 11 14];
plot(gridCapitalComplete,gridCapitalComplete,'--k','LineWidth',0.7), hold on
for ii = 1:5
    plot(gridCapitalComplete,policyFunction(:,j2(ii)))
    hold on
end
xlim([capitalMin, capitalMax])
ylim([capitalMin, capitalMax])
xlabel('Current Capital')
ylabel('Capital Next Period')
tit=title('Policy Function for Capital (\alpha = 0.30)');
leg=legend('45° line',strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','southeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'KFunc_030_Cheby.png')
figure
j3 = [3 6 9 12 15];
plot(gridCapitalComplete,gridCapitalComplete,'--k','LineWidth',0.7), hold on
for ii = 1:5
    plot(gridCapitalComplete,policyFunction(:,j3(ii)))
    hold on
end
xlim([capitalMin, capitalMax])
ylim([capitalMin, capitalMax])
xlabel('Current Capital')
ylabel('Capital Next Period')
tit=title('Policy Function for Capital (\alpha = 0.35)');
leg=legend('45° line',strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','southeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'KFunc_035_Cheby.png')

%% (POLICY FUNCTIONS FOR LABOR)

figure
j1 = [1 4 7 10 13];
for ii = 1:5
    plot(gridCapitalComplete,laborFunction(:,j1(ii)))
    hold on
end
xlim([capitalMin, capitalMax])
ylim([45, 55])
xlabel('Capital')
ylabel('Labor')
tit=title('Policy Function for Labor (\alpha = 0.25)');
leg=legend(strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','northeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'LFunc_025_Cheby.png')
figure
j2 = [2 5 8 11 14];
for ii = 1:5
    plot(gridCapitalComplete,laborFunction(:,j2(ii)))
    hold on
end
xlim([capitalMin, capitalMax])
ylim([45, 55])
xlabel('Capital')
ylabel('Labor')
tit=title('Policy Function for Labor (\alpha = 0.30)');
leg=legend(strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','northeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'LFunc_030_Cheby.png')
figure
j3 = [3 6 9 12 15];
for ii = 1:5
    plot(gridCapitalComplete,laborFunction(:,j3(ii)))
    hold on
end
xlim([capitalMin, capitalMax])
ylim([45, 55])
xlabel('Capital')
ylabel('Labor')
tit=title('Policy Function for Labor (\alpha = 0.35)');
leg=legend(strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','northeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'LFunc_035_Cheby.png')

%% (POLICY FUNCTIONS FOR CONSUMPTION)

figure
j1 = [1 4 7 10 13];
for ii = 1:5
    plot(gridCapitalComplete,consumptionFunction(:,j1(ii)))
    hold on
end
xlim([capitalMin, capitalMax])
ylim([35, 75])
xlabel('Capital')
ylabel('Consumption')
tit=title('Policy Function for Consumption (\alpha = 0.25)');
leg=legend(strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','southeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'CFunc_025_Cheby.png')
figure
j2 = [2 5 8 11 14];
for ii = 1:5
    plot(gridCapitalComplete,consumptionFunction(:,j2(ii)))
    hold on
end
xlim([capitalMin, capitalMax])
ylim([35, 75])
xlabel('Capital')
ylabel('Consumption')
tit=title('Policy Function for Consumption (\alpha = 0.30)');
leg=legend(strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','southeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'CFunc_030_Cheby.png')
figure
j3 = [3 6 9 12 15];
for ii = 1:5
    plot(gridCapitalComplete,consumptionFunction(:,j3(ii)))
    hold on
end
xlim([capitalMin, capitalMax])
ylim([35, 75])
xlabel('Capital')
ylabel('Consumption')
tit=title('Policy Function for Consumption (\alpha = 0.35)');
leg=legend(strcat('z = ', num2str(round(gridProd(1),2))),strcat('z = ', num2str(round(gridProd(2),2))),...
strcat('z = ', num2str(round(gridProd(3),2))),strcat('z = ', num2str(round(gridProd(4),2))),strcat('z = ', num2str(round(gridProd(5),2))),'location','southeast');
set(leg,'FontSize',14);
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'CFunc_035_Cheby.png')

%% Euler Equation Error on the Grid
figure(2)
plot(gridCapitalComplete,eulerError)
title('Log10 Euler Error')
xlim([capitalMin, capitalMax])
xlabel('Capital')
saveas(gcf,'Euler_errors.png')

%% Distribution of simulated variables

[f_c,x_c]   = ksdensity(cSeries);
[f_k,x_k]   = ksdensity(kSeries);
[f_y,x_y]   = ksdensity(ySeries);
[f_l,x_l]   = ksdensity(lSeries);
[f_rb,x_rb] = ksdensity(rbSeries);
[f_rk,x_rk] = ksdensity(rkSeries);

figure(3)
subplot(3,2,1)
plot(x_c,f_c)
title('Density of Consumption')
subplot(3,2,2)
plot(x_l,f_l)
title('Density of Labor')
subplot(3,2,3)
plot(x_k,f_k)
title('Density of Capital')
subplot(3,2,4)
plot(x_y,f_y)
title('Density of Output')
subplot(3,2,5)
plot(x_rb,f_rb)
title('Density of Return on Risk Free bond')
saveas(gcf,'densities.png')
subplot(3,2,6)
plot(x_rk,f_rk)
title('Density of Return on Equity')
saveas(gcf,'densities.png')