%% Problem Set 2
% Estimation of the Real Business Cycle model with Epstein-Zin
% Preferences, Productivity Shocks and Capital Share Shocks
% (5) Value Function Iteration with Accelerator
%
% Felipe Ruiz Mazin
% December 10, 2018

%% 0. Housekeeping

clear all
close all
clc

global inputs

tic

%%  1. Calibration

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
y = fsolve(@(x) sstate(x,laborSteadyState), [30 30]);

capitalSteadyState     = y(1);
eeta                   = y(2);
inputs.params(5)       = eeta; % eta
outputSteadyState      = capitalSteadyState^0.3 * laborSteadyState^0.7;
consumptionSteadyState = capitalSteadyState^0.3 * laborSteadyState^0.7 - (1 - (1 - inputs.params(2)))*capitalSteadyState;
utilitySteadyState     = log(capitalSteadyState^0.3 * laborSteadyState^0.7 - inputs.params(2)*capitalSteadyState) - y(2)*laborSteadyState^2/2;

fprintf(' Output = %2.6f, Capital = %2.6f, Consumption = %2.6f\n', outputSteadyState, capitalSteadyState, consumptionSteadyState); 
fprintf('\n')

inputs.simulParams = NaN(3,1);

% Generate the grid for capital
sizeCapital           = 250; % number of points in grid for capital
inputs.simulParams(1) = sizeCapital; 
capitalMin            = 0.5 * capitalSteadyState;
capitalMax            = 1.5 * capitalSteadyState;
gridCapital           = linspace(capitalMin, capitalMax, inputs.simulParams(1));
inputs.gridCapital    = gridCapital;

gridProd = [-0.0673; -0.0336; 0; 0.0336; 0.0673];
sizeProd = length(gridProd);
inputs.gridProd       = gridProd;
inputs.simulParams(2) = sizeProd;

gridAlpha = [0.25; 0.3; 0.35];
sizeAlpha = length(gridAlpha);
inputs.gridAlpha      = gridAlpha;
inputs.simulParams(3) = sizeAlpha;

sizeShocks = sizeProd*sizeAlpha;

%% 3. Required matrices and vectors

valueFunction0        = utilitySteadyState.*ones(sizeCapital,sizeShocks);
inputs.valueFunction0 = valueFunction0;
valueFunction         = ones(sizeCapital,sizeShocks);
policyFunction        = zeros(sizeCapital,sizeShocks);
laborFunction         = laborSteadyState.*ones(sizeCapital,sizeShocks);
consumptionFunction   = zeros(sizeCapital,sizeShocks);
inputs.laborFunction  = laborFunction;

policyFunction0       = repmat(interp1([capitalMin capitalMax],[capitalMin capitalMax],linspace(capitalMin, capitalMax, sizeCapital))',1,sizeShocks);%0.9*capitalSteadyState.*ones(sizeCapital,sizeShocks);

%% 4. Main iteration
maxIter = 10000;
tol     = 1e-6;

iter       = 1;
diff(iter) = 10;

options = optimset('Display', 'off');

% Finding value and policy functions numerically
while diff(iter) > tol && (iter <= maxIter)
    expectedValue0        = (valueFunction0.^ppsi) * transitionMatrix';
    inputs.expectedValue0 = expectedValue0;
    for indexShocks = 1:sizeShocks
        for indexCapital = 1:sizeCapital
            if mod(iter,10) == 1
                [policyFunction(indexCapital,indexShocks), vAux] = fminbnd(@(capitalPrime) ...
                    -value_function(capitalPrime,indexCapital,indexShocks),...
                    max(0.8*gridCapital(indexCapital),gridCapital(1)),min(1.2*gridCapital(indexCapital),gridCapital(end)),options);
                valueFunction(indexCapital,indexShocks) = -vAux;
                [~, laborFunction(indexCapital,indexShocks), consumptionFunction(indexCapital,indexShocks)] = ...
                    value_function(policyFunction(indexCapital,indexShocks),indexCapital,indexShocks);
            else
                indexProd  = ceil(indexShocks/sizeAlpha);
                indexAlpha = mod(indexShocks-1,sizeAlpha)+1;
                
                z            = gridProd(indexProd);
                aalpha       = gridAlpha(indexAlpha);
                capital      = gridCapital(indexCapital);
                capitalPrime = policyFunction(indexCapital,indexShocks);
                
                labor        = laborFunction(indexCapital,indexShocks);
                consumption  = consumptionFunction(indexCapital,indexShocks);
                
                indexCapitalLow  = max(sum(capitalPrime > gridCapital),1);
                indexCapitalHigh = indexCapitalLow + 1;

                expectedValue = expectedValue0(indexCapitalLow,indexShocks) + (expectedValue0(indexCapitalHigh,indexShocks) - ...
                    expectedValue0(indexCapitalLow,indexShocks))./(gridCapital(indexCapitalHigh) - gridCapital(indexCapitalLow)).*...
                    (capitalPrime - gridCapital(indexCapitalLow));
                valueFunction(indexCapital,indexShocks) =  ((1-bbeta)*(log(consumption) - eeta * labor^2/2)^rrho + bbeta * expectedValue^(rrho/ppsi))^(1/rrho);
                if isreal(valueFunction(indexCapital,indexShocks)) == 0
                    valueFunction(indexCapital,indexShocks) = valueFunction0(indexCapital,indexShocks);
                end
            end
        end
    end
    iter = iter + 1;
    if mod(iter,100) == 1
        fprintf(' Iteration: %2.0f, Sup diff: %2.6f\n', iter-1, diff(iter)); 
    end
    diff(iter) = max(abs(valueFunction - valueFunction0),[],'all');
    valueFunction0         = valueFunction;
    inputs.valueFunction0  = valueFunction0;
    policyFunction0        = policyFunction;
    inputs.laborFunction   = laborFunction;
end

toc

save PS2_accelerator2

%% (VALUE FUNCTIONS)

figure
j1 = [1 4 7 10 13];
for ii = 1:5
plot(gridCapital,valueFunction(:,j1(ii)))
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
saveas(gcf,'VFunc_025_acc.png')
figure
j2 = [2 5 8 11 14];
for ii = 1:5
plot(gridCapital,valueFunction(:,j2(ii)))
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
saveas(gcf,'VFunc_030_acc.png')
figure
j3 = [3 6 9 12 15];
for ii = 1:5
plot(gridCapital,valueFunction(:,j3(ii)))
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
saveas(gcf,'VFunc_035_acc.png')

%% (POLICY FUNCTIONS FOR CAPITAL)

figure
j1 = [1 4 7 10 13];
plot(gridCapital,gridCapital,'--k','LineWidth',0.7), hold on
for ii = 1:5
    plot(gridCapital,policyFunction(:,j1(ii)))
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
saveas(gcf,'KFunc_025_acc.png')
figure
j2 = [2 5 8 11 14];
plot(gridCapital,gridCapital,'--k','LineWidth',0.7), hold on
for ii = 1:5
    plot(gridCapital,policyFunction(:,j2(ii)))
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
saveas(gcf,'KFunc_030_acc.png')
figure
j3 = [3 6 9 12 15];
plot(gridCapital,gridCapital,'--k','LineWidth',0.7), hold on
for ii = 1:5
    plot(gridCapital,policyFunction(:,j3(ii)))
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
saveas(gcf,'KFunc_035_acc.png')

%% (POLICY FUNCTIONS FOR LABOR)

figure
j1 = [1 4 7 10 13];
for ii = 1:5
    plot(gridCapital,laborFunction(:,j1(ii)))
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
saveas(gcf,'LFunc_025_acc.png')
figure
j2 = [2 5 8 11 14];
for ii = 1:5
    plot(gridCapital,laborFunction(:,j2(ii)))
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
saveas(gcf,'LFunc_030_acc.png')
figure
j3 = [3 6 9 12 15];
for ii = 1:5
    plot(gridCapital,laborFunction(:,j3(ii)))
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
saveas(gcf,'LFunc_035_acc.png')

%% (POLICY FUNCTIONS FOR CONSUMPTION)

figure
j1 = [1 4 7 10 13];
for ii = 1:5
    plot(gridCapital,consumptionFunction(:,j1(ii)))
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
saveas(gcf,'CFunc_025_acc.png')
figure
j2 = [2 5 8 11 14];
for ii = 1:5
    plot(gridCapital,consumptionFunction(:,j2(ii)))
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
saveas(gcf,'CFunc_030_acc.png')
figure
j3 = [3 6 9 12 15];
for ii = 1:5
    plot(gridCapital,consumptionFunction(:,j3(ii)))
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
saveas(gcf,'CFunc_035_acc.png')

%% Calculate impulse response

%load PS2
T = 40;

% Find ergodic steady state when z=0,alpha=0.3
position=policyFunction(:,8)>gridCapital'; % find first time policy function crosses the 45° line!
index = find(position==1,1,'last');

% Initial position in ergodic steady-state
% one time transitory shock to productivity
impulseRespCapitalProd = zeros(1,T);
impulseRespLaborProd   = zeros(1,T);
impulseRespConsProd    = zeros(1,T);

impulseRespCapitalProd(1) = policyFunction(index,8);
impulseRespLaborProd(1)   = laborFunction(index,8);
impulseRespConsProd(1)    = consumptionFunction(index,8);
impulseRespCapitalProd(2) = policyFunction(index,14);
impulseRespLaborProd(2)   = laborFunction(index,14);
impulseRespConsProd(2)    = consumptionFunction(index,14);

% Interpolate
for t = 3:T
    
    capitalLowIR = max(sum(impulseRespCapitalProd(t-1)>gridCapital));
    capitalHighIR = capitalLowIR+1;
    
    impulseRespCapitalProd(t) = interp1([gridCapital(capitalLowIR),...
        gridCapital(capitalHighIR)],[policyFunction(capitalLowIR,8),...
        policyFunction(capitalHighIR,8)],impulseRespCapitalProd(t-1));
    
    impulseRespLaborProd(t)  = interp1([gridCapital(capitalLowIR),...
        gridCapital(capitalHighIR)],[laborFunction(capitalLowIR,8),...
        laborFunction(capitalHighIR,8)],impulseRespCapitalProd(t-1));
    
    impulseRespConsProd(t)  = interp1([gridCapital(capitalLowIR),...
        gridCapital(capitalHighIR)],[consumptionFunction(capitalLowIR,8),...
        consumptionFunction(capitalHighIR,8)],impulseRespCapitalProd(t-1));    
end

    
% Change to percentage deviations from ss
impulseRespCapitalProdPD = zeros(1,T);
impulseRespLaborProdPD   = zeros(1,T);
impulseRespConsProdPD    = zeros(1,T);
    
for t = 1:T
    impulseRespCapitalProdPD(t) = (log(impulseRespCapitalProd(t))...
        -log(impulseRespCapitalProd(1)))*100;
    impulseRespLaborProdPD(t) = (log(impulseRespLaborProd(t))...
        -log(impulseRespLaborProd(1)))*100;
    impulseRespConsProdPD(t) = (log(impulseRespConsProd(t))...
        -log(impulseRespConsProd(1)))*100;
end

% Initial position in ergodic steady-state
% one time transitory shock to capital share
impulseRespCapitalShare = zeros(1,T);
impulseRespLaborShare = zeros(1,T);
impulseRespConsShare = zeros(1,T);

impulseRespCapitalShare(1)  = policyFunction(index,8);
impulseRespLaborShare(1)    = laborFunction(index,8);
impulseRespConsShare(1)     = consumptionFunction(index,8);
impulseRespCapitalShare(2)  = policyFunction(index,9);
impulseRespLaborShare(2)    = laborFunction(index,9);
impulseRespConsShare(2)     = consumptionFunction(index,9);

% Interpolate
for t = 3:T
    
    capitalLowIR = max(sum(impulseRespCapitalShare(t-1)>gridCapital));
    capitalHighIR = capitalLowIR+1;
    
    impulseRespCapitalShare(t) = interp1([gridCapital(capitalLowIR),...
        gridCapital(capitalHighIR)],[policyFunction(capitalLowIR,8),...
        policyFunction(capitalHighIR,8)],impulseRespCapitalShare(t-1));
    
    impulseRespLaborShare(t)  = interp1([gridCapital(capitalLowIR),...
        gridCapital(capitalHighIR)],[laborFunction(capitalLowIR,8),...
        laborFunction(capitalHighIR,8)],impulseRespCapitalShare(t-1));
    
    impulseRespConsShare(t)  = interp1([gridCapital(capitalLowIR),...
        gridCapital(capitalHighIR)],[consumptionFunction(capitalLowIR,8),...
        consumptionFunction(capitalHighIR,8)],impulseRespCapitalShare(t-1));    
end

% Change to percentage deviations from ss
impulseRespCapitalSharePD = zeros(1,T);
impulseRespLaborSharePD   = zeros(1,T);
impulseRespConsSharePD    = zeros(1,T);
    
for t = 1:T
    impulseRespCapitalSharePD(t) = (log(impulseRespCapitalShare(t))...
        -log(impulseRespCapitalShare(1)))*100;
    impulseRespLaborSharePD(t) = (log(impulseRespLaborShare(t))...
        -log(impulseRespLaborShare(1)))*100;
    impulseRespConsSharePD(t) = (log(impulseRespConsShare(t))...
        -log(impulseRespConsShare(1)))*100;
end

%% Plot impulse response functions

zPath = zeros(1,T)+0;
zPath(2) = gridProd(5);
alphaPath = zeros(1,T)+0.3;
alphaPath(2) = gridAlpha(3);

figure
subplot(2,2,1)
plot(impulseRespCapitalProdPD)
hold on
plot(zeros(1,T),'k')
title('Capital')
%xlabel('time')
ylabel('% deviation from ss')
subplot(2,2,2)
plot(impulseRespLaborProdPD)
title('Labor')
%xlabel('time')
ylabel('% deviation from ss')
hold on
plot(zeros(1,T),'k')
subplot(2,2,3)
plot(impulseRespConsProdPD)
title('Consumption')
xlabel('time')
ylabel('% deviation from ss')
hold on
plot(zeros(1,T),'k')
subplot(2,2,4)
plot(zPath)
hold on
plot(zeros(1,T)+0,'k')
title('Productivity')
xlabel('time')
ylabel('z (level)')
saveas(gcf,'IRFz_acc.png')

figure
subplot(2,2,1)
plot(impulseRespCapitalSharePD)
hold on
plot(zeros(1,T),'k')
title('Capital')
%xlabel('time')
ylabel('% deviation from ss')
subplot(2,2,2)
plot(impulseRespLaborSharePD)
hold on
plot(zeros(1,T),'k')
title('Labor')
%xlabel('time')
ylabel('% deviation from ss')
subplot(2,2,3)
plot(impulseRespConsSharePD)
hold on
plot(zeros(1,T),'k')
title('Consumption')
xlabel('time')
ylabel('% deviation from ss')
subplot(2,2,4)
plot(alphaPath)
hold on
plot(zeros(1,T)+0.3,'k')
title('Capital Share')
xlabel('time')
ylabel('\alpha (level)')
saveas(gcf,'IRFalpha_acc.png')

%% Convergence
figure
plot(log10(diff(2:end)))
xlim([1, length(diff)])
tit=title('Log10 Sup Difference');
set(tit,'FontSize',14);
ax  = gca;
set(gca,'FontSize',14)
saveas(gcf,'error_acc.png')