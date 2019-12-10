%% Problem Set 2
% Estimation of the Real Business Cycle model with Epstein-Zin
% Preferences, Productivity Shocks and Capital Share Shocks
% (4) Value Function Iteration with an Endogenous Grid
%
% Felipe Ruiz Mazin
% December 10, 2018

%% 0. Housekeeping

clear all
close all
clc

tic

%% 1. Calibration

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
for indexCapital = 1:sizeCapital
    valueFunction0(indexCapital,:) = (0.5*utilitySteadyState + (indexCapital - sizeCapital/2)/sizeCapital*utilitySteadyState)...
        .*ones(1,sizeShocks);
end
inputs.valueFunction0 = valueFunction0;
policyFunction        = zeros(sizeCapital,sizeShocks);
laborFunction         = laborSteadyState.*ones(sizeCapital,sizeShocks);
consumptionFunction   = zeros(sizeCapital,sizeShocks);
inputs.laborFunction  = laborFunction;

policyFunction0       = 0.8*capitalSteadyState.*ones(sizeCapital,sizeShocks);


%% 4. Endogenous Grid Method 1
% Implement EGM with labor fixed to steady state value until convergence
% is achieved

%Setting up Vtilde at n=0
valueFunctionTilde0 = bbeta*(valueFunction0).^rrho;

% Preallocate
endogenousValueFunction   = zeros(sizeCapital,sizeShocks);
endogenousMarketResources = zeros(sizeCapital,sizeShocks);
valueFunction             = zeros(sizeCapital,sizeShocks);
valueFunctionTilde        = zeros(sizeCapital,sizeShocks);
endogenousConsumption     = zeros(sizeCapital,sizeShocks);

% Setting the Market Resources
gridOutput = zeros(sizeCapital,sizeShocks);
for indexShocks = 1:sizeShocks
    indexProd  = ceil(indexShocks/sizeAlpha);
    indexAlpha = mod(indexShocks-1,sizeAlpha)+1;

    aalpha     = gridAlpha(indexAlpha);
    z          = gridProd(indexProd);

    gridOutput(:,indexShocks) = exp(z)*(gridCapital'.^aalpha*laborSteadyState.^(1-aalpha))+(1-ddelta)*gridCapital';
end

% Tolerance and max iterations
maxIter = 10000;
tol     = 1e-6;

diff     = 10;
iter     = 1;

tic
while diff > tol && iter <= maxIter
        for indexShocks = 1:sizeShocks
             for indexCapital = 1:sizeCapital
                 
% Take the derivative of ValueFunctionTilde
% Use average of the slopes of the linearly interpolated value function
        
                if indexCapital == sizeCapital
                    
                    derivativeValueFunctionTilde = (valueFunctionTilde0(indexCapital,indexShocks)-...
                    valueFunctionTilde0(indexCapital-1,indexShocks))/...
                    (gridCapital(indexCapital)-gridCapital(indexCapital-1));
                
                elseif indexCapital == 1
                    
                    derivativeValueFunctionTilde = (valueFunctionTilde0(indexCapital+1,indexShocks)-...
                    valueFunctionTilde0(indexCapital,indexShocks))/...
                    (gridCapital(indexCapital+1)-gridCapital(indexCapital));
                                                   
                else
                    
                    derivativeValueFunctionTilde = ((valueFunctionTilde0(indexCapital,indexShocks)...
                        -valueFunctionTilde0(indexCapital-1,indexShocks))/...
                       (gridCapital(indexCapital)-gridCapital(indexCapital-1))+...
                       (valueFunctionTilde0(indexCapital+1,indexShocks)-...
                       valueFunctionTilde0(indexCapital,indexShocks))/...
                       (gridCapital(indexCapital+1)-gridCapital(indexCapital)))/2;
                end
                
% Compute the optimal level of consumption 

                consumption_function = @(x) x - rrho*(1-bbeta)*(log(x)-eeta*(laborSteadyState^2/2))^(rrho-1)*...
                    derivativeValueFunctionTilde^(-1);
                consumption = fzero(consumption_function, [20 1000]); 
                
% Compute the value of the endogenously determined market resources
                endogenousMarketResources(indexCapital,indexShocks)=consumption+gridCapital(indexCapital);

% Update the value function
                endogenousValueFunction(indexCapital,indexShocks) =...
                 ((1-beta)*(log(Consumption) - eta*(laborSteadyState^2)/2)^(1-rrho) ...
                 + valueFunctionTilde0(indexCapital,indexShocks))^(1/(1-rrho));
             
                   if isreal(endogenousMarketResources(indexCapital,indexShocks)) == 0|| ...
                           isreal(endogenousValueFunction(indexCapital,indexShocks)) ==0
                        break;
                   else
%Interpolating so that we define the value function on the
% grid for capital established before.
                
                    for indexOutput = 1:sizeCapital
                        outputLow = max(sum(gridOutput(indexOutput,indexShocks)>endogenousMarketResources(:,indexShocks)),1);
                        if outputLow == sizeCapital
                            outputLow = sizeCapital-1;
                        end
                        outputHigh = outputLow + 1;
                        valueFunction(indexOutput,indexShocks) = endogenousValueFunction(outputLow,indexShocks)...
                            +(gridOutput(indexOutput,indexShocks)-endogenousMarketResources(outputLow,indexShocks))...
                            *(endogenousValueFunction(outputHigh,indexShocks)-endogenousValueFunction(outputLow,indexShocks))...
                            /(endogenousMarketResources(outputHigh,indexShocks)-endogenousMarketResources(outputLow,indexShocks));
                    end           
                   end
             end % end capital
        end % end shocks
        
% Calculate Vtilde in t+1
        valueFunctionTilde0New = beta*(valueFunction.^ppsi*transitionMatrix').^(rrho/ppsi);
% Criteria for convergence
        diff = max(abs(valueFunctionTilde0New(:)-valueFunctionTilde0(:)));
        valueFunctionTilde0 = valueFunctionTilde0New;
    
        if (mod(iter,10)==0 || iter ==1)
            fprintf(' Iteration EGM 1 = %d, Sup Diff = %2.8f\n', iter, diff); 
        end
        iter = iter+1;
end
fprintf(' Total iterations EGM 1 = %d\n', iter-1)

%% Retrieve K_endogenous for all values of the shocks
endogenousCapitalGrid = zeros(sizeCapital,sizeShocks);

for indexShocks = 1:sizeShocks
    indexProd  = ceil(indexShocks/sizeAlpha);
    indexAlpha = mod(indexShocks-1,sizeAlpha)+1;

    capital    = gridCapital(indexCapital);
    aalpha     = gridAlpha(indexAlpha);
    z          = gridProd(indexProd);
        
    capital0 = capitalSteadyState;
    for indexCapital = 1:sizeCapital

        KendEGM1 = @(x) endogenousMarketResources(indexCapital,indexShocks)...
        -exp(z)*x^aalpha*laborSteadyState^(1-aalpha)-(1-ddelta)*x;

        endogenousCapitalGrid(indexCapital,indexShocks) = fsolve(KendEGM1,capitalSteadyState,options);

    end
end

%% Interpolated value function on the capital grid to do VFI

valueFunctionInterpolatedVFI1 = zeros(sizeCapital,sizeShocks);

for indexShocks = 1:sizeShocks
        
    for indexCapital = 1:sizeCapital

        capitalLow = max(sum(gridCapital(indexCapital)>endogenousCapitalGrid(:,indexShocks)),1);

        if capitalLow == sizeCapital
           capitalLow = sizeCapital-1;
        end

        capitalHigh = capitalLow + 1;

        valueFunctionInterpolatedVFI1(indexCapital,indexShocks) = valueFunction(capitalLow,indexShocks)...
            +(gridCapital(indexCapital)-endogenousCapitalGrid(capitalLow,indexShocks))...
            *(valueFunction(capitalHigh,indexShocks)-valueFunction(capitalLow,indexShocks))...
            /(endogenousCapitalGrid(capitalHigh,indexShocks)-endogenousCapitalGrid(capitalLow,indexShocks));

    end
end

%% 6. Generalized EGM
valueFunction0 = valueFunctionInterpolatedVFI1;

tol     = 1e-06;
maxIter = 3;
diff = 1;
tic
while diff > tol
    %% 6a. VFI1 algorithm to recover the policy functions.

    [valueFunction, policyCapital, policyLabor, policyConsumption, diff] = ...
        VFIinterpolation(valueFunction0,gridCapital,...
        vShocks,transitionMatrix,capitalMin,capitalMax,40*ones(1,sizeShocks),...
    400*ones(1,sizeShocks),functions,parameters,tol,maxIter);

    fprintf('\n Error = %2.8f\n', diff)
   
    
    %% 6b. EGM algorithm to get next value function guess
if diff > 10^(-5)
    valueFunctionEGM = EGM(parameters,gridCapital,vShocks,...
        mTransition,c_ss,k_ss,policyCapital,policyLabor,valueFunction0,10*tol,options);
else
        valueFunctionEGM = EGM(parameters,gridCapital,vShocks,...
        mTransition,c_ss,k_ss,policyCapital,policyLabor,valueFunction0,tol,options);
end

    valueFunction0 = valueFunctionEGM;
    
end
toc

%% 7. Euler errors analysis

marginalReturnCapital=zeros(1,sizeShocks);
marginalUtility=zeros(1,sizeShocks);
value=zeros(1,sizeShocks);
utility=zeros(1,sizeShocks);

eulerRHS = zeros(sizeCapital,sizeShocks);

for indexShocks = 1:sizeShocks        
    for indexCapital = 1:sizeCapital
               
         capitalTomorrow = policyCapital(indexCapital,indexShocks);
         
         [~, iCapitalTomorrow] = min(abs(gridCapital-capitalTomorrow));   % start on the grid from the value closest to k1

         for iShocksTomorrow = 1:sizeShocks
             z = vShocks(1,iShocksTomorrow);
             aalpha = vShocks(2,iShocksTomorrow);             
            laborTomorrow =  policyLabor(iCapitalTomorrow,iShocksTomorrow);
            consumptionTomorrow = policyConsumption(iCapitalTomorrow,iShocksTomorrow);
            valueTomorrow = valueFunction(iCapitalTomorrow,iShocksTomorrow);

            utility(iShocksTomorrow) = (log(consumptionTomorrow)...
                -eta*laborTomorrow.^2/2).^rrho;
            
            marginalUtility(iShocksTomorrow) = rrho*...
                utility(iShocksTomorrow)^((rrho-1)/rrho)...
             /consumptionTomorrow;
         
            marginalReturnCapital(iShocksTomorrow) = 1 - delta + (aalpha*z*...
            capitalTomorrow^(aalpha-1)*...
            laborTomorrow^(1-aalpha));
        
            value(iShocksTomorrow) = valueTomorrow;
         end
         
         eulerRHS(indexCapital,indexShocks) = beta*(value.^ppsi*mTransition(indexShocks,:)').^((rrho-ppsi)/ppsi).* ...
                    (marginalUtility.*marginalReturnCapital.*value.^(ppsi-rrho))*mTransition(indexShocks,:)';
         
    end
   
end

mUtilityToday = (log(policyConsumption)-eta*policyLabor.^2/2).^rrho;
mMarginalUtilityToday = rrho*mUtilityToday.^((rrho-1)/rrho)./policyConsumption;

EulerLHS =  mMarginalUtilityToday;               

mEulerErrorsBasic = abs(1 - eulerRHS./EulerLHS);
        
mLogEulerErrorsBasic=log10(mEulerErrorsBasic);

avgLogEulerErrorsBasic=mean(mLogEulerErrorsBasic,2);

totAvgLogEulerErrorsBasic=mean(avgLogEulerErrorsBasic);

maxEulerErrorsBasic=max(max(mLogEulerErrorsBasic));

maxErrorsAcrossShocksBasic=max(mLogEulerErrorsBasic,[],2);

%% 8. Plot value and policy functions

% All together
figure
subplot(2,1,1)
plot(gridCapital,valueFunction)
title('Value function Endogenous Grid')
xlabel('k')
ylabel('V')
xlim([minCapital-10, maxCapital+10])
subplot(2,1,2)
h1 = plot(gridCapital,policyCapital);
hold on
h2 = plot(gridCapital,gridCapital);
title('Capital Decision Rule Endogenous Grid')
xlabel('k')
ylabel('k prime')
legend(h2,{'45° line'},'location','northwest')
xlim([minCapital-10, maxCapital+10])

figure
subplot(2,1,1)
plot(gridCapital,policyLabor)
title('Labor Decision Rule')
xlabel('k')
ylabel('l')
xlim([minCapital-10, maxCapital+10])
subplot(2,1,2)
plot(gridCapital,policyConsumption);
title('Consumption Decision Rule')
xlabel('k')
ylabel('c')
xlim([minCapital-10, maxCapital+10])

%% 9. Plot Euler errors

%figure
%plot(gridCapital,mLogEulerErrorsBasic)

%Plotting the location of the mean and worst Euler Errors across the grid.
figure
plot(gridCapital,avgLogEulerErrorsBasic)
hold on
plot(gridCapital,maxErrorsAcrossShocksBasic)
xla=xlabel('Capital');
yla=ylabel('Log10(Euler Errors)');
tit=title('Location of the Errors on the Grid');
legend('Mean Errors','Max Errors')