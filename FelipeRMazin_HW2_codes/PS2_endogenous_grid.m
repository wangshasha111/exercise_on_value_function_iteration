%% Problem Set 2
% Estimation of the Neoclassical Growth Model with Epstein-Zin
% Preferences, Productivity Shocks and Capital Share Shocks
% (4) Value Function Iteration with Endogenous Grid
%
% Felipe Ruiz Mazin
% December 10, 2018

%% 0. Housekeeping

clear all
close all
clc

global bbeta ddelta eeta ppsi rrho gridCapital gridProd sizeProd gridAlpha sizeAlpha transitionProd transitionAlpha valueFunction0

tic

%%  1. Calibration

bbeta  = 0.99;  % Discount factor
ddelta = 0.1 ;  % Depreciation rate
rrho   = 0.5;   % Elasticity of intertemporal substitution factor: 1/(1-rrho)
ppsi   = -9;    % Risk aversion factor

% Productivity values


% Transition matrices
transitionProd  = [0.9727 0.0273 0 0 0; ...
                   0.0041 0.9806 0.0153 0 0; ...
                   0 0.0082 0.9836 0.0082 0; ...
                   0 0 0.0153 0.9806 0.0041; ...
                   0 0 0 0.0273 0.9727];

transitionAlpha = [0.9 0.07 0.03; ...
                   0.05 0.9 0.05; ...
                   0.03 0.07 0.9];

%% 2. Steady State

laborSteadyState = 5;
y = fsolve(@(x) sstate(x,laborSteadyState), [4 5]);

capitalSteadyState     = y(1);
eeta                   = y(2);
outputSteadyState      = capitalSteadyState^0.3 * laborSteadyState^0.7;
consumptionSteadyState = capitalSteadyState^0.3 * laborSteadyState^0.7 - (1 - (1 - ddelta))*capitalSteadyState;
utilitySteadyState     = log(capitalSteadyState^0.3 * laborSteadyState^0.7 - 0.1*capitalSteadyState) - eeta*laborSteadyState^2/2;

fprintf(' Output = %2.6f, Capital = %2.6f, Consumption = %2.6f\n', outputSteadyState, capitalSteadyState, consumptionSteadyState); 
fprintf('\n')

% Generate the grid for capital
sizeCapital = 50;
capitalMin  = 0.7 * capitalSteadyState;
capitalMax  = 1.3 * capitalSteadyState;
gridCapital = linspace(capitalMin, capitalMax, sizeCapital);
stepCapital = (capitalMax - capitalMin)/sizeCapital;

gridProd = [-0.0673; -0.0336; 0; 0.0336; 0.0673];
sizeProd = length(gridProd);

gridAlpha = [0.25; 0.3; 0.35];
sizeAlpha = length(gridAlpha);               

%% 3. Required matrices and vectors

mktResourcesGrid        = zeros(sizeK,sizeZ,sizeAlpha);
valueFunctionTilde      = zeros(sizeK,sizeZ,sizeAlpha);
valueFunctionTildeDeriv = zeros(sizeK,sizeZ,sizeAlpha);
uMin = 0.1 * uSS^(-1/2);
uMax = 3 * uSS^(-1/2);
valueFunctionTilde0 = repmat(repmat(linspace(uMin,uMax,sizeK)',1,sizeZ),1,1,sizeAlpha);
for indexAlpha = 1:sizeAlpha
    for indexZ = 1:sizeZ
        mktResourcesGrid(:,indexZ,indexAlpha) = exp(zGrid(indexZ)) * kGrid.^(alphaGrid(indexAlpha)) * lSS.^(1 - alphaGrid(indexAlpha)) + 0.9 .* kGrid;
    end
end

valueFunction       = ones(sizeK,sizeZ,sizeAlpha);
consumptionFunction = zeros(sizeK,sizeZ,sizeAlpha);
indexFunction       = ones(sizeK,sizeZ,sizeAlpha);

maxit = 1000;
tol   = 1e-6;

it    = 1;
diff  = 10;

consumption = @(c,v) (1/2*(log(c) - eeta * lSS^2/2)^(-1/2) - c * v);

% Finding value and policy functions numerically
while it <= maxit && diff > tol
    for indexAlpha = 1:sizeAlpha
        for indexZ = 1:sizeZ
            for indexK = 1:sizeK
                it = it + 1;
                
                if indexK == 1 
                    valueFunctionTildeDeriv(indexK,indexZ,indexAlpha) = ...
                        (valueFunctionTilde0(indexK+1,indexZ,indexAlpha)-valueFunctionTilde0(indexK,indexZ,indexAlpha))/stepK;
                elseif indexK == sizeK
                    valueFunctionTildeDeriv(indexK,indexZ,indexAlpha) = ...
                        (valueFunctionTilde0(indexK,indexZ,indexAlpha)-valueFunctionTilde0(indexK-1,indexZ,indexAlpha))/stepK;
                else
                    valueFunctionTildeDeriv(indexK,indexZ,indexAlpha) = ...
                        mean([(valueFunctionTilde0(indexK,indexZ,indexAlpha)-valueFunctionTilde0(indexK-1,indexZ,indexAlpha))/stepK, ...
                        (valueFunctionTilde0(indexK+1,indexZ,indexAlpha)-valueFunctionTilde0(indexK,indexZ,indexAlpha))/stepK]);
                end
                
                
                consumptionFunction(indexK,indexZ,indexAlpha) = fzero(@(c) consumption(c,valueFunctionTildeDeriv(indexK,indexZ,indexAlpha)),50);
                mktResources(indexK,indexZ,indexAlpha)        = consumptionFunction(indexK,indexZ,indexAlpha) + kGrid(indexK);
            
                valueFunction(indexK,indexZ,indexAlpha) = ((log(consumptionFunction(indexK,indexZ,indexAlpha)) - eeta*lSS^2/2)^(1/2) + ...
                    valueFunctionTilde(indexK,indexZ,indexAlpha))^2;
                
                indexYLow  = max(sum(mktResources(indexK,indexZ,indexAlpha) > mktResourcesGrid(:,indexZ,indexAlpha)),1);
                indexYHigh = indexYLow + 1;
                
                valueFunction(indexK,indexZ,indexAlpha) = valueFunction(indexYLow,indexZ,indexAlpha) + (valueFunction(indexYHigh,indexZ,indexAlpha) - ...
                    valueFunction(indexYLow,indexZ,indexAlpha))./(mktResourcesGrid(indexYHigh,indexZ,indexAlpha) - mktResourcesGrid(indexYLow,indexZ,indexAlpha)).*...
                    (mktResources(indexK,indexZ,indexAlpha) - mktResourcesGrid(indexYLow,indexZ,indexAlpha));
                        
            end
        end
    end
    for indexAlpha = 1:sizeAlpha
        for indexZ = 1:sizeZ
            probZ      = repmat(transitionZ(indexZ,:),sizeAlpha,1);
            probAlpha  = repmat(transitionAlpha(indexAlpha,:)',1,sizeZ);

            probMatrix = probZ .* probAlpha;
            for indexK = 1:sizeK
                valueFunctionTilde(indexK,indexZ,indexAlpha) = bbeta * sum(probMatrix .* reshape(valueFunction(indexK,:,:),sizeZ,sizeAlpha).^(-9),'all').^(-0.5/9);
            end
        end
    end
    diff = max(abs(valueFunctionTilde - valueFunctionTilde0),[],'all');
    valueFunctionTilde0 = valueFunctionTilde;
end

toc