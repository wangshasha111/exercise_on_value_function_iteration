% Homework 2
% ECON 714 - Prof. Jesus Fernandez-Villaverde
% Sara Casella
% December 2018

% ENDOGENOUS GRID METHOD

clear
%clc
close all

%% 1. Calibration

% Parameters
theta = 0.5;
gamma = 10;
beta  = 0.99;
delta = 0.1;

% save parameters in an array
% can be passed to functions later
parameters.theta = theta;
parameters.gamma = gamma;
parameters.beta  = beta;
parameters.delta = delta;

% Number of points in capital grid
nGridCapital = 250;

% Shocks processes
vProductivity = exp([-0.0673 -0.0336 0 0.0336 0.0673]);
nGridProductivity = length(vProductivity);

mTransitionProductivity  = [0.9727 0.0273 0      0      0;
                            0.0041 0.9806 0.0153 0      0; 
                            0      0.0082 0.9836 0.0082 0;	
                            0      0      0.0153 0.9806 0.0041;
                            0      0      0      0.0273 0.9727];
                             
vCapitalShare = [0.25 0.3 0.35];
nGridCapitalShare= length(vCapitalShare);

mTransitionCapitalShare =  [0.9  0.07 0.03;
                            0.05 0.9 0.05;
                            0.03 0.07 0.9];
                        
% Collapse the two shocks into one
% 5*3 = 15
vShocks = zeros(2, length(vProductivity)*length(vCapitalShare));
nGridShocks = length(vShocks);

for i=1:length(vProductivity)
    for j=1:length(vCapitalShare)
        
      vShocks(:,j+(i-1)*length(vCapitalShare))=...
          [vProductivity(i), vCapitalShare(j)]';
                                                   
    end 
end

mTransition = kron(mTransitionProductivity,mTransitionCapitalShare); 

options = optimset('Display','off');

%% 2. Deterministic Steady State

z_ss = 1;        % steady state productivity
alpha_ss = 0.3;  % steady state share of capital

l_ss = 100;       % labor normalized to 100
k_ss = (beta*alpha_ss*z_ss*l_ss^(1-alpha_ss)/(1+beta*(delta-1)))^(1/(1-alpha_ss));
c_ss = -delta*k_ss + z_ss*k_ss^alpha_ss*l_ss^(1-alpha_ss);
y_ss = z_ss*k_ss^alpha_ss*l_ss^(1-alpha_ss);

% Normalization in the utility function
eta = (1-alpha_ss)*z_ss*k_ss^alpha_ss*l_ss^(-1-alpha_ss)/c_ss;

% Steady state value. Scale up utility by (1-beta)
v_ss = (log(c_ss)-eta*l_ss^2/2);

% Save eta in parameters array
parameters.eta = eta;

%% 3. Generate grid for capital

minCapital = 0.7*k_ss;
maxCapital = 1.5*k_ss;
vGridCapital = linspace(minCapital,maxCapital,nGridCapital);

%% 4. Calculate period utility in advance

% Functions to get labor consumption and utility
getLabor = @(z,k,kprime,alpha,lowerBound,upperBound) fminbnd(@(l)...
     ((1-alpha)*z*k^alpha*l.^(-(1+alpha))-eta*z*k^alpha*l.^(1-alpha)-eta*(1-delta)*k+eta*kprime).^2,...
     lowerBound, upperBound); % from mkt clearing
getConsumption = @(z,k,l,alpha) ((1-alpha)*z*k^alpha)/(eta*l^(1+alpha)); % from consumption/labor tradeoff
getPeriodUtility = @(c,l) log(c)-eta*l^2/2; % from utility function

% Save functions in an array
functions.getLabor = getLabor;
functions.getConsumption = getConsumption;
functions.getPeriodUtility = getPeriodUtility; 

% Get some good lower and upper bounds for function getLabor
laborMax = zeros(1,nGridShocks);
laborMin = zeros(1,nGridShocks);

for iShocks = 1:nGridShocks
    productivity = vShocks(1,iShocks);
    alpha        = vShocks(2,iShocks);
    laborMin(iShocks) = getLabor(productivity,maxCapital,minCapital,alpha,0.0001,5*l_ss);
    laborMax(iShocks) = getLabor(productivity,minCapital,maxCapital,alpha,0.0001,5*l_ss);
end

% Calculate labor, consumption and utility in advance
mLabor       = zeros(nGridCapital,nGridCapital,nGridShocks);
mConsumption = zeros(nGridCapital,nGridCapital,nGridShocks);
mPeriodUtility     = zeros(nGridCapital,nGridCapital,nGridShocks);

tic
for iShocks = 1:nGridShocks
        
        productivity = vShocks(1,iShocks);
        alpha        = vShocks(2,iShocks);
        
        for iCapital = 1:nGridCapital
            
            capitalToday = vGridCapital(iCapital);
            
            for iCapitalNextPeriod = 1:nGridCapital
                
                capitalNextPeriod = vGridCapital(iCapitalNextPeriod);
                
                labor  = getLabor(productivity,capitalToday,capitalNextPeriod,...
                    alpha,laborMin(iShocks),laborMax(iShocks));
                consumption  = getConsumption(productivity,capitalToday,labor,alpha); 
                periodUtility = getPeriodUtility(consumption,labor);
                
                if periodUtility < 0
                    periodUtility = 0;
                end
                
                mLabor(iCapital,iCapitalNextPeriod,iShocks) = labor;
                mConsumption(iCapital,iCapitalNextPeriod,iShocks) = consumption;              
                mPeriodUtility(iCapital,iCapitalNextPeriod,iShocks) = periodUtility;
                

            end
        end
            if (mod(iShocks,3)==0)
                fprintf('Progress calculating period utility = %d%% \n', iShocks/nGridShocks*100); 
            end
end
toc
elapsed_value_period_utility = toc/60;

%% 5. Do VFI without interpolation to initialize close to the solution
% Value function iteration with a fixed grid using monotonicity
% NO INTERPOLATION (very fast)

% Initial guess for value function: steady state value
mValueFunctionGuess = zeros(nGridCapital,nGridShocks)+v_ss; 

% Preallocation for required matrices
mValueFunctionNew     = zeros(nGridCapital,nGridShocks);

% Tolerance and max iterations
tolerance     = 1e-07;
maxIterations = 1000;
error         = 1;
iteration     = 1;

tic
while error > tolerance && iteration <= maxIterations
    
    % Continuation Value with Epstein Zin Preferences
     expectedValueFunction = mValueFunctionGuess.^(1-gamma)*mTransition';
    
    for iShocks = 1:nGridShocks
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for iCapital = 1:nGridCapital
                                    
            valueHighSoFar = -1000.0;
            
            for iCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital
                
                periodUtility = mPeriodUtility(iCapital,iCapitalNextPeriod,iShocks);    

                valueProvisional = ((1-beta)*periodUtility^theta+...
                     beta*expectedValueFunction(iCapitalNextPeriod,iShocks)^(theta/(1-gamma)))^(1/theta);

                if (valueProvisional>valueHighSoFar)
                    valueHighSoFar = valueProvisional;
                    gridCapitalNextPeriod = iCapitalNextPeriod;
                else
                    break; % We break when we have achieved the max
                end    

            end

            mValueFunctionNew(iCapital,iShocks)  = valueHighSoFar;
        
        end
    end
    
    error = max(abs(mValueFunctionNew(:)-mValueFunctionGuess(:)));
    mValueFunctionGuess = mValueFunctionNew;
    
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration Guess = %d, Sup Diff = %2.8f\n', iteration, error); 
    end
    
    iteration = iteration+1;
 

end
toc;
fprintf(' Total iterations Guess = %d\n', iteration-1)
elapsed_value_fixed_grid_monotonicity = toc;


%% 5. EGM 1
% Implement EGM with labor fixed to steady state value until convergence
% is achieved

%Setting up Vtilde at n=0
mValueFunctionTilde = beta*(mValueFunctionGuess).^(theta);

% Preallocate
mEndogenousValueFunction   = zeros(nGridCapital,nGridShocks);
mEndogenousMarketResources = zeros(nGridCapital,nGridShocks);
mValueFunctionNew          = zeros(nGridCapital,nGridShocks);
mValueFunctionTildeNew     = zeros(nGridCapital,nGridShocks);
mEndogenousConsumption     = zeros(nGridCapital,nGridShocks);

% Setting the Market Resources
mGridOutput = zeros(nGridCapital,nGridShocks);
for iShocks = 1:nGridShocks
        
        productivity = vShocks(1,iShocks);
        alpha        = vShocks(2,iShocks);
        mGridOutput(:,iShocks) = productivity*(vGridCapital'.^alpha*l_ss.^(1-alpha))+(1-delta)*vGridCapital';
end

% Tolerance and max iterations
tolerance     = 1e-06;
maxIterations = 1000;
error         = 1;
iteration     = 1;

tic
while error > tolerance && iteration <= maxIterations
        for iShocks = 1:nGridShocks
             for iCapital = 1:nGridCapital
                 
% Take the derivative of ValueFunctionTilde
% Use average of the slopes of the linearly interpolated value function
        
                if iCapital == nGridCapital
                    
                    derivativeValueFunctionTilde = (mValueFunctionTilde(iCapital,iShocks)-...
                    mValueFunctionTilde(iCapital-1,iShocks))/...
                    (vGridCapital(iCapital)-vGridCapital(iCapital-1));
                
                elseif iCapital == 1
                    
                    derivativeValueFunctionTilde = (mValueFunctionTilde(iCapital+1,iShocks)-...
                    mValueFunctionTilde(iCapital,iShocks))/...
                    (vGridCapital(iCapital+1)-vGridCapital(iCapital));
                                                   
                else
                    
                    derivativeValueFunctionTilde = ((mValueFunctionTilde(iCapital,iShocks)...
                        -mValueFunctionTilde(iCapital-1,iShocks))/...
                       (vGridCapital(iCapital)-vGridCapital(iCapital-1))+...
                       (mValueFunctionTilde(iCapital+1,iShocks)-...
                       mValueFunctionTilde(iCapital,iShocks))/...
                       (vGridCapital(iCapital+1)-vGridCapital(iCapital)))/2;
                end
                
% Compute the optimal level of consumption 

                ConsumptionFunction = @(x) x - theta*(1-beta)*...
                    (log(x)-eta*(l_ss^2/2))^(theta-1)*...
                    derivativeValueFunctionTilde^(-1);
                Consumption = fzero(ConsumptionFunction, [20 1000]); 
                
% Compute the value of the endogenously determined market resources
                mEndogenousMarketResources(iCapital,iShocks)=Consumption+vGridCapital(iCapital);

% Update the value function
                mEndogenousValueFunction(iCapital,iShocks) =...
                 ((1-beta)*(log(Consumption) - eta*(l_ss^2)/2)^(1-theta) ...
                 + mValueFunctionTilde(iCapital,iShocks))^(1/(1-theta));
             
                   if isreal(mEndogenousMarketResources(iCapital,iShocks)) == 0|| ...
                           isreal(mEndogenousValueFunction(iCapital,iShocks)) ==0
                        break;
                   else
%Interpolating so that we define the value function on the
% grid for capital established before.
                
                    for iOutput = 1:nGridCapital
                        outputLow = max(sum(mGridOutput(iOutput,iShocks)>mEndogenousMarketResources(:,iShocks)),1);
                        if outputLow == nGridCapital
                            outputLow = nGridCapital-1;
                        end
                        outputHigh = outputLow + 1;
                        mValueFunctionNew(iOutput,iShocks) = mEndogenousValueFunction(outputLow,iShocks)...
                            +(mGridOutput(iOutput,iShocks)-mEndogenousMarketResources(outputLow,iShocks))...
                            *(mEndogenousValueFunction(outputHigh,iShocks)-mEndogenousValueFunction(outputLow,iShocks))...
                            /(mEndogenousMarketResources(outputHigh,iShocks)-mEndogenousMarketResources(outputLow,iShocks));
                    end           
                   end
             end % end capital
        end % end shocks
        
% Calculate Vtilde in t+1
        mValueFunctionTildeNew = beta*(mValueFunctionNew.^(1-gamma)*mTransition').^(theta/(1-gamma));
% Criteria for convergence
        error = max(abs(mValueFunctionTildeNew(:)-mValueFunctionTilde(:)));
        mValueFunctionTilde = mValueFunctionTildeNew;
    
        if (mod(iteration,10)==0 || iteration ==1)
            fprintf(' Iteration EGM 1 = %d, Sup Diff = %2.8f\n', iteration, error); 
        end
        iteration = iteration+1;
end
fprintf(' Total iterations EGM 1 = %d\n', iteration-1)

%% Retrieve K_endogenous for all values of the shocks
mEndogenousCapitalGrid = zeros(nGridCapital,nGridShocks);

for iShocks = 1:nGridShocks
        productivity = vShocks(1,iShocks);
        alpha        = vShocks(2,iShocks);
        
          CapInitialGuess = k_ss;
        for iCapital = 1:nGridCapital
            
            KendEGM1 = @(x) mEndogenousMarketResources(iCapital,iShocks)...
            -productivity*x^alpha*...
            l_ss^(1-alpha)-(1-delta)*x;
            
            mEndogenousCapitalGrid(iCapital,iShocks) = fsolve(KendEGM1,k_ss,options);
          
        end
end

%% Interpolated value function on the capital grid to do VFI

mValueFunctionInterpolatedVFI1 = zeros(nGridCapital,nGridShocks);

for iShocks = 1:nGridShocks
        
    for iCapital = 1:nGridCapital

        capitalLow = max(sum(vGridCapital(iCapital)>mEndogenousCapitalGrid(:,iShocks)),1);

        if capitalLow == nGridCapital
           capitalLow = nGridCapital-1;
        end

        capitalHigh = capitalLow + 1;

        mValueFunctionInterpolatedVFI1(iCapital,iShocks) = mValueFunctionNew(capitalLow,iShocks)...
            +(vGridCapital(iCapital)-mEndogenousCapitalGrid(capitalLow,iShocks))...
            *(mValueFunctionNew(capitalHigh,iShocks)-mValueFunctionNew(capitalLow,iShocks))...
            /(mEndogenousCapitalGrid(capitalHigh,iShocks)-mEndogenousCapitalGrid(capitalLow,iShocks));

    end
end

%% 6. Generalized EGM
mValueFunctionGuess = mValueFunctionInterpolatedVFI1;

tolerance     = 1e-06;
maxIterations = 3;
error = 1;
tic
while error > tolerance
    %% 6a. VFI1 algorithm to recover the policy functions.

    [mValueFunction, mPolicyCapital, mPolicyLabor, mPolicyConsumption, error] = ...
        VFIinterpolation(mValueFunctionGuess,vGridCapital,...
        vShocks,mTransition,minCapital,maxCapital,laborMin,...
    laborMax,functions,parameters,tolerance,maxIterations);

    fprintf('\n Error = %2.8f\n', error)
   
    
    %% 6b. EGM algorithm to get next value function guess
if error > 10^(-5)
    mValueFunctionEGM = EGM(parameters,vGridCapital,vShocks,...
        mTransition,c_ss,k_ss,mPolicyCapital,mPolicyLabor,mValueFunctionGuess,10*tolerance,options);
else
        mValueFunctionEGM = EGM(parameters,vGridCapital,vShocks,...
        mTransition,c_ss,k_ss,mPolicyCapital,mPolicyLabor,mValueFunctionGuess,tolerance,options);
end

    mValueFunctionGuess = mValueFunctionEGM;
    
end
toc

%% 7. Euler errors analysis

mMarginalReturnCapital=zeros(1,nGridShocks);
mMarginalUtility=zeros(1,nGridShocks);
mValue=zeros(1,nGridShocks);
mUtility=zeros(1,nGridShocks);

EulerRHS = zeros(nGridCapital,nGridShocks);

for iShocks = 1:nGridShocks        
    for iCapital = 1:nGridCapital
               
         capitalTomorrow = mPolicyCapital(iCapital,iShocks);
         
         [~, iCapitalTomorrow] = min(abs(vGridCapital-capitalTomorrow));   % start on the grid from the value closest to k1

         for iShocksTomorrow = 1:nGridShocks
             z = vShocks(1,iShocksTomorrow);
             alpha = vShocks(2,iShocksTomorrow);             
            laborTomorrow =  mPolicyLabor(iCapitalTomorrow,iShocksTomorrow);
            consumptionTomorrow = mPolicyConsumption(iCapitalTomorrow,iShocksTomorrow);
            valueTomorrow = mValueFunction(iCapitalTomorrow,iShocksTomorrow);

            mUtility(iShocksTomorrow) = (log(consumptionTomorrow)...
                -eta*laborTomorrow.^2/2).^theta;
            
            mMarginalUtility(iShocksTomorrow) = theta*...
                mUtility(iShocksTomorrow)^((theta-1)/theta)...
             /consumptionTomorrow;
         
            mMarginalReturnCapital(iShocksTomorrow) = 1 - delta + (alpha*z*...
            capitalTomorrow^(alpha-1)*...
            laborTomorrow^(1-alpha));
        
            mValue(iShocksTomorrow) = valueTomorrow;
         end
         
         EulerRHS(iCapital,iShocks) = beta*(mValue.^(1-gamma)*mTransition(iShocks,:)').^((theta+gamma-1)/(1-gamma)).* ...
                    (mMarginalUtility.*mMarginalReturnCapital.*mValue.^(1-gamma-theta))*mTransition(iShocks,:)';
         
    end
   
end

mUtilityToday = (log(mPolicyConsumption)-eta*mPolicyLabor.^2/2).^theta;
mMarginalUtilityToday = theta*mUtilityToday.^((theta-1)/theta)./mPolicyConsumption;

EulerLHS =  mMarginalUtilityToday;               

mEulerErrorsBasic = abs(1 - EulerRHS./EulerLHS);
        
mLogEulerErrorsBasic=log10(mEulerErrorsBasic);

avgLogEulerErrorsBasic=mean(mLogEulerErrorsBasic,2);

totAvgLogEulerErrorsBasic=mean(avgLogEulerErrorsBasic);

maxEulerErrorsBasic=max(max(mLogEulerErrorsBasic));

maxErrorsAcrossShocksBasic=max(mLogEulerErrorsBasic,[],2);

%% 8. Plot value and policy functions

% All together
figure
subplot(2,1,1)
plot(vGridCapital,mValueFunction)
title('Value function Endogenous Grid')
xlabel('k')
ylabel('V')
xlim([minCapital-10, maxCapital+10])
subplot(2,1,2)
h1 = plot(vGridCapital,mPolicyCapital);
hold on
h2 = plot(vGridCapital,vGridCapital);
title('Capital Decision Rule Endogenous Grid')
xlabel('k')
ylabel('k prime')
legend(h2,{'45° line'},'location','northwest')
xlim([minCapital-10, maxCapital+10])

figure
subplot(2,1,1)
plot(vGridCapital,mPolicyLabor)
title('Labor Decision Rule')
xlabel('k')
ylabel('l')
xlim([minCapital-10, maxCapital+10])
subplot(2,1,2)
plot(vGridCapital,mPolicyConsumption);
title('Consumption Decision Rule')
xlabel('k')
ylabel('c')
xlim([minCapital-10, maxCapital+10])

%% 9. Plot Euler errors

%figure
%plot(vGridCapital,mLogEulerErrorsBasic)

%Plotting the location of the mean and worst Euler Errors across the grid.
figure
plot(vGridCapital,avgLogEulerErrorsBasic)
hold on
plot(vGridCapital,maxErrorsAcrossShocksBasic)
xla=xlabel('Capital');
yla=ylabel('Log10(Euler Errors)');
tit=title('Location of the Errors on the Grid');
legend('Mean Errors','Max Errors')