% Homework 2
% ECON 714 - Prof. Jesus Fernandez-Villaverde
% Sara Casella
% December 2018

% MULTIGRID SCHEME

clear
close all

%% 1. Calibration

% Parameters

theta = 0.5;
gamma = 10;
beta  = 0.99;
delta = 0.1;

nGridCapital = 100;

% save parameters in an array
parameters.theta = theta;
parameters.gamma = gamma;
parameters.beta  = beta;
parameters.delta = delta;

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
                        
% Collapse the two shocks into one made of all the possible combinations
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

%% 2. Deterministic Steady State

z_ss = 1;        % productivity
alpha_ss = 0.3;  % share of capital

l_ss = 100;       % labor normalized to 10
k_ss = (beta*alpha_ss*z_ss*l_ss^(1-alpha_ss)/(1+beta*(delta-1)))^(1/(1-alpha_ss));
c_ss = -delta*k_ss + z_ss*k_ss^alpha_ss*l_ss^(1-alpha_ss);
y_ss = z_ss*k_ss^alpha_ss*l_ss^(1-alpha_ss);

% Normalization in the utility function
eta = (1-alpha_ss)*z_ss*k_ss^alpha_ss*l_ss^(-1-alpha_ss)/c_ss;

% Steady state value. Note we scale up utility by (1-beta)
v_ss = (log(c_ss)-eta*l_ss^2/2);

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

% % Calculate labor, consumption and utility in advance
% mLabor       = zeros(nGridCapital,nGridCapital,nGridShocks);
% mConsumption = zeros(nGridCapital,nGridCapital,nGridShocks);
% mPeriodUtility     = zeros(nGridCapital,nGridCapital,nGridShocks);
% 
% tic
% for iShocks = 1:nGridShocks
%         
%         productivity = vShocks(1,iShocks);
%         alpha        = vShocks(2,iShocks);
%         
%         for iCapital = 1:nGridCapital
%             
%             capitalToday = vGridCapital(iCapital);
%             
%             for iCapitalNextPeriod = 1:nGridCapital
%                 
%                 capitalNextPeriod = vGridCapital(iCapitalNextPeriod);
%                 
%                 labor  = getLabor(productivity,capitalToday,capitalNextPeriod,...
%                     alpha,laborMin(iShocks),laborMax(iShocks));
%                 consumption  = getConsumption(productivity,capitalToday,labor,alpha); 
%                 periodUtility = getPeriodUtility(consumption,labor);
%                 
%                 if periodUtility < 0
%                     periodUtility = 0;
%                 end
%                 
%                 mLabor(iCapital,iCapitalNextPeriod,iShocks) = labor;
%                 mConsumption(iCapital,iCapitalNextPeriod,iShocks) = consumption;              
%                 mPeriodUtility(iCapital,iCapitalNextPeriod,iShocks) = periodUtility;
%                 
% 
%             end
%         end
%             if (mod(iShocks,3)==0)
%                 fprintf('Progress calculating period utility = %d%% \n', iShocks/nGridShocks*100); 
%             end
% end
% toc
% elapsed_value_period_utility = toc/60;
% 
% %% 5. Good initial guess 
% % Value function iteration with a fixed grid using monotonicity
% % NO INTERPOLATION (very fast)
% 
% % Initial guess for value function: steady state value
% mValueFunctionGuess = zeros(nGridCapital,nGridShocks)+v_ss; 
% 
% % Preallocation for required matrices
% mValueFunctionNew     = zeros(nGridCapital,nGridShocks);
% 
% % Tolerance and max iterations
% tolerance     = 1e-06;
% maxIterations = 1000;
% error         = 1;
% iteration     = 1;
% 
% tic
% while error > tolerance && iteration <= maxIterations
%     
%     % Continuation Value with Epstein Zin Preferences
%      expectedValueFunction = mValueFunctionGuess.^(1-gamma)*mTransition';
%     
%     for iShocks = 1:nGridShocks
%         
%         % We start from previous choice (monotonicity of policy function)
%         gridCapitalNextPeriod = 1;
%         
%         for iCapital = 1:nGridCapital
%                                     
%             valueHighSoFar = -1000.0;
%             
%             for iCapitalNextPeriod = gridCapitalNextPeriod:nGridCapital
%                 
%                 periodUtility = mPeriodUtility(iCapital,iCapitalNextPeriod,iShocks);    
% 
%                 valueProvisional = ((1-beta)*periodUtility^theta+...
%                      beta*expectedValueFunction(iCapitalNextPeriod,iShocks)^(theta/(1-gamma)))^(1/theta);
% 
%                 if (valueProvisional>valueHighSoFar)
%                     valueHighSoFar = valueProvisional;
%                     gridCapitalNextPeriod = iCapitalNextPeriod;
%                 else
%                     break; % We break when we have achieved the max
%                 end    
% 
%             end
% 
%             mValueFunctionNew(iCapital,iShocks)  = valueHighSoFar;
%         
%         end
%     end
%     
%     error = max(abs(mValueFunctionNew(:)-mValueFunctionGuess(:)));
%     mValueFunctionGuess = mValueFunctionNew;
%     
%     if (mod(iteration,10)==0 || iteration ==1)
%         fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, error); 
%     end
%     
%     iteration = iteration+1;
%  
% 
% end
% toc;
% fprintf(' Total iterations = %d\n', iteration-1)
% elapsed_value_fixed_grid_monotonicity = toc;


%% 6. First round value function guess
% WITH INTERPOLATION 

%Initial guess for value function: previously found guess
%mValueFunctionMultiGrid1 = mValueFunctionGuess; 
mValueFunctionMultiGrid1 = zeros(nGridCapital,nGridShocks)+v_ss; 

% Preallocation for required matrices
mValueFunctionNew     = zeros(nGridCapital,nGridShocks);
mPolicyCapital        = zeros(nGridCapital,nGridShocks);
expectedValueFunction = zeros(nGridCapital,nGridShocks);

% Tolerance and max iterations
tolerance     = 1e-06;
maxIterations = 1000;
error         = 1;
iteration     = 1;


tic
while error > tolerance && iteration <= maxIterations
    
    % Continuation Value with Epstein Zin Preferences
     expectedValueFunction = mValueFunctionMultiGrid1.^(1-gamma)*mTransition';

    for iShocks = 1:nGridShocks
        
        productivity = vShocks(1,iShocks);
        alpha        = vShocks(2,iShocks);
                
         for iCapital = 1:nGridCapital
             
             capitalToday = vGridCapital(iCapital);

                    interk = vGridCapital;
                    interv = expectedValueFunction(:,iShocks)';

                    [capitalChoice, value] = fminbnd(@(kprime) (-getValue(productivity,...
                        alpha,capitalToday,kprime,parameters,functions,laborMin(iShocks)...
                        ,laborMax(iShocks),...
                        interk,interv)),minCapital,maxCapital);
                
                mValueFunctionNew(iCapital,iShocks) = -value;
                mPolicyCapital(iCapital,iShocks) = capitalChoice;
        
         end   
        
    end
    
    error = max(abs(mValueFunctionNew(:)-mValueFunctionMultiGrid1(:)));
    mValueFunctionMultiGrid1 = mValueFunctionNew;
    
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, error); 
    end
    
    iteration = iteration+1;

end
toc;
fprintf(' Total iterations First Round = %d\n', iteration-1)
elapsed_value_multi_grid1 = toc/60;

%% 7a. Second round grid generation

% Generate grid for capital
nGridCapital = 500;
vGridCapital2 = linspace(minCapital,maxCapital,nGridCapital);

%% 7b. Second round value function guess

% Use previously calculated value function

mValueFunctionMultiGrid2 = zeros(nGridCapital,nGridShocks);

for i = 1:nGridShocks
    mValueFunctionMultiGrid2(:,i) = interpn(vGridCapital,mValueFunctionMultiGrid1(:,i),...
        vGridCapital2);
end

%% 7c. Second round Value function iteration

% Preallocation for required matrices
mValueFunctionNew     = zeros(nGridCapital,nGridShocks);
mPolicyCapital        = zeros(nGridCapital,nGridShocks);

% Tolerance and max iterations
error         = 1;
iteration     = 1;

tic
while error > tolerance && iteration <= maxIterations
    
    % Continuation Value with Epstein Zin Preferences
     expectedValueFunction = mValueFunctionMultiGrid2.^(1-gamma)*mTransition';

    for iShocks = 1:nGridShocks
        
        productivity = vShocks(1,iShocks);
        alpha        = vShocks(2,iShocks);
                
         for iCapital = 1:nGridCapital
             
             capitalToday = vGridCapital2(iCapital);

                    interk = vGridCapital2;
                    interv = expectedValueFunction(:,iShocks)';

                    [capitalChoice, value] = fminbnd(@(kprime) (-getValue(productivity,...
                        alpha,capitalToday,kprime,parameters,functions,...
                        laborMin(iShocks),laborMax(iShocks),...
                        interk,interv)),minCapital,maxCapital);
                
                mValueFunctionNew(iCapital,iShocks) = -value;
                mPolicyCapital(iCapital,iShocks) = capitalChoice;
        
         end   
        
    end
    
    error = max(abs(mValueFunctionNew(:)-mValueFunctionMultiGrid2(:)));
    mValueFunctionMultiGrid2 = mValueFunctionNew;
    
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, error); 
    end
    
    iteration = iteration+1;
end
toc;
fprintf(' Total iterations Second Round = %d\n', iteration-1)
elapsed_value_multi_grid2 = toc/60;

%% 8a. Third round grid generation

% Generate grid for capital
nGridCapital = 5000;
vGridCapital3 = linspace(minCapital,maxCapital,nGridCapital);

%% 8b. Third round value function guess

% Use previously calculated value function

mValueFunctionMultiGrid3 = zeros(nGridCapital,nGridShocks);

for i = 1:nGridShocks
    mValueFunctionMultiGrid3(:,i) = interpn(vGridCapital2,mValueFunctionMultiGrid2(:,i),...
        vGridCapital3);
end

%% 8c. Third round Value function iteration

% Preallocation for required matrices
mValueFunctionNew     = zeros(nGridCapital,nGridShocks);
mPolicyCapital        = zeros(nGridCapital,nGridShocks);
mPolicyLabor          = zeros(nGridCapital,nGridShocks);
mPolicyConsumption    = zeros(nGridCapital,nGridShocks);

% Tolerance and max iterations
error         = 1;
iteration     = 1;

tic
while error > tolerance && iteration <= maxIterations
    
    % Continuation Value with Epstein Zin Preferences
     expectedValueFunction = mValueFunctionMultiGrid3.^(1-gamma)*mTransition';

    for iShocks = 1:nGridShocks
        
        productivity = vShocks(1,iShocks);
        alpha        = vShocks(2,iShocks);
                
         for iCapital = 1:nGridCapital
             
             capitalToday = vGridCapital3(iCapital);

                    interk = vGridCapital3;
                    interv = expectedValueFunction(:,iShocks)';

                    capitalChoice = fminbnd(@(kprime) (-getValue(productivity,...
                        alpha,capitalToday,kprime,parameters,functions,...
                        laborMin(iShocks),laborMax(iShocks),...
                        interk,interv)),minCapital,maxCapital);
                    [value, laborChoice, consumptionChoice] = ...
                        getValue(productivity,alpha,capitalToday,capitalChoice,...
                        parameters,functions,laborMin(iShocks),laborMax(iShocks),...
                        interk,interv);
                
                mValueFunctionNew(iCapital,iShocks) = value;
                mPolicyCapital(iCapital,iShocks) = capitalChoice;
                mPolicyLabor(iCapital,iShocks) = laborChoice;
                mPolicyConsumption(iCapital,iShocks) = consumptionChoice;
         end   
        
    end
    
    error = max(abs(mValueFunctionNew(:)-mValueFunctionMultiGrid3(:)));
    mValueFunctionMultiGrid3 = mValueFunctionNew;
    
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, error); 
    end
    
    iteration = iteration+1;

end
toc;
fprintf(' Total iterations Third Round = %d\n', iteration-1)
elapsed_value_multi_grid3 = toc/60;

%% 9. Euler errors analysis

mMarginalReturnCapital=zeros(1,nGridShocks);
mMarginalUtility=zeros(1,nGridShocks);
mValue=zeros(1,nGridShocks);
mUtility=zeros(1,nGridShocks);

EulerRHS = zeros(nGridCapital,nGridShocks);

for iShocks = 1:nGridShocks        
    for iCapital = 1:nGridCapital
               
         capitalTomorrow = mPolicyCapital(iCapital,iShocks);
         
         [~, iCapitalTomorrow] = min(abs(vGridCapital3-capitalTomorrow));   % start on the grid from the value closest to k1

         for iShocksTomorrow = 1:nGridShocks
             z = vShocks(1,iShocksTomorrow);
             alpha = vShocks(2,iShocksTomorrow);             
            laborTomorrow =  mPolicyLabor(iCapitalTomorrow,iShocksTomorrow);
            consumptionTomorrow = mPolicyConsumption(iCapitalTomorrow,iShocksTomorrow);
            valueTomorrow = mValueFunctionMultiGrid3(iCapitalTomorrow,iShocksTomorrow);

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

mEulerErrorsMulti = abs(1 - EulerRHS./EulerLHS);
        
mLogEulerErrorsMulti=log10(mEulerErrorsMulti);

%figure
%plot(vGridCapital,mLogEulerErrorsBasic)

avgLogEulerErrorsMulti=mean(mLogEulerErrorsMulti,2);

totAvgLogEulerErrorsMulti=mean(avgLogEulerErrorsMulti);

maxEulerErrorsMulti=max(max(mLogEulerErrorsMulti));

maxErrorsAcrossShocksMulti=max(mLogEulerErrorsMulti,[],2);

save multigrid2.mat

%% Plot

% All together
figure
subplot(2,1,1)
plot(vGridCapital3,mValueFunctionMultiGrid3)
title('Value Function Multi Grid')
xlabel('k')
ylabel('V')
xlim([minCapital-10, maxCapital+10])
subplot(2,1,2)
h1 = plot(vGridCapital3,mPolicyCapital);
hold on
h2 = plot(vGridCapital3,vGridCapital3);
title('Capital Decision Rule Multi Grid')
xlabel('k')
ylabel('k prime')
legend(h2,{'45° line'},'location','northwest')
xlim([minCapital-10, maxCapital+10])

% All together
figure
subplot(2,1,1)
plot(vGridCapital3,mPolicyLabor)
title('Labor Decision Rule Multi Grid')
xlabel('k')
ylabel('V')
xlim([minCapital-10, maxCapital+10])
subplot(2,1,2)
plot(vGridCapital3,mPolicyConsumption);
title('Consumption Decision Rule Multi Grid')
xlabel('k')
ylabel('k prime')
xlim([minCapital-10, maxCapital+10])

figure
plot(vGridCapital3,avgLogEulerErrorsMulti)
hold on
plot(vGridCapital3,maxErrorsAcrossShocksMulti)
xla=xlabel('Capital');
yla=ylabel('Log10(Euler Errors)');
tit=title('Location of the Errors on the Grid Multigrid');
legend('Mean Errors','Max Errors')
