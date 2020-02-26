
% Homework 2
% ECON 714 - Prof. Jesus Fernandez-Villaverde
% Sara Casella
% November 2018

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

%% Value function iteration with Accelerator

% Initial guess for value function: steady state value
mValueFunction = zeros(nGridCapital,nGridShocks)+v_ss; 

% Preallocation for required matrices
mValueFunctionNew     = zeros(nGridCapital,nGridShocks);
mPolicyCapital        = zeros(nGridCapital,nGridShocks);
mPolicyLabor          = zeros(nGridCapital,nGridShocks);
mPolicyConsumption    = zeros(nGridCapital,nGridShocks);

expectedValueFunction = zeros(nGridCapital,nGridShocks);

mPeriodUtilityAccelerator = zeros(nGridCapital,nGridShocks);

% Tolerance and max iterations
tolerance     = 1e-06;
maxIterations = 1000;
error         = 1;
iteration     = 1;


tic
while error > tolerance && iteration <= maxIterations
    
    % Continuation Value with Epstein Zin Preferences
     expectedValueFunction = mValueFunction.^(1-gamma)*mTransition';
    
    for iShocks = 1:nGridShocks
        
         productivity = vShocks(1,iShocks);
         alpha        = vShocks(2,iShocks);       
        
        % We start from previous choice (monotonicity of policy function)
        gridCapitalNextPeriod = 1;
        
        for iCapital = 1:nGridCapital
            
            capitalToday = vGridCapital(iCapital);
                        
            
            if (rem((iteration),10)==0 || iteration ==1)
            
               interk = vGridCapital;
               interv = expectedValueFunction(:,iShocks)';

             capitalChoice = fminbnd(@(kprime) (-getValue(productivity,...
                alpha,capitalToday,kprime,parameters,functions,...
                laborMin(iShocks),laborMax(iShocks),...
                interk,interv)),minCapital,maxCapital);
             [value, laborChoice, consumptionChoice] = getValue(productivity,...
                alpha,capitalToday,capitalChoice,...
                parameters,functions,laborMin(iShocks),laborMax(iShocks),...
                interk,interv);


                mPeriodUtilityAccelerator(iCapital,iShocks)  = getPeriodUtility(consumptionChoice,laborChoice);
                
                mPolicyCapital(iCapital,iShocks)     = capitalChoice;
                mPolicyLabor(iCapital,iShocks)       = laborChoice;
                mPolicyConsumption(iCapital,iShocks) = consumptionChoice;
            
            end
             
              mValueFunctionNew(iCapital,iShocks) = ((1-beta)*mPeriodUtilityAccelerator(iCapital,iShocks)^theta+...
                beta*interp1(vGridCapital,expectedValueFunction(:,iShocks), mPolicyCapital(iCapital,iShocks))^(theta/(1-gamma)))^(1/theta);           

        end
    end
    
    error = max(abs(mValueFunctionNew(:)-mValueFunction(:)));
    mValueFunction = mValueFunctionNew;
    
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, error); 
    end
    
    iteration = iteration+1;
 

end
toc;
fprintf(' Total iterations     = %d\n', iteration-1)

elapsed_value_accelerator = toc;
