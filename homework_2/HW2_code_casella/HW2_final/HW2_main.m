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
nGridCapital = 10;

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
tolerance     = 1e-06;
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
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, error); 
    end
    
    iteration = iteration+1;
 

end
toc;
fprintf(' Total iterations = %d\n', iteration-1)
elapsed_value_fixed_grid_monotonicity = toc;

%% 6. Value function iteration with a fixed grid and interpolation
%INTERPOLATING ON THE WHOLE GRID
%INTERPOLATING ONLY ON EXPECTED CONTINUATION VALUE 

%Initial guess for value function: previously found guess
mValueFunction = mValueFunctionGuess; 
%mValueFunction = zeros(nGridCapital,nGridShocks)+v_ss; 

% Preallocation for required matrices
mValueFunctionNew     = zeros(nGridCapital,nGridShocks);
mPolicyCapital        = zeros(nGridCapital,nGridShocks);
mPolicyLabor          = zeros(nGridCapital,nGridShocks);
mPolicyConsumption    = zeros(nGridCapital,nGridShocks);
expectedValueFunction = zeros(nGridCapital,nGridShocks);

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
                
         for iCapital = 1:nGridCapital
             
             capitalToday = vGridCapital(iCapital);

             interk = vGridCapital;
             interv = expectedValueFunction(:,iShocks)';

             [capitalChoice, value] = fminbnd(@(kprime) (-getValue(productivity,...
                alpha,capitalToday,kprime,parameters,functions,laborMin(iShocks),laborMax(iShocks),...
                interk,interv)),minCapital,maxCapital);
            if error < tolerance*10 % only retrieve policy functions when close to convergence
             [~, laborChoice, consumptionChoice] = getValue(productivity,...
                alpha,capitalToday,capitalChoice,...
                parameters,functions,laborMin(iShocks),laborMax(iShocks),...
                interk,interv);
                mPolicyCapital(iCapital,iShocks)    = capitalChoice;
                mPolicyLabor(iCapital,iShocks)      = laborChoice;
                mPolicyConsumption(iCapital,iShocks) = consumptionChoice;
            end
                
             mValueFunctionNew(iCapital,iShocks) = -value;
        
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
elapsed_value_fixed_grid_interpol = toc;

%save basic.mat

%save basic2.mat

%load basic.mat
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
title('Value function')
xlabel('k')
ylabel('V')
xlim([minCapital-10, maxCapital+10])
subplot(2,1,2)
h1 = plot(vGridCapital,mPolicyCapital);
hold on
h2 = plot(vGridCapital,vGridCapital);
title('Capital Decision Rule')
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

% Shock by shock
mValueFunction = reshape(mValueFunction,nGridCapital,nGridCapitalShare,nGridProductivity);
mPolicyCapital = reshape(mPolicyCapital,nGridCapital,nGridCapitalShare,nGridProductivity);
mPolicyLabor = reshape(mPolicyLabor,nGridCapital,nGridCapitalShare,nGridProductivity);
mPolicyConsumption = reshape(mPolicyConsumption,nGridCapital,nGridCapitalShare,nGridProductivity);

%%
figure
subplot(2,3,1)
plot(vGridCapital,squeeze(mValueFunction(:,2,:)))
legend(strcat('z = ', num2str(round(vProductivity(1),2))),strcat('z = ', ...
    num2str(round(vProductivity(2),2))),...
    strcat('z = ', num2str(round(vProductivity(3),2))),strcat('z = ', ...
    num2str(round(vProductivity(4),2))),...
    strcat('z = ', num2str(round(vProductivity(5),2))),'location','northwest')
xlabel('k')
ylabel('V')
title('Value function for mean \alpha')
xlim([minCapital-10, maxCapital+10])

subplot(2,3,2)
plot(vGridCapital,squeeze(mPolicyCapital(:,2,:)))
%hold on
%plot(vGridCapital,vGridCapital)
legend(strcat('z = ', num2str(round(vProductivity(1),2))),strcat('z = ', num2str(round(vProductivity(2),2))),...
    strcat('z = ', num2str(round(vProductivity(3),2))),strcat('z = ', num2str(round(vProductivity(4),2))),...
    strcat('z = ', num2str(round(vProductivity(5),2))),'45° line','location','northwest')
xlabel('k')
ylabel('kprime')
title('Policy function capital for mean \alpha')
xlim([minCapital-10, maxCapital+10])

subplot(2,3,3)
plot(vGridCapital,squeeze(mPolicyLabor(:,2,:)))
%hold on
%plot(vGridCapital,vGridCapital)
legend(strcat('z = ', num2str(round(vProductivity(1),2))),strcat('z = ', num2str(round(vProductivity(2),2))),...
    strcat('z = ', num2str(round(vProductivity(3),2))),strcat('z = ', num2str(round(vProductivity(4),2))),...
    strcat('z = ', num2str(round(vProductivity(5),2))),'45° line','location','southwest')
xlabel('k')
ylabel('kprime')
title('Policy function labor for mean \alpha')
xlim([minCapital-10, maxCapital+10])

subplot(2,3,4)
plot(vGridCapital,squeeze(mValueFunction(:,:,3)))
legend(strcat('\alpha = ', num2str(round(vCapitalShare(1),2))),strcat('\alpha = ', num2str(round(vCapitalShare(2),2))),...
    strcat('\alpha = ', num2str(round(vCapitalShare(3),2))),'location','northwest')
xlabel('k')
ylabel('V')
title('Value function for mean z')
xlim([minCapital-10, maxCapital+10])


subplot(2,3,5)
plot(vGridCapital,squeeze(mPolicyCapital(:,:,3)))
%hold on
%plot(vGridCapital,vGridCapital)
legend(strcat('\alpha = ', num2str(round(vCapitalShare(1),2))),strcat('\alpha = ', num2str(round(vCapitalShare(2),2))),...
    strcat('\alpha = ', num2str(round(vCapitalShare(3),2))),'45° line','location','northwest')
xlabel('k')
ylabel('kprime')
title('Policy function capital for mean z')
xlim([minCapital-10, maxCapital+10])

subplot(2,3,6)
plot(vGridCapital,squeeze(mPolicyLabor(:,:,3)))
%hold on
%plot(vGridCapital,vGridCapital)
legend(strcat('\alpha = ', num2str(round(vCapitalShare(1),2))),strcat('\alpha = ', num2str(round(vCapitalShare(2),2))),...
    strcat('\alpha = ', num2str(round(vCapitalShare(3),2))),'45° line','location','southwest')
xlabel('k')
ylabel('kprime')
title('Policy function labor for mean z')
xlim([minCapital-10, maxCapital+10])

%% 9. Plot Euler errors

figure
plot(vGridCapital,mLogEulerErrorsBasic)

%Plotting the location of the mean and worst Euler Errors across the grid.
figure
plot(vGridCapital,avgLogEulerErrorsBasic)
hold on
plot(vGridCapital,maxErrorsAcrossShocksBasic)
xla=xlabel('Capital');
yla=ylabel('Log10(Euler Errors)');
tit=title('Location of the Errors on the Grid Basic');
legend('Mean Errors','Max Errors')

%% 10. Calculate impulse response functions
nPeriods = 40;

% Find ergodic steady state when z=1,alpha=0.3
pos=mPolicyCapital(:,2,3)>vGridCapital'; % find first time policy function crosses 45° line
ind = find(pos==1,1,'last');

% Initial position in steady-state
% one time transitory shock to productivity
impulseRespCapitalProd = zeros(1,nPeriods);
impulseRespLaborProd   = zeros(1,nPeriods);
impulseRespConsProd    = zeros(1,nPeriods);

impulseRespCapitalProd(1) = mPolicyCapital(ind,2,3);
impulseRespLaborProd(1)   = mPolicyLabor(ind,2,3);
impulseRespConsProd(1)    = mPolicyConsumption(ind,2,3);
impulseRespCapitalProd(2) = mPolicyCapital(ind,2,5);
impulseRespLaborProd(2)   = mPolicyLabor(ind,2,5);
impulseRespConsProd(2)    = mPolicyConsumption(ind,2,5);

% Interpolate
for t = 3:nPeriods
    
    capitalLowIR = max(sum(impulseRespCapitalProd(t-1)>vGridCapital));
    capitalHighIR = capitalLowIR+1;
    
    impulseRespCapitalProd(t) = interp1([vGridCapital(capitalLowIR),...
        vGridCapital(capitalHighIR)],[mPolicyCapital(capitalLowIR,2,3),...
        mPolicyCapital(capitalHighIR,2,3)],impulseRespCapitalProd(t-1));
    
    impulseRespLaborProd(t)  = interp1([vGridCapital(capitalLowIR),...
        vGridCapital(capitalHighIR)],[mPolicyLabor(capitalLowIR,2,3),...
        mPolicyLabor(capitalHighIR,2,3)],impulseRespCapitalProd(t-1));
    
    impulseRespConsProd(t)  = interp1([vGridCapital(capitalLowIR),...
        vGridCapital(capitalHighIR)],[mPolicyConsumption(capitalLowIR,2,3),...
        mPolicyConsumption(capitalHighIR,2,3)],impulseRespCapitalProd(t-1));    
end

    
% Change to percentage deviations from ss
impulseRespCapitalProdPD = zeros(1,nPeriods);
impulseRespLaborProdPD   = zeros(1,nPeriods);
impulseRespConsProdPD    = zeros(1,nPeriods);
    
for t = 1:nPeriods
    impulseRespCapitalProdPD(t) = (log(impulseRespCapitalProd(t))...
        -log(impulseRespCapitalProd(1)))*100;
    impulseRespLaborProdPD(t) = (log(impulseRespLaborProd(t))...
        -log(impulseRespLaborProd(1)))*100;
    impulseRespConsProdPD(t) = (log(impulseRespConsProd(t))...
        -log(impulseRespConsProd(1)))*100;
end

% Initial position in steady-state
% one time transitory shock to capital share
impulseRespCapitalShare = zeros(1,nPeriods);
impulseRespLaborShare = zeros(1,nPeriods);
impulseRespConsShare = zeros(1,nPeriods);

impulseRespCapitalShare(1)  = mPolicyCapital(ind,2,3);
impulseRespLaborShare(1)    = mPolicyLabor(ind,2,3);
impulseRespConsShare(1)     = mPolicyConsumption(ind,2,3);
impulseRespCapitalShare(2)  = mPolicyCapital(ind,3,3);
impulseRespLaborShare(2)    = mPolicyLabor(ind,3,3);
impulseRespConsShare(2)     = mPolicyConsumption(ind,3,3);

% Interpolate
for t = 3:nPeriods
    
    capitalLowIR = max(sum(impulseRespCapitalShare(t-1)>vGridCapital));
    capitalHighIR = capitalLowIR+1;
    
    impulseRespCapitalShare(t) = interp1([vGridCapital(capitalLowIR),...
        vGridCapital(capitalHighIR)],[mPolicyCapital(capitalLowIR,2,3),...
        mPolicyCapital(capitalHighIR,2,3)],impulseRespCapitalShare(t-1));
    
    impulseRespLaborShare(t)  = interp1([vGridCapital(capitalLowIR),...
        vGridCapital(capitalHighIR)],[mPolicyLabor(capitalLowIR,2,3),...
        mPolicyLabor(capitalHighIR,2,3)],impulseRespCapitalShare(t-1));
    
    impulseRespConsShare(t)  = interp1([vGridCapital(capitalLowIR),...
        vGridCapital(capitalHighIR)],[mPolicyConsumption(capitalLowIR,2,3),...
        mPolicyConsumption(capitalHighIR,2,3)],impulseRespCapitalShare(t-1));    
end

% Change to percentage deviations from ss
impulseRespCapitalSharePD = zeros(1,nPeriods);
impulseRespLaborSharePD   = zeros(1,nPeriods);
impulseRespConsSharePD    = zeros(1,nPeriods);
    
for t = 1:nPeriods
    impulseRespCapitalSharePD(t) = (log(impulseRespCapitalShare(t))...
        -log(impulseRespCapitalShare(1)))*100;
    impulseRespLaborSharePD(t) = (log(impulseRespLaborShare(t))...
        -log(impulseRespLaborShare(1)))*100;
    impulseRespConsSharePD(t) = (log(impulseRespConsShare(t))...
        -log(impulseRespConsShare(1)))*100;
end

%% 11. Plot impulse response functions

zPath = zeros(1,nPeriods)+z_ss;
zPath(2) = vProductivity(5);
alphaPath = zeros(1,nPeriods)+alpha_ss;
alphaPath(2) = vCapitalShare(3);

figure
subplot(2,2,1)
plot(impulseRespCapitalProdPD)
hold on
plot(zeros(1,nPeriods),'k')
title('Capital')
ylabel('% deviation from ss')
subplot(2,2,2)
plot(impulseRespLaborProdPD)
title('Labor')
ylabel('% deviation from ss')
hold on
plot(zeros(1,nPeriods),'k')
subplot(2,2,3)
plot(impulseRespConsProdPD)
title('Consumption')
xlabel('time')
ylabel('% deviation from ss')
hold on
plot(zeros(1,nPeriods),'k')
subplot(2,2,4)
plot(zPath)
hold on
plot(zeros(1,nPeriods)+z_ss,'k')
title('Productivity (level)')
xlabel('time')
ylabel('z ')

figure
subplot(2,2,1)
plot(impulseRespCapitalSharePD)
hold on
plot(zeros(1,nPeriods),'k')
title('Capital')
ylabel('% deviation from ss')
subplot(2,2,2)
plot(impulseRespLaborSharePD)
hold on
plot(zeros(1,nPeriods),'k')
title('Labor')
ylabel('% deviation from ss')
subplot(2,2,3)
plot(impulseRespConsSharePD)
hold on
plot(zeros(1,nPeriods),'k')
title('Consumption')
xlabel('time')
ylabel('% deviation from ss')
subplot(2,2,4)
plot(alphaPath)
hold on
plot(zeros(1,nPeriods)+alpha_ss,'k')
title('Capital Share')
xlabel('time')
ylabel('\alpha (level)')

