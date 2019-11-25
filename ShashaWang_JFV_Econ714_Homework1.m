% Homework 1 Econ714 JFV

% A two good production economy
% Written by Shasha Wang

% This code includes value function iteration with 

% 1. Fixed Grid
% 2. Accelerator
% 3. Multigrid
% 4. Multigrid and Accelerator

close all;
clear;

cd 'E:\Dropbox\fall 19-20\jesus\homework_1';

%% Calibration
bbeta = 0.96;

mmu_1 = 0.5; % utility weight for good 1
mmu_2 = 0.5; % utility weight for good 2

ddelta = 0.1; % depreciation for capital used to produce good 1
% ddelta_2 = 0; % depreciation for capital used to produce good 2, which is irrelevant in this setting since we don't use capital to produce good 2

aalphaK = 0.33;
aalphaL = 0.67;

% aalphaK_2 = 0;
% aalphaL_2 = 1;

%% Shocks
vGrid_a1 = exp([-0.0673; -0.0336; 0; 0.0336; 0.0673]);
Na_1 = length(vGrid_a1);
mProb_a1 = [0.9727 0.0273 0 0 0;...
            0.0041 0.9806 0.0153 0 0;...
            0 0.0082 0.9836 0.0082 0;...
            0 0 0.0153 0.9806 0.0041;...
            0 0 0 0.0273 0.9727];
        
vGrid_a2 = [0.9;1;1.1];
Na_2 = length(vGrid_a2);
mProb_a2 = [0.9 0.1 0;...
            0.05 0.9 0.05;...
            0 0.1 0.9];

% Combine the two shocks into one shock
mProb_a1a2 = kron(mProb_a1,mProb_a2);
[A1,A2] = meshgrid(vGrid_a2,vGrid_a1);
temp=cat(2,A2',A1');
mGrid_a1a2=reshape(temp,[],2);
Na = Na_1 * Na_2;

% 
% mGrid_a1a2 = zeros(2,length(vGrid_a1)*length(vGrid_a2)); 
% for i = 1:length(vGrid_a1)
%     for j = 1:length(vGrid_a2)
% %         mGrid_a1a2(1,(i-1)*length(vGrid_a1)+1:(i-1)*length(vGrid_a1)+length(vGrid_a2))= 
%         (i-1)*length(vGrid_a1)+1:length(vGrid_a2)

global inputs;
inputs.vGrid_a1 = vGrid_a1;
inputs.vGrid_a2 = vGrid_a2;
inputs.mProb_a1 = mProb_a1;
inputs.mProb_a2 = mProb_a2;
inputs.mGrid_a1a2 = mGrid_a1a2;
inputs.mProb_a1a2 = mProb_a1a2;

%% 2. Steady State

% Initially I try to "fsolve" for labor_1 and then solve for the other two
% inputs. Then for robustness I use "fsolve" to work on the system of
% equations of the steady state. Since the results come out significantly
% different, I take the result of the second approach.

% The first approach
% labor_1_SteadyState_function = @(labor_1)labor_1 + 1.14486 * labor_1.^0.67 - 0.9346 * labor_1.^0.165;
% labor_1_y=labor_1_SteadyState_function([0:0.01:6]);
% figure(1);
% plot([0:0.01:6],labor_1_y);
% yline(0); % We can see the nonzero solution lies between 0 and 1, less than 0.5, around 0.3.
% 
% % We could use grid search to find labor_1. But now let's use fsolve.
% opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','iter');
% labor_1_initial = 0.3; 
% display('Start Solving the Function');
% labor_1_SteadyState = fsolve(@(labor_1) labor_1_SteadyState_function(labor_1), labor_1_initial,opts1);
% display('Solution Obtained.');
% 
% % Print results
% fprintf('The Candidate Solution Is Found to Be: %2.4f \n', labor_1_SteadyState);
% fprintf('The Function Value At this Candidate Solution Is: %2.6f \n', ...
%         labor_1_SteadyState_function(labor_1_SteadyState));
% 
% k_SteadyState = 3.53290 * labor_1_SteadyState;
% labor_2_SteadyState = 1.14486 * labor_1_SteadyState^0.67;
% 
% T = table(k_SteadyState,labor_1_SteadyState,labor_2_SteadyState)

% The second approach: Use fsolve to solve the system of equations
    
input_ss_initial=[0.9,0.2,0.5];
opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','iter');

display('Start Solving the Function');
SS = fsolve(@(input_ss) steadyStateFunction(input_ss,bbeta,ddelta,aalphaK,aalphaL,mmu_1,mmu_2), input_ss_initial,opts1);
display('Solution Obtained.');
kSteadyState = SS(1);
labor_1_SteadyState = SS(2);
labor_2_SteadyState = SS(3);
fprintf('The Candidate Solution Is Found to Be: %2.4f \n', SS);
fprintf('The Function Value At this Candidate Solution Is: %2.6f \n', ...
        steadyStateFunction(SS,bbeta,ddelta,aalphaK,aalphaL,mmu_1,mmu_2));
T = table(kSteadyState,labor_1_SteadyState,labor_2_SteadyState)

%% 3. Value function iteration with a Fixed Grid
% Iterate on the Value function implied by the Social Planner’s Problem using linear interpolation
% until the change in the sup norm between to iterations is less than 1e-6.
% Compute the Policy function.
% Describe the responses of the economy to a technology shock.
kSpread = 0.3;
kMax = kSteadyState * (1 + kSpread);
kMin = kSteadyState * (1 - kSpread);
Nk = 250;
vGrid_k = linspace(kMin,kMax,Nk)';
inputs.vGrid_k = vGrid_k;

% We can see from the question that once given a1,a2,k,and chosen k',
% labor_1, labor_2 are both chosen by optimization, and consumption_1,
% consumption_2 are both chosen immediately.

% So we need to write a function laborFunction.m from a1,a2,k,k' to labor_1,labor_2
% And we need to write consumption functions consumptionFunction1.m,
% consumptionFunction2.m

%% Efficiency Matrices to save time
% Because fsolve takes a long time in iteration to solve for labor_1,
% labor_2, I create the matrices for them as well as for consumption_1 and
% consumption_2 to retrieve from later.

% Note we have to put "kPrime" and "k" on the first two dimensions in order not
% to cause any problem of dimensionality during interpolation

mLabor_1Fsolve = zeros(Nk,Nk,Na); % kPrime,k,[a_1,a_2]
mLabor_2Fsolve = zeros(Nk,Nk,Na);
mConsumption_1Fsolve = zeros(Nk,Nk,Na);
mConsumption_2Fsolve = zeros(Nk,Nk,Na);
mCurrentUtilityFsolve = zeros(Nk,Nk,Na);
% mMarginalUtilityTodayFsolve = zeros(Nk,Na,Nk);
% mMarginalUtilityTomorrowFsolve = ;
laborInitial=[labor_1_SteadyState,labor_2_SteadyState];
opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','off');

tic
for ia = 1:Na
    a_1 = mGrid_a1a2(ia,1);
    a_2 = mGrid_a1a2(ia,2);
    for ik = 1:Nk
        
        k = vGrid_k(ik);
        
        for ikPrime = 1:Nk
            kPrime = vGrid_k(ikPrime);
            
            vLaborFsolve = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);

            mLabor_1Fsolve(ikPrime,ik,ia) = vLaborFsolve(1);
            mLabor_2Fsolve(ikPrime,ik,ia) = vLaborFsolve(2);
            mConsumption_1Fsolve(ikPrime,ik,ia) = consumptionFunction1(a_1,k,kPrime,vLaborFsolve(1),aalphaK,aalphaL,ddelta);
            mConsumption_2Fsolve(ikPrime,ik,ia) = consumptionFunction2(a_2,vLaborFsolve(2));
            mCurrentUtilityFsolve(ikPrime,ik,ia) = utilityFunction(mConsumption_1Fsolve(ikPrime,ik,ia),mConsumption_2Fsolve(ikPrime,ik,ia),vLaborFsolve(1),vLaborFsolve(2),mmu_1,mmu_2);
%             mMarginalUtilityTodayFsolve(ikPrime,ia,ik) = mmu_1 * (mConsumption_2Fsolve(ikPrime,ia,ik))^mmu_2 * (mConsumption_1Fsolve(ikPrime,ia,ik))^(mmu_1-1);
            laborInitial=[vLaborFsolve(1),vLaborFsolve(2)];
        end
        
    end
end
toc

inputs.mLabor_1Fsolve = mLabor_1Fsolve;
inputs.mLabor_2Fsolve = mLabor_2Fsolve;
inputs.mConsumption_1Fsolve = mConsumption_1Fsolve;
inputs.mConsumption_2Fsolve = mConsumption_2Fsolve;
inputs.mCurrentUtilityFsolve = mCurrentUtilityFsolve;
% inputs.mMarginalUtilityTodayFsolve = mMarginalUtilityTodayFsolve;
% 
% mLabor_1Fsolve=permute(mLabor_1Fsolve,[3,2,1]);
% mLabor_2Fsolve=permute(mLabor_2Fsolve,[3,2,1]);
% mConsumption_1Fsolve=permute(mConsumption_1Fsolve,[3,2,1]);
% mConsumption_2Fsolve=permute(mConsumption_2Fsolve,[3,2,1]);
% mCurrentUtilityFsolve=permute(mCurrentUtilityFsolve,[3,2,1]);

save('efficiencyMatricesNk250','mLabor_1Fsolve','mLabor_2Fsolve','mConsumption_1Fsolve','mConsumption_2Fsolve','mCurrentUtilityFsolve')
% 历时 2516.447092 秒。

%% Required matrices and vectors

% First compute utility at steady state to be the initial guess
consumption_1_SteadyState = consumptionFunction1(1,kSteadyState,kSteadyState,labor_1_SteadyState,aalphaK,aalphaL,ddelta);
consumption_2_SteadyState = consumptionFunction2(1,labor_2_SteadyState);
utilitySteadyState = utilityFunction(consumption_1_SteadyState,consumption_2_SteadyState,labor_1_SteadyState,labor_2_SteadyState,mmu_1,mmu_2);

mValue0        = utilitySteadyState.*ones(Nk,Na);
inputs.mValue0 = mValue0;
mValue         = zeros(Nk,Na);
mKPolicy        = zeros(Nk,Na);
mLaborPolicy_1 = zeros(Nk,Na);
mLaborPolicy_2 = zeros(Nk,Na); 
mConsumptionPolicy_1   = zeros(Nk,Na);
mConsumptionPolicy_2   = zeros(Nk,Na);

%% Main iteration
maxIter = 10000;
tolerance     = 1e-6;

iteration       = 1;
mDifference = zeros(maxIter,1);
mDifference(iteration) = 100;

options = optimset('Display', 'off');
opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','off');
% laborInitial=[labor_1_SteadyState,labor_2_SteadyState];

tic

% Finding value and policy functions numerically
while iteration <= maxIter && mDifference(iteration) > tolerance
    expectedValue0 = mValue0 * mProb_a1a2'; % row: kPrime, column: a, i.e., expected value if today's shock is a and I choose k as tomorrow's capital stock
%     inputs.expectedValue0 = expectedValue0;

    for ia = 1:Na
%         a = mGrid_a1a2(ia,:);
        a_1 = mGrid_a1a2(ia,1);
        a_2 = mGrid_a1a2(ia,2);
        
        for ik = 1:Nk
            k = vGrid_k(ik);
            
            if ik ==1
                        
                [kPrime, fval] = fminbnd(@(kPrime) ...
                    -valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                    vGrid_k(1),min(1.2*vGrid_k(ik),vGrid_k(end)),options);       
            else
                [kPrime, fval] = fminbnd(@(kPrime) ...
                    -valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                    mKPolicy((ik-1),ia),min(1.2*vGrid_k(ik),vGrid_k(end)),options);
            end
            
            mKPolicy(ik,ia) = kPrime;
            mValue(ik,ia) = -fval;     
            
 %% If we can retrieve other outputs of the objective function of fminbnd, it would be much easier and FASTER to get the labor and consumption policy function
            if mDifference(iteration) < tolerance*10 % Store the policy functions in case it's the last iteration
%                 vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
%                 mLaborPolicy_1(ik,ia) = vLabor(1);
%                 mLaborPolicy_2(ik,ia) = vLabor(2);
%                 mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
%                 mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,vLabor(2));
%                 laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
%                 inputs.laborInitial = laborInitial;
                mLaborPolicy_1(ik,ia) = interp1(vGrid_k,mLabor_1Fsolve(:,ia,ik),kPrime);
                mLaborPolicy_2(ik,ia) = interp1(vGrid_k,mLabor_2Fsolve(:,ia,ik),kPrime);
                mConsumptionPolicy_1(ik,ia) = interp1(vGrid_k,mConsumption_1Fsolve(:,ia,ik),kPrime);
                mConsumptionPolicy_2(ik,ia) = interp1(vGrid_k,mConsumption_2Fsolve(:,ia,ik),kPrime);
            end
        end
    end
    
    iteration = iteration + 1;
    mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
    mValue0         = mValue;
    
    if mod(iteration,10) == 2
        fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 
    end
    
%     inputs.valueFunction0 = mValue0;
%     inputs.laborFunction  = mLaborPolicy_1;
%     policyFunction0       = mKPolicy;
end
toc

fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 

%  Iteration:  1, Sup diff: 0.003807
%  Iteration: 11, Sup diff: 0.000832
%  Iteration: 21, Sup diff: 0.000485
%  Iteration: 31, Sup diff: 0.000300
%  Iteration: 41, Sup diff: 0.000191
%  Iteration: 51, Sup diff: 0.000122
%  Iteration: 61, Sup diff: 0.000079
%  Iteration: 71, Sup diff: 0.000051
%  Iteration: 81, Sup diff: 0.000033
%  Iteration: 91, Sup diff: 0.000022
%  Iteration: 101, Sup diff: 0.000014
%  Iteration: 111, Sup diff: 0.000009
%  Iteration: 121, Sup diff: 0.000006
%  Iteration: 131, Sup diff: 0.000004
%  Iteration: 141, Sup diff: 0.000003
%  Iteration: 151, Sup diff: 0.000002
%  Iteration: 161, Sup diff: 0.000001
% toc
% 历时 1590.565320 秒。

%% figures for Value Function Iteration with a Fixed Grid

figure;
[kk,aa]=meshgrid(vGrid_k, mGrid_a1a2(:,1));
mesh(kk, aa, mValue');

title('Value - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Value','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_value_fixed_grid')

figure;
mesh(kk, aa, mKPolicy');

title('Policy for Next Period Capital - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_kPolicy_fixed_grid')

figure;
mesh(kk, aa, mLaborPolicy_1');

title('Policy for Good 1 Labor - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_1_fixed_grid')

figure;
mesh(kk, aa, mLaborPolicy_2');

title('Policy for Good 2 Labor - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_2_fixed_grid')

figure;
mesh(kk, aa, mConsumptionPolicy_1');

title('Policy for Good 1 Consumption - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_1_fixed_grid')

figure;
mesh(kk, aa, mConsumptionPolicy_2');

title('Policy for Good 2 Consumption - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_2_fixed_grid')

% save ShashaWang_JFV_PS1_5_capital_grid_points

%  Iteration:  1, Sup diff: 0.003807
%  Iteration: 11, Sup diff: 0.000836
%  Iteration: 21, Sup diff: 0.000488
%  Iteration: 31, Sup diff: 0.000302
%  Iteration: 41, Sup diff: 0.000192
%  Iteration: 51, Sup diff: 0.000123
%  Iteration: 61, Sup diff: 0.000080
%  Iteration: 71, Sup diff: 0.000052
%  Iteration: 81, Sup diff: 0.000034
%  Iteration: 91, Sup diff: 0.000022
%  Iteration: 101, Sup diff: 0.000014
%  Iteration: 111, Sup diff: 0.000009
%  Iteration: 121, Sup diff: 0.000006
%  Iteration: 131, Sup diff: 0.000004
%  Iteration: 141, Sup diff: 0.000003
%  Iteration: 151, Sup diff: 0.000002
%  Iteration: 161, Sup diff: 0.000001
% time elapsed 2831.227060 s
% 166 iterations
% Only 5 grid points for capital
% If 250 points for capital, it would take up to 2831*50s, i.e.,39.3 hours,
% or 1.6 days to finish iteration.

% save ShashaWang_JFV_PS1_10_capital_grid_points

%  Iteration:  1, Sup diff: 0.003807
%  Iteration: 11, Sup diff: 0.000835
%  Iteration: 21, Sup diff: 0.000487
%  Iteration: 31, Sup diff: 0.000301
%  Iteration: 41, Sup diff: 0.000191
%  Iteration: 51, Sup diff: 0.000123
%  Iteration: 61, Sup diff: 0.000079
%  Iteration: 71, Sup diff: 0.000051
%  Iteration: 81, Sup diff: 0.000033
%  Iteration: 91, Sup diff: 0.000022
%  Iteration: 101, Sup diff: 0.000014
%  Iteration: 111, Sup diff: 0.000009
%  Iteration: 121, Sup diff: 0.000006
%  Iteration: 131, Sup diff: 0.000004
%  Iteration: 141, Sup diff: 0.000003
%  Iteration: 151, Sup diff: 0.000002
%  Iteration: 161, Sup diff: 0.000001
% time elapsed  6566.472522 s。

% save ShashaWang_JFV_PS1_50_capital_grid_points
%  Iteration:  1, Sup diff: 0.003807
%  Iteration: 11, Sup diff: 0.000832
%  Iteration: 21, Sup diff: 0.000486
%  Iteration: 31, Sup diff: 0.000300
%  Iteration: 41, Sup diff: 0.000191
%  Iteration: 51, Sup diff: 0.000122
%  Iteration: 61, Sup diff: 0.000079
%  Iteration: 71, Sup diff: 0.000051
%  Iteration: 81, Sup diff: 0.000033
%  Iteration: 91, Sup diff: 0.000022
%  Iteration: 101, Sup diff: 0.000014
%  Iteration: 111, Sup diff: 0.000009
%  Iteration: 121, Sup diff: 0.000006
%  Iteration: 131, Sup diff: 0.000004
%  Iteration: 141, Sup diff: 0.000003
%  Iteration: 151, Sup diff: 0.000002
%  Iteration: 161, Sup diff: 0.000001
% time elapsed  25248.327232 s。

% save ShashaWang_JFV_PS1_100_capital_grid_points_fixed_grid

% toc
%  Iteration:  1, Sup diff: 0.004372
%  Iteration: 11, Sup diff: 0.002907
%  Iteration: 21, Sup diff: 0.001932
%  Iteration: 31, Sup diff: 0.001285
%  Iteration: 41, Sup diff: 0.000854
%  Iteration: 51, Sup diff: 0.000568
%  Iteration: 61, Sup diff: 0.000378
%  Iteration: 71, Sup diff: 0.000251
%  Iteration: 81, Sup diff: 0.000167
%  Iteration: 91, Sup diff: 0.000111
%  Iteration: 101, Sup diff: 0.000074
%  Iteration: 111, Sup diff: 0.000049
%  Iteration: 121, Sup diff: 0.000033
%  Iteration: 131, Sup diff: 0.000022
%  Iteration: 141, Sup diff: 0.000014
%  Iteration: 151, Sup diff: 0.000010
%  Iteration: 161, Sup diff: 0.000006
%  Iteration: 171, Sup diff: 0.000004
%  Iteration: 181, Sup diff: 0.000003
%  Iteration: 191, Sup diff: 0.000002
%  Iteration: 201, Sup diff: 0.000001
% time elapsed 2731.808359 s, i.e., 46min。

%% For accuracy test, compute the euler equation error

% % First let's not use linear interpolation for tomorrow's values
% 
% leftHandSideOfEulerEquation = zeros(Nk,Na); % Marginal utility of today
% rightHandSideOfEulerEquation = zeros(Nk,Na); % expected Marginal utility of tomorrow
% % laborInitial=[labor_1_SteadyState,labor_2_SteadyState];
% 
% % Compute marginal utility of today
% for ia = 1:Na
%     a = mGrid_a1a2(ia,:);
%     a_1 = mGrid_a1a2(ia,1);
%     a_2 = mGrid_a1a2(ia,2);        
%     
%     for ik = 1:Nk
% %         k = vGrid_k(ik);
%         kPrime = mKPolicy(ik,ia);
% %         kPrimePrime = depend on aPrime
% %         vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
%         [v,ikPrime]=min(abs(kPrime - vGrid_k));
%         
%         consumption_1 = mConsumptionPolicy_1(ik,ia);
%         consumption_2 = mConsumptionPolicy_2(ik,ia);
% %         labor_1 = mLaborPolicy_1(ik,ia);
% %         labor_2 = mLaborPolicy_2(ik,ia);
% %         laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
% %         leftHandSideOfEulerEquation(ik,ia) = mmu_1 * (a_2 * labor_2)^mmu_2 * ...
% %             (a_1 * k^aalphaK * labor_1^aalphaL + (1-ddelta) * k - kPrime)^(mmu_1-1);
%        
%         marginalUtilityToday = mmu_1 * (consumption_2)^mmu_2 * (consumption_1)^(mmu_1-1);
%         leftHandSideOfEulerEquation(ik,ia) = marginalUtilityToday;
%         
% % Compute expected marginal utility of tomorrow
%         expected = zeros(Na,1);
%         
%         for iaPrime = 1:Na
%             aPrime = mGrid_a1a2(iaPrime,:);
%             aPrime_1 = mGrid_a1a2(iaPrime,1);
%             aPrime_2 = mGrid_a1a2(iaPrime,2); 
%             
% %             kPrimePrime = mKPolicy(ikPrime,iaPrime);
%             consumptionPrime_1 = mConsumptionPolicy_1(ikPrime,iaPrime);
%             consumptionPrime_2 = mConsumptionPolicy_2(ikPrime,iaPrime);
%             laborPrime_1 = mLaborPolicy_1(ikPrime,iaPrime);
%             
%             marginalUtilityTomorrow = mmu_1 * (consumptionPrime_2)^mmu_2 * (consumptionPrime_1)^(mmu_1-1);
%             returnOnCapital = 1 - ddelta + aPrime_1 * laborPrime_1^aalphaL * aalphaK * (kPrime)^(aalphaK-1);
%             
%             unexpected = bbeta * marginalUtilityTomorrow * returnOnCapital;
%             expected(iaPrime) = unexpected * mProb_a1a2(ia,iaPrime);
%         end
%         rightHandSideOfEulerEquation(ik,ia) = sum(expected);
%     end
% end
% 
% errorEulerEquationNoLinearInterpolation = abs(leftHandSideOfEulerEquation - rightHandSideOfEulerEquation);
% errorEulerEquationNoLinearInterpolationDecimalLog = log10( errorEulerEquationNoLinearInterpolation );
% 
% figure(8);
% mesh(kk, aa, errorEulerEquationNoLinearInterpolationDecimalLog');
% 
% title('Euler Equation Error $log_{10}$ No Linear Interpolation','interpreter','latex')
% xlabel('Capital Stock $k$','interpreter','latex')
% ylabel('shocks $z_1$ $z_2$','interpreter','latex')
% zlabel('error','interpreter','latex')
% xlim([min(vGrid_k),max(vGrid_k)])
% ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
% savefig('q3_eulerEquationErrorNoLinearInterpolation_3D')

errorEulerEquationLinearInterpolationFixedGrid = eulerEquationErrorFunction(Nk,vGrid_k,mKPolicy,mLaborPolicy_1,mConsumptionPolicy_1,mConsumptionPolicy_2,Na,mGrid_a1a2,mProb_a1a2,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL);
errorEulerEquationLinearInterpolationFixedGridDecimalLog = log10( errorEulerEquationLinearInterpolationFixedGrid );

figure;
mesh(kk, aa, errorEulerEquationLinearInterpolationFixedGridDecimalLog');

title('Euler Equation Error $log_{10}$ Linear Interpolation - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('error','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_eulerEquationErrorLinearInterpolation_fixed_grid')

% save ShashaWang_JFV_PS1_60_capital_grid_points
%  Iteration:  1, Sup diff: 0.003807
%  Iteration: 11, Sup diff: 0.000832
%  Iteration: 21, Sup diff: 0.000485
%  Iteration: 31, Sup diff: 0.000300
%  Iteration: 41, Sup diff: 0.000191
%  Iteration: 51, Sup diff: 0.000122
%  Iteration: 61, Sup diff: 0.000079
%  Iteration: 71, Sup diff: 0.000051
%  Iteration: 81, Sup diff: 0.000033
%  Iteration: 91, Sup diff: 0.000022
%  Iteration: 101, Sup diff: 0.000014
%  Iteration: 111, Sup diff: 0.000009
%  Iteration: 121, Sup diff: 0.000006
%  Iteration: 131, Sup diff: 0.000004
%  Iteration: 141, Sup diff: 0.000003
%  Iteration: 151, Sup diff: 0.000002
%  Iteration: 161, Sup diff: 0.000001
% toc
% 历时 1590.565320 秒。

% After I take advantage of the monotonicity of the policy function wrt
% capital, the minimazation takes a lot less time
%  Iteration:  1, Sup diff: 0.006262
%  Iteration: 11, Sup diff: 0.000826
%  Iteration: 21, Sup diff: 0.000481
%  Iteration: 31, Sup diff: 0.000298
%  Iteration: 41, Sup diff: 0.000189
%  Iteration: 51, Sup diff: 0.000121
%  Iteration: 61, Sup diff: 0.000078
%  Iteration: 71, Sup diff: 0.000051
%  Iteration: 81, Sup diff: 0.000033
%  Iteration: 91, Sup diff: 0.000022
%  Iteration: 101, Sup diff: 0.000014
%  Iteration: 111, Sup diff: 0.000009
%  Iteration: 121, Sup diff: 0.000006
%  Iteration: 131, Sup diff: 0.000004
%  Iteration: 141, Sup diff: 0.000003
%  Iteration: 151, Sup diff: 0.000002
%  Iteration: 161, Sup diff: 0.000001
% 历时 729.139481 秒。

save ShashaWang_JFV_PS1_250_capital_grid_points
%  Iteration:  1, Sup diff: 0.006056
%  Iteration: 11, Sup diff: 0.000830
%  Iteration: 21, Sup diff: 0.000484
%  Iteration: 31, Sup diff: 0.000300
%  Iteration: 41, Sup diff: 0.000190
%  Iteration: 51, Sup diff: 0.000122
%  Iteration: 61, Sup diff: 0.000079
%  Iteration: 71, Sup diff: 0.000051
%  Iteration: 81, Sup diff: 0.000033
%  Iteration: 91, Sup diff: 0.000022
%  Iteration: 101, Sup diff: 0.000014
%  Iteration: 111, Sup diff: 0.000009
%  Iteration: 121, Sup diff: 0.000006
%  Iteration: 131, Sup diff: 0.000004
%  Iteration: 141, Sup diff: 0.000003
%  Iteration: 151, Sup diff: 0.000002
%  Iteration: 161, Sup diff: 0.000001
% 历时 2263.827634 秒。

%% 6. Value function iteration using Accelerator

%% Required matrices and vectors

mValue0        = utilitySteadyState.*ones(Nk,Na);
inputs.mValue0 = mValue0;
mValue         = zeros(Nk,Na);
mKPolicy        = zeros(Nk,Na);
mLaborPolicy_1 = zeros(Nk,Na);
mLaborPolicy_2 = zeros(Nk,Na); 
mConsumptionPolicy_1   = zeros(Nk,Na);
mConsumptionPolicy_2   = zeros(Nk,Na);

%% Main iteration
maxIter = 10000;
tolerance     = 1e-6;

iteration       = 1;
mDifference = zeros(maxIter,1);
mDifference(iteration) = 100;

options = optimset('Display', 'off');
opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','off');
% laborInitial=[labor_1_SteadyState,labor_2_SteadyState];

tic

% Finding value and policy functions numerically
while iteration <= maxIter  ...
        && ( (mDifference(iteration) > tolerance)  | (mDifference(iteration) <= tolerance && mod(iteration,10)~=2)) % make sure the last iteration does the maximization
    
    expectedValue0 = mValue0 * mProb_a1a2'; % row: kPrime, column: a, i.e., expected value if today's shock is a and I choose k as tomorrow's capital stock
    if (mDifference(iteration) > tolerance)
        for ia = 1:Na
    %         a = mGrid_a1a2(ia,:);
            a_1 = mGrid_a1a2(ia,1);
            a_2 = mGrid_a1a2(ia,2);

            for ik = 1:Nk
                k = vGrid_k(ik);

                if mod(iteration,10) == 1

                    if ik ==1

                        [kPrime, fval] = fminbnd(@(kPrime) ...
                            -valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                            vGrid_k(1),min(1.2*vGrid_k(ik),vGrid_k(end)),options);

                    else
                        [kPrime, fval] = fminbnd(@(kPrime) ...
                        -valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                        mKPolicy((ik-1),ia),min(1.2*vGrid_k(ik),vGrid_k(end)),options);

                    end
                    mKPolicy(ik,ia) = kPrime;
                    mValue(ik,ia) = -fval;    

    %             if mDifference(iteration) < tolerance*10 % Store the policy functions in case it's the last iteration
        %                 vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
        %                 mLaborPolicy_1(ik,ia) = vLabor(1);
        %                 mLaborPolicy_2(ik,ia) = vLabor(2);
        %                 mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
        %                 mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,vLabor(2));
        %                 laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
        %                 inputs.laborInitial = laborInitial;
                    mLaborPolicy_1(ik,ia) = interp1(vGrid_k,mLabor_1Fsolve(:,ia,ik),kPrime);
                    mLaborPolicy_2(ik,ia) = interp1(vGrid_k,mLabor_2Fsolve(:,ia,ik),kPrime);
                    mConsumptionPolicy_1(ik,ia) = interp1(vGrid_k,mConsumption_1Fsolve(:,ia,ik),kPrime);
                    mConsumptionPolicy_2(ik,ia) = interp1(vGrid_k,mConsumption_2Fsolve(:,ia,ik),kPrime);
    %             end
                else
                    currentUtility = interp1(vGrid_k,mCurrentUtilityFsolve(:,ia,ik),mKPolicy(ik,ia));
                    expectedValue = interp1(vGrid_k,expectedValue0(:,ia),mKPolicy(ik,ia));
                    value = (1-bbeta)*currentUtility + bbeta * expectedValue;
                    mValue(ik,ia) = value;
                end
            end
        end
        
        iteration = iteration + 1;
        mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
        mValue0         = mValue;

        if mod(iteration,10) == 2
            fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 
        end
        
    else
        for ia = 1:Na
    %         a = mGrid_a1a2(ia,:);
            a_1 = mGrid_a1a2(ia,1);
            a_2 = mGrid_a1a2(ia,2);

            for ik = 1:Nk
                k = vGrid_k(ik);


                if ik ==1

                    [kPrime, fval] = fminbnd(@(kPrime) ...
                        -valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                        vGrid_k(1),min(1.2*vGrid_k(ik),vGrid_k(end)),options);

                else
                    [kPrime, fval] = fminbnd(@(kPrime) ...
                    -valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                    mKPolicy((ik-1),ia),min(1.2*vGrid_k(ik),vGrid_k(end)),options);

                end
                mKPolicy(ik,ia) = kPrime;
                mValue(ik,ia) = -fval;    

%             if mDifference(iteration) < tolerance*10 % Store the policy functions in case it's the last iteration
    %                 vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
    %                 mLaborPolicy_1(ik,ia) = vLabor(1);
    %                 mLaborPolicy_2(ik,ia) = vLabor(2);
    %                 mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
    %                 mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,vLabor(2));
    %                 laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
    %                 inputs.laborInitial = laborInitial;
                mLaborPolicy_1(ik,ia) = interp1(vGrid_k,mLabor_1Fsolve(:,ia,ik),kPrime);
                mLaborPolicy_2(ik,ia) = interp1(vGrid_k,mLabor_2Fsolve(:,ia,ik),kPrime);
                mConsumptionPolicy_1(ik,ia) = interp1(vGrid_k,mConsumption_1Fsolve(:,ia,ik),kPrime);
                mConsumptionPolicy_2(ik,ia) = interp1(vGrid_k,mConsumption_2Fsolve(:,ia,ik),kPrime);
%             end
            end
        end
        
        iteration = iteration + 1;
        mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
        mValue0         = mValue;

        fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 
        break
    end
    
    

end

toc

fprintf(' Convergence achieved. Total Number of Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 

%% figures for Value Function Iteration with Accelerator

figure;
[kk,aa]=meshgrid(vGrid_k, mGrid_a1a2(:,1));
mesh(kk, aa, mValue');

title('Value - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Value','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_value_accelerator')

figure;
mesh(kk, aa, mKPolicy');

title('Policy for Next Period Capital - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_kPolicy_accelerator')

figure;
mesh(kk, aa, mLaborPolicy_1');

title('Policy for Good 1 Labor - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_1_accelerator')

figure;
mesh(kk, aa, mLaborPolicy_2');

title('Policy for Good 2 Labor - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_2_accelerator')

figure;
mesh(kk, aa, mConsumptionPolicy_1');

title('Policy for Good 1 Consumption - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_1_accelerator')

figure;
mesh(kk, aa, mConsumptionPolicy_2');

title('Policy for Good 2 Consumption - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_2_accelerator')


%% For accuracy test, compute the euler equation error
errorEulerEquationLinearInterpolationAccelerator = eulerEquationErrorFunction(Nk,vGrid_k,mKPolicy,mLaborPolicy_1,mConsumptionPolicy_1,mConsumptionPolicy_2,Na,mGrid_a1a2,mProb_a1a2,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL);
errorEulerEquationLinearInterpolationAcceleratorDecimalLog = log10( errorEulerEquationLinearInterpolationAccelerator );

figure;
mesh(kk, aa, errorEulerEquationLinearInterpolationAcceleratorDecimalLog');

title('Euler Equation Error $log_{10}$ Linear Interpolation - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('error','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_eulerEquationErrorLinearInterpolation_accelerator')
% save ShashaWang_JFV_PS1_60_capital_grid_points_accelerator
%  Iteration:  1, Sup diff: 0.003807
%  Iteration: 11, Sup diff: 0.001170
%  Iteration: 21, Sup diff: 0.000566
%  Iteration: 31, Sup diff: 0.000351
%  Iteration: 41, Sup diff: 0.000223
%  Iteration: 51, Sup diff: 0.000144
%  Iteration: 61, Sup diff: 0.000093
%  Iteration: 71, Sup diff: 0.000060
%  Iteration: 81, Sup diff: 0.000039
%  Iteration: 91, Sup diff: 0.000026
%  Iteration: 101, Sup diff: 0.000017
%  Iteration: 111, Sup diff: 0.000011
%  Iteration: 121, Sup diff: 0.000007
%  Iteration: 131, Sup diff: 0.000005
%  Iteration: 141, Sup diff: 0.000003
%  Iteration: 151, Sup diff: 0.000002
%  Iteration: 161, Sup diff: 0.000001
% 历时 183.795166 秒。

% save ShashaWang_JFV_PS1_100_capital_grid_points_accelerator
%  Iteration: 41, Sup diff: 0.000854
%  Iteration: 51, Sup diff: 0.000568
%  Iteration: 61, Sup diff: 0.000378
%  Iteration: 71, Sup diff: 0.000251
%  Iteration: 81, Sup diff: 0.000167
%  Iteration: 91, Sup diff: 0.000111
%  Iteration: 101, Sup diff: 0.000074
%  Iteration: 111, Sup diff: 0.000049
%  Iteration: 121, Sup diff: 0.000033
%  Iteration: 131, Sup diff: 0.000022
%  Iteration: 141, Sup diff: 0.000014
%  Iteration: 151, Sup diff: 0.000010
%  Iteration: 161, Sup diff: 0.000006
%  Iteration: 171, Sup diff: 0.000004
%  Iteration: 181, Sup diff: 0.000003
%  Iteration: 191, Sup diff: 0.000002
%  Iteration: 201, Sup diff: 0.000001
% 历时 327.309912 秒。

% After I take advantage of the monotonicity of capital policy function wrt
% capital, the process speeds up a lot
%  Iteration:  1, Sup diff: 0.006262
%  Iteration: 11, Sup diff: 0.001010
%  Iteration: 21, Sup diff: 0.000705
%  Iteration: 31, Sup diff: 0.000379
%  Iteration: 41, Sup diff: 0.000242
%  Iteration: 51, Sup diff: 0.000157
%  Iteration: 61, Sup diff: 0.000102
%  Iteration: 71, Sup diff: 0.000066
%  Iteration: 81, Sup diff: 0.000043
%  Iteration: 91, Sup diff: 0.000028
%  Iteration: 101, Sup diff: 0.000019
%  Iteration: 111, Sup diff: 0.000012
%  Iteration: 121, Sup diff: 0.000008
%  Iteration: 131, Sup diff: 0.000005
%  Iteration: 141, Sup diff: 0.000003
%  Iteration: 151, Sup diff: 0.000002
%  Iteration: 161, Sup diff: 0.000002
%  Iteration: 171, Sup diff: 0.000001
% 历时 131.402012 秒。

save ShashaWang_JFV_PS1_250_capital_grid_points_accelerator
%  Iteration:  1, Sup diff: 0.006056
%  Iteration: 11, Sup diff: 0.000997
%  Iteration: 21, Sup diff: 0.000681
%  Iteration: 31, Sup diff: 0.000371
%  Iteration: 41, Sup diff: 0.000237
%  Iteration: 51, Sup diff: 0.000153
%  Iteration: 61, Sup diff: 0.000099
%  Iteration: 71, Sup diff: 0.000065
%  Iteration: 81, Sup diff: 0.000042
%  Iteration: 91, Sup diff: 0.000028
%  Iteration: 101, Sup diff: 0.000018
%  Iteration: 111, Sup diff: 0.000012
%  Iteration: 121, Sup diff: 0.000008
%  Iteration: 131, Sup diff: 0.000005
%  Iteration: 141, Sup diff: 0.000003
%  Iteration: 151, Sup diff: 0.000002
%  Iteration: 161, Sup diff: 0.000001
%  Iteration: 171, Sup diff: 0.000001
% 历时 392.471287 秒。

%% 7. Multigrid

%% 7.1 Multigrid withOUT Accelerator

% Use multigrid method to speed up iteration
% % I didn't choose [100 500 5000] because I choose Nk = 100 for fixed grid
% % instead of 250 as required.

kGridLength = [100,500,5000]; % number of points in grid for capital 
Nk = max(kGridLength);

%% Create cells to restore results
% For consistency I use cells, even though for iteration and mDifference I
% could use matrices.
% Now I realized I could've used cell to combine the two exogenous shocks

cValueMultigrid = cell(length(kGridLength),1);
cKPolicyMultigrid = cell(length(kGridLength),1);
cLaborPolicy_1Multigrid = cell(length(kGridLength),1);
cLaborPolicy_2Multigrid = cell(length(kGridLength),1);
cConsumptionPolicy_1Multigrid = cell(length(kGridLength),1);
cConsumptionPolicy_2Multigrid = cell(length(kGridLength),1);
cIterationMultigrid  = cell(length(kGridLength),1);
cDifferenceMultigrid  = cell(length(kGridLength),1);

%% Required matrices and vectors

mValue0        = utilitySteadyState.*ones(kGridLength(1),Na);
inputs.mValue0 = mValue0;
mValue         = zeros(kGridLength(1),Na);
mKPolicy        = zeros(kGridLength(1),Na);
mLaborPolicy_1 = zeros(kGridLength(1),Na);
mLaborPolicy_2 = zeros(kGridLength(1),Na); 
mConsumptionPolicy_1   = zeros(kGridLength(1),Na);
mConsumptionPolicy_2   = zeros(kGridLength(1),Na);

tic


for i=1:length(kGridLength)
    vGrid_kMultigrid           = linspace(kMin,kMax,kGridLength(i))';
    
    %% Compute efficiency matrices
    % Admittedly and ideally, I should compute efficiency matrices to retrieve from
    % later for EVERY grid number of k. But since it already takes more
    % than 40 minutes when Nk=250, it would take enormous amount of time
    % when Nk is larger. So I use the efficiency matrices for Nk=250, and 
    % use linear interpolation "interp2" to retrieve fsolve results.  
    
    % So I need to rewrite valueFunction.m as valueFunctionMultigrid.m to
    % change interp1 to interp2.
    
    % Same change much be made when filling in the policy function matrices.
            
    %% Main iteration
    maxIter = 10000;
    tolerance     = 1e-6;

    iteration       = 1;
    mDifference = zeros(maxIter,1);
    mDifference(iteration) = 100;

    options = optimset('Display', 'off');
    opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','off');
    % laborInitial=[labor_1_SteadyState,labor_2_SteadyState];

    

    % Finding value and policy functions numerically
    while iteration <= maxIter && mDifference(iteration) > tolerance
        expectedValue0 = mValue0 * mProb_a1a2'; % row: kPrime, column: a, i.e., expected value if today's shock is a and I choose k as tomorrow's capital stock
    %     inputs.expectedValue0 = expectedValue0;

        for ia = 1:Na
    %         a = mGrid_a1a2(ia,:);
            a_1 = mGrid_a1a2(ia,1);
            a_2 = mGrid_a1a2(ia,2);

            for ik = 1:kGridLength(i)
                k = vGrid_kMultigrid(ik);

                if ik ==1

                    [kPrime, fval] = fminbnd(@(kPrime) ...
                        -valueFunctionMultigrid(kPrime,k,ik,ia,a_1,a_2,expectedValue0,vGrid_kMultigrid,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                        vGrid_kMultigrid(1),min(1.2*vGrid_kMultigrid(ik),vGrid_kMultigrid(end)),options);       
                else
                    [kPrime, fval] = fminbnd(@(kPrime) ...
                        -valueFunctionMultigrid(kPrime,k,ik,ia,a_1,a_2,expectedValue0,vGrid_kMultigrid,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                        mKPolicy((ik-1),ia),min(1.2*vGrid_kMultigrid(ik),vGrid_kMultigrid(end)),options);
                end

                mKPolicy(ik,ia) = kPrime;
                mValue(ik,ia) = -fval;     

     %% If we can retrieve other outputs of the objective function of fminbnd, it would be much easier and FASTER to get the labor and consumption policy function
                if mDifference(iteration) < tolerance*10 % Store the policy functions in case it's the last iteration
    %                 vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
    %                 mLaborPolicy_1(ik,ia) = vLabor(1);
    %                 mLaborPolicy_2(ik,ia) = vLabor(2);
    %                 mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
    %                 mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,vLabor(2));
    %                 laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
    %                 inputs.laborInitial = laborInitial;
                    mLaborPolicy_1(ik,ia) = interp2(vGrid_k,vGrid_k,mLabor_1Fsolve(:,:,ia),kPrime,k);
                    mLaborPolicy_2(ik,ia) = interp2(vGrid_k,vGrid_k,mLabor_2Fsolve(:,:,ia),kPrime,k);
                    mConsumptionPolicy_1(ik,ia) = interp2(vGrid_k,vGrid_k,mConsumption_1Fsolve(:,:,ia),kPrime,k);
                    mConsumptionPolicy_2(ik,ia) = interp2(vGrid_k,vGrid_k,mConsumption_2Fsolve(:,:,ia),kPrime,k);
                end
            end
        end

        iteration = iteration + 1;
        mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
        mValue0         = mValue;

%         if mod(iteration,10) == 2
            fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 
%         end
        
    %     inputs.valueFunction0 = mValue0;
    %     inputs.laborFunction  = mLaborPolicy_1;
    %     policyFunction0       = mKPolicy;
    end
    
    fprintf(' Convergence Achieved. \n Number of Grid Points: %2.0f\n Iteration: %2.0f, Sup difference: %2.8f\n ', kGridLength(i),iteration-1, mDifference(iteration)); 
    
    % Store result for this grid point
    cValueMultigrid{i}= mValue;
    cKPolicyMultigrid{i} = mKPolicy;
    cLaborPolicy_1Multigrid{i} = mLaborPolicy_1;
    cLaborPolicy_2Multigrid{i} = mLaborPolicy_2;
    cConsumptionPolicy_1Multigrid{i} = mConsumptionPolicy_1;
    cConsumptionPolicy_2Multigrid{i} = mConsumptionPolicy_2;
    cIteration{i} = iteration - 1;
    cDifference{i} = mDifference(2:iteration);

    if i ~= length(kGridLength)
        mValue0 = interp1(vGrid_kMultigrid, mValue,linspace(kMin, kMax, kGridLength(i+1))');% 这里不知道linspace后面要不要加'变为列向量。试试吧
        mValue  = mValue0; % Just for consistency of dimensionality. mValue can be anything you want, as long as the dimensionality is right
    %         kPolicy = interp1(grid_k,kPolicy,linspace(kMin, kMax, kGridLength(i+1)));
        mKPolicy         = zeros(kGridLength(i+1),Na);
        mLaborPolicy_1         = zeros(kGridLength(i+1),Na);
        mLaborPolicy_2         = zeros(kGridLength(i+1),Na);
        mConsumptionPolicy_1         = zeros(kGridLength(i+1),Na);
        mConsumptionPolicy_2         = zeros(kGridLength(i+1),Na);
        
        
    end
    
end

toc

%% For accuracy test, compute the euler equation error
errorEulerEquationLinearInterpolationAccelerator = eulerEquationErrorFunction(Nk,vGrid_kMultigrid,mKPolicy,mLaborPolicy_1,mConsumptionPolicy_1,mConsumptionPolicy_2,Na,mGrid_a1a2,mProb_a1a2,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL);
errorEulerEquationLinearInterpolationAcceleratorDecimalLog = log10( errorEulerEquationLinearInterpolationAccelerator );

figure;
[kk,aa]=meshgrid(vGrid_kMultigrid, mGrid_a1a2(:,1));
mesh(kk, aa, errorEulerEquationLinearInterpolationAcceleratorDecimalLog');

title('Euler Equation Error $log_{10}$ Linear Interpolation - Multigrid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('error','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_eulerEquationErrorLinearInterpolation_multigrid')

save ShashaWang_JFV_PS1_250_capital_grid_points_multigrid

%% Figures for multigrid WITHOUT accelerator
figure;
mesh(kk, aa, mValue');

title('Value - Multigrid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Value','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_value_multigrid')

figure;
mesh(kk, aa, mKPolicy');

title('Policy for Next Period Capital - Multigrid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_kPolicy_multigrid')

figure;
mesh(kk, aa, mLaborPolicy_1');

title('Policy for Good 1 Labor - Multigrid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Labor','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_1_multigrid')

figure;
mesh(kk, aa, mLaborPolicy_2');

title('Policy for Good 2 Labor - Multigrid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Labor','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_2_multigrid')

figure;
mesh(kk, aa, mConsumptionPolicy_1');

title('Policy for Good 1 Consumption - Multigrid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Consumption','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_1_multigrid')

figure;
mesh(kk, aa, mConsumptionPolicy_2');

title('Policy for Good 2 Consumption - Multigrid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Consumption','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_2_multigrid')

%% 7.2 Multigrid with Accelerator
% Use multigrid method to speed up iteration
% % I didn't choose [100 500 5000] because I choose Nk = 100 for fixed grid
% % instead of 250 as required.

kGridLength = [100,500,5000]; % number of points in grid for capital 
Nk = max(kGridLength);

%% Create cells to restore results
% For consistency I use cells, even though for iteration and mDifference I
% could use matrices.
% Now I realized I could've used cell to combine the two exogenous shocks

cValueMultigrid = cell(length(kGridLength),1);
cKPolicyMultigrid = cell(length(kGridLength),1);
cLaborPolicy_1Multigrid = cell(length(kGridLength),1);
cLaborPolicy_2Multigrid = cell(length(kGridLength),1);
cConsumptionPolicy_1Multigrid = cell(length(kGridLength),1);
cConsumptionPolicy_2Multigrid = cell(length(kGridLength),1);
cIterationMultigrid  = cell(length(kGridLength),1);
cDifferenceMultigrid  = cell(length(kGridLength),1);

%% Required matrices and vectors

mValue0        = utilitySteadyState.*ones(kGridLength(1),Na);
inputs.mValue0 = mValue0;
mValue         = zeros(kGridLength(1),Na);
mKPolicy        = zeros(kGridLength(1),Na);
mLaborPolicy_1 = zeros(kGridLength(1),Na);
mLaborPolicy_2 = zeros(kGridLength(1),Na); 
mConsumptionPolicy_1   = zeros(kGridLength(1),Na);
mConsumptionPolicy_2   = zeros(kGridLength(1),Na);

tic

for i=1:length(kGridLength)
    vGrid_kMultigrid           = linspace(kMin,kMax,kGridLength(i))';
    
    %% Compute efficiency matrices
    % Admittedly and ideally, I should compute efficiency matrices to retrieve from
    % later for EVERY grid number of k. But since it already takes more
    % than 40 minutes when Nk=250, it would take enormous amount of time
    % when Nk is larger. So I use the efficiency matrices for Nk=250, and 
    % use linear interpolation "interp2" to retrieve fsolve results.  
    
    % So I need to rewrite valueFunction.m as valueFunctionMultigrid.m to
    % change interp1 to interp2.
    
    % Same change much be made when filling in the policy function matrices.
            
    %% Main iteration
    maxIter = 10000;
    tolerance     = 1e-6;

    iteration       = 1;
    mDifference = zeros(maxIter,1);
    mDifference(iteration) = 100;

    options = optimset('Display', 'off');
    opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','off');
    % laborInitial=[labor_1_SteadyState,labor_2_SteadyState];

    tic

    % Finding value and policy functions numerically
    while iteration <= maxIter ...
            && ( (mDifference(iteration) > tolerance)  | (mDifference(iteration) <= tolerance && mod(iteration,10)~=2)) % make sure the last iteration does the maximization
        
        if mDifference(iteration) > tolerance % do the value function iteration with the accelerator
            
            expectedValue0 = mValue0 * mProb_a1a2'; % row: kPrime, column: a, i.e., expected value if today's shock is a and I choose k as tomorrow's capital stock
        %     inputs.expectedValue0 = expectedValue0;

            for ia = 1:Na
        %         a = mGrid_a1a2(ia,:);
                a_1 = mGrid_a1a2(ia,1);
                a_2 = mGrid_a1a2(ia,2);

                for ik = 1:kGridLength(i)
                    k = vGrid_kMultigrid(ik);

                    if mod(iteration,10) == 1

                        if ik ==1

                            [kPrime, fval] = fminbnd(@(kPrime) ...
                                -valueFunctionMultigrid(kPrime,k,ik,ia,a_1,a_2,expectedValue0,vGrid_kMultigrid,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                                vGrid_kMultigrid(1),min(1.2*vGrid_kMultigrid(ik),vGrid_kMultigrid(end)),options);       
                        else
                            [kPrime, fval] = fminbnd(@(kPrime) ...
                                -valueFunctionMultigrid(kPrime,k,ik,ia,a_1,a_2,expectedValue0,vGrid_kMultigrid,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                                mKPolicy((ik-1),ia),min(1.2*vGrid_kMultigrid(ik),vGrid_kMultigrid(end)),options);
                        end

                        % Policy function has to be filled in every 10th
                        % time 
                        mKPolicy(ik,ia) = kPrime;
                        mValue(ik,ia) = -fval;     
                        mLaborPolicy_1(ik,ia) = interp2(vGrid_k,vGrid_k,mLabor_1Fsolve(:,:,ia),kPrime,k);
                        mLaborPolicy_2(ik,ia) = interp2(vGrid_k,vGrid_k,mLabor_2Fsolve(:,:,ia),kPrime,k);
                        mConsumptionPolicy_1(ik,ia) = interp2(vGrid_k,vGrid_k,mConsumption_1Fsolve(:,:,ia),kPrime,k);
                        mConsumptionPolicy_2(ik,ia) = interp2(vGrid_k,vGrid_k,mConsumption_2Fsolve(:,:,ia),kPrime,k);
                    else
                        currentUtility = interp2(vGrid_k,vGrid_k,mCurrentUtilityFsolve(:,:,ia),mKPolicy(ik,ia),k);
                        expectedValue = interp1(vGrid_kMultigrid,expectedValue0(:,ia),mKPolicy(ik,ia));
                        value = (1-bbeta)*currentUtility + bbeta * expectedValue;
                        mValue(ik,ia) = value;
                    end
                end
            end

            iteration = iteration + 1;
            mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
            mValue0         = mValue;

    %         if mod(iteration,10) == 2
                fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 
    %         end
        else % do the value function iteration without the accelerator
            expectedValue0 = mValue0 * mProb_a1a2'; % row: kPrime, column: a, i.e., expected value if today's shock is a and I choose k as tomorrow's capital stock

            for ia = 1:Na
                a_1 = mGrid_a1a2(ia,1);
                a_2 = mGrid_a1a2(ia,2);

                for ik = 1:kGridLength(i)
                    k = vGrid_kMultigrid(ik);


                    if ik ==1

                        [kPrime, fval] = fminbnd(@(kPrime) ...
                            -valueFunctionMultigrid(kPrime,k,ik,ia,a_1,a_2,expectedValue0,vGrid_kMultigrid,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                            vGrid_kMultigrid(1),min(1.2*vGrid_kMultigrid(ik),vGrid_kMultigrid(end)),options);       
                    else
                        [kPrime, fval] = fminbnd(@(kPrime) ...
                            -valueFunctionMultigrid(kPrime,k,ik,ia,a_1,a_2,expectedValue0,vGrid_kMultigrid,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                            mKPolicy((ik-1),ia),min(1.2*vGrid_kMultigrid(ik),vGrid_kMultigrid(end)),options);
                    end

                    mKPolicy(ik,ia) = kPrime;
                    mValue(ik,ia) = -fval;     
                    mLaborPolicy_1(ik,ia) = interp2(vGrid_k,vGrid_k,mLabor_1Fsolve(:,:,ia),kPrime,k);
                    mLaborPolicy_2(ik,ia) = interp2(vGrid_k,vGrid_k,mLabor_2Fsolve(:,:,ia),kPrime,k);
                    mConsumptionPolicy_1(ik,ia) = interp2(vGrid_k,vGrid_k,mConsumption_1Fsolve(:,:,ia),kPrime,k);
                    mConsumptionPolicy_2(ik,ia) = interp2(vGrid_k,vGrid_k,mConsumption_2Fsolve(:,:,ia),kPrime,k);
                end
            end
            iteration = iteration + 1;
            mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
            mValue0         = mValue;

            fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 
            
            
            break
        
        end
        
        

    end
    
    fprintf(' Convergence Achieved. \n Number of Grid Points: %2.0f\n Iteration: %2.0f, Sup difference: %2.8f\n ', kGridLength(i),iteration-1, mDifference(iteration)); 
    
    % Store result for this grid point
    cValueMultigrid{i}= mValue;
    cKPolicyMultigrid{i} = mKPolicy;
    cLaborPolicy_1Multigrid{i} = mLaborPolicy_1;
    cLaborPolicy_2Multigrid{i} = mLaborPolicy_2;
    cConsumptionPolicy_1Multigrid{i} = mConsumptionPolicy_1;
    cConsumptionPolicy_2Multigrid{i} = mConsumptionPolicy_2;
    cIteration{i} = iteration - 1;
    cDifference{i} = mDifference(2:iteration);

    if i ~= length(kGridLength)
        mValue0 = interp1(vGrid_kMultigrid, mValue,linspace(kMin, kMax, kGridLength(i+1))');
        mValue  = mValue0; % Just for consistency of dimensionality. mValue can be anything you want, as long as the dimensionality is right
    %         kPolicy = interp1(grid_k,kPolicy,linspace(kMin, kMax, kGridLength(i+1)));
        mKPolicy         = zeros(kGridLength(i+1),Na);
        mLaborPolicy_1         = zeros(kGridLength(i+1),Na);
        mLaborPolicy_2         = zeros(kGridLength(i+1),Na);
        mConsumptionPolicy_1         = zeros(kGridLength(i+1),Na);
        mConsumptionPolicy_2         = zeros(kGridLength(i+1),Na);
        
        
    end
    
end

toc

%% For accuracy test, compute the euler equation error
errorEulerEquationLinearInterpolationAccelerator = eulerEquationErrorFunction(Nk,vGrid_kMultigrid,mKPolicy,mLaborPolicy_1,mConsumptionPolicy_1,mConsumptionPolicy_2,Na,mGrid_a1a2,mProb_a1a2,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL);
errorEulerEquationLinearInterpolationAcceleratorDecimalLog = log10( errorEulerEquationLinearInterpolationAccelerator );

[kk,aa]=meshgrid(vGrid_kMultigrid, mGrid_a1a2(:,1));

figure;
mesh(kk, aa, errorEulerEquationLinearInterpolationAcceleratorDecimalLog');

title('Euler Equation Error $log_{10}$ Linear Interpolation - Multigrid & Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('error','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_eulerEquationErrorLinearInterpolation_multigrid&accelerator')

save ShashaWang_JFV_PS1_250_capital_grid_points_multigrid&accelerator

%% Figures for multigrid WITH accelerator
figure;
mesh(kk, aa, mValue');

title('Value - Multigrid & Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Value','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_value_multigrid&accelerator')

figure;
mesh(kk, aa, mKPolicy');

title('Policy for Next Period Capital - Multigrid & Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_kPolicy_multigrid&accelerator')

figure;
mesh(kk, aa, mLaborPolicy_1');

title('Policy for Good 1 Labor - Multigrid & Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Labor','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_1_multigrid&accelerator')

figure;
mesh(kk, aa, mLaborPolicy_2');

title('Policy for Good 2 Labor - Multigrid & Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Labor','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_2_multigrid&accelerator')

figure;
mesh(kk, aa, mConsumptionPolicy_1');

title('Policy for Good 1 Consumption - Multigrid & Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Consumption','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_1_multigrid&accelerator')

figure;
mesh(kk, aa, mConsumptionPolicy_2');

title('Policy for Good 2 Consumption - Multigrid & Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Consumption','interpreter','latex')
xlim([min(vGrid_kMultigrid),max(vGrid_kMultigrid)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_2_multigrid&accelerator')
