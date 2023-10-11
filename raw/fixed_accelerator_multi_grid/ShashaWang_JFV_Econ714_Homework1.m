% Homework 1 Econ714 JFV
% A two good production economy
% Written by Shasha Wang

close all;
clear;

cd 'E:\Dropbox\fall 19-20\jesus\homework_1\set_labor_to_steady_state_first_then_standard_iteration';

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
Nk = 50;
vGrid_k = linspace(kMin,kMax,Nk)';
inputs.vGrid_k = vGrid_k;

% We can see from the question that once given a1,a2,k,and chosen k',
% labor_1, labor_2 are both chosen by optimization, and consumption_1,
% consumption_2 are both chosen immediately.

% So we need to write a function laborFunction.m from a1,a2,k,k' to labor_1,labor_2
% And we need to write consumption functions consumptionFunction1.m,
% consumptionFunction2.m

% consumptionFunction1 = @(a_1,k,kPrime,labor_1,aalphaK,aalphaL,ddelta)a_1 * k^aalphaK * labor_1^aalphaL + (1-ddelta) * k - kPrime;
% consumptionFunction2 = @(a_2,labor_2)a_2 * labor_2;

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
laborInitial=[labor_1_SteadyState,labor_2_SteadyState];

tic

% Finding value and policy functions numerically
while iteration <= maxIter && mDifference(iteration) > tolerance
    expectedValue0 = mValue0 * mProb_a1a2';
%     inputs.expectedValue0 = expectedValue0;

    for ia = 1:Na
%         a = mGrid_a1a2(ia,:);
        a_1 = mGrid_a1a2(ia,1);
        a_2 = mGrid_a1a2(ia,2);
        
        for ik = 1:Nk
            k = vGrid_k(ik);
                        
            [kPrime, vAux] = fminbnd(@(kPrime) ...
                -valueFunction(kPrime,k,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                max(0.8*vGrid_k(ik),vGrid_k(1)),min(1.2*vGrid_k(ik),vGrid_k(end)),options);
            mKPolicy(ik,ia) = kPrime;
            mValue(ik,ia) = -vAux;            
            
 %% If we can retrieve other outputs of the objective function of fminbnd, it would be much easier and FASTER to get the labor and consumption policy function
            if mDifference(iteration) < tolerance*10 % Store the policy functions in case it's the last iteration
                vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
                mLaborPolicy_1(ik,ia) = vLabor(1);
                mLaborPolicy_2(ik,ia) = vLabor(2);
                mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
                mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,vLabor(2));
                laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
%                 inputs.laborInitial = laborInitial;

            end
        end
    end
    
    iteration = iteration + 1;
    mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
    mValue0         = mValue;
    
    if mod(iteration,10) == 2
        fprintf(' Iteration: %2.0f, Sup diff: %2.6f\n', iteration-1, mDifference(iteration)); 
    end
    
%     inputs.valueFunction0 = mValue0;
%     inputs.laborFunction  = mLaborPolicy_1;
%     policyFunction0       = mKPolicy;
end

toc

%% figures for Value Function Iteration with a Fixed Grid

figure(2);
[kk,aa]=meshgrid(vGrid_k, mGrid_a1a2(:,1));
mesh(kk, aa, mValue');

title('Value Under Different Shocks','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Value','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_value_3D')

figure(3);
mesh(kk, aa, mKPolicy');

title('Policy for Next Period Capital Under Different Shocks','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_kPolicy_3D')

figure(4);
mesh(kk, aa, mLaborPolicy_1');

title('Policy for Good 1 Labor Under Different Shock','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_1_3D')

figure(5);
mesh(kk, aa, mLaborPolicy_2');

title('Policy for Good 2 Labor Under Different Shock','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_2_3D')

figure(6);
mesh(kk, aa, mConsumptionPolicy_1');

title('Policy for Good 1 Consumption Under Different Shock','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_1_3D')

figure(7);
mesh(kk, aa, mConsumptionPolicy_2');

title('Policy for Good 2 Consumption Under Different Shocks','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_2_3D')

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
% 历时 6566.472522 秒。

save ShashaWang_JFV_PS1_50_capital_grid_points
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
% 历时 25248.327232 秒。

%% For accuracy test, compute the euler equation error

% First let's not use linear interpolation for tomorrow's values

leftHandSideOfEulerEquation = zeros(Nk,Na); % Marginal utility of today
rightHandSideOfEulerEquation = zeros(Nk,Na); % expected Marginal utility of tomorrow
% laborInitial=[labor_1_SteadyState,labor_2_SteadyState];

% Compute marginal utility of today
for ia = 1:Na
    a = mGrid_a1a2(ia,:);
    a_1 = mGrid_a1a2(ia,1);
    a_2 = mGrid_a1a2(ia,2);        
    
    for ik = 1:Nk
%         k = vGrid_k(ik);
        kPrime = mKPolicy(ik,ia);
%         kPrimePrime = depend on aPrime
%         vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
        [v,ikPrime]=min(abs(kPrime - vGrid_k));
        
        consumption_1 = mConsumptionPolicy_1(ik,ia);
        consumption_2 = mConsumptionPolicy_2(ik,ia);
%         labor_1 = mLaborPolicy_1(ik,ia);
%         labor_2 = mLaborPolicy_2(ik,ia);
%         laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
%         leftHandSideOfEulerEquation(ik,ia) = mmu_1 * (a_2 * labor_2)^mmu_2 * ...
%             (a_1 * k^aalphaK * labor_1^aalphaL + (1-ddelta) * k - kPrime)^(mmu_1-1);
       
        marginalUtilityToday = mmu_1 * (consumption_2)^mmu_2 * (consumption_1)^(mmu_1-1);
        leftHandSideOfEulerEquation(ik,ia) = marginalUtilityToday;
        
% Compute expected marginal utility of tomorrow
        expected = zeros(Na,1);
        
        for iaPrime = 1:Na
            aPrime = mGrid_a1a2(iaPrime,:);
            aPrime_1 = mGrid_a1a2(iaPrime,1);
            aPrime_2 = mGrid_a1a2(iaPrime,2); 
            
%             kPrimePrime = mKPolicy(ikPrime,iaPrime);
            consumptionPrime_1 = mConsumptionPolicy_1(ikPrime,iaPrime);
            consumptionPrime_2 = mConsumptionPolicy_2(ikPrime,iaPrime);
            laborPrime_1 = mLaborPolicy_1(ikPrime,iaPrime);
            
            marginalUtilityTomorrow = mmu_1 * (consumptionPrime_2)^mmu_2 * (consumptionPrime_1)^(mmu_1-1);
            returnOnCapital = 1 - ddelta + aPrime_1 * laborPrime_1^aalphaL * aalphaK * (kPrime)^(aalphaK-1);
            
            unexpected = bbeta * marginalUtilityTomorrow * returnOnCapital;
            expected(iaPrime) = unexpected * mProb_a1a2(ia,iaPrime);
        end
        rightHandSideOfEulerEquation(ik,ia) = sum(expected);
    end
end

errorEulerEquationNoLinearInterpolation = abs(leftHandSideOfEulerEquation - rightHandSideOfEulerEquation);
errorEulerEquationNoLinearInterpolationDecimalLog = log10( errorEulerEquationNoLinearInterpolation );

figure(8);
mesh(kk, aa, errorEulerEquationNoLinearInterpolationDecimalLog');

title('Euler Equation Error $log_{10}$ No Linear Interpolation','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('error','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_eulerEquationErrorNoLinearInterpolation_3D')

% Then let's use LINEAR INTERPOLATION for tomorrow's values

errorEulerEquationLinearInterpolation = eulerEquationErrorFunction(Nk,vGrid_k,mKPolicy,mLaborPolicy_1,mConsumptionPolicy_1,mConsumptionPolicy_2,Na,mGrid_a1a2,mProb_a1a2,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL);
errorEulerEquationLinearInterpolationDecimalLog = log10( errorEulerEquationLinearInterpolation );

figure(9);
mesh(kk, aa, errorEulerEquationLinearInterpolationDecimalLog');

title('Euler Equation Error $log_{10}$ Linear Interpolation','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('error','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_eulerEquationErrorLinearInterpolation_3D')













% %% 6. Value function iteration using Accelerator
% mValue0        = utilitySteadyState.*ones(Nk,Na);
% inputs.mValue0 = mValue0;
% mValue         = zeros(Nk,Na);
% mKPolicy        = zeros(Nk,Na);
% mLaborPolicy_1 = zeros(Nk,Na);
% mLaborPolicy_2 = zeros(Nk,Na); 
% mConsumptionPolicy_1   = zeros(Nk,Na);
% mConsumptionPolicy_2   = zeros(Nk,Na);
% 
% %% Main iteration
% maxIter = 10000;
% tolerance     = 1e-6;
% 
% iteration       = 1;
% mDifference = zeros(maxIter,1);
% mDifference(iteration) = 100;
% 
% options = optimset('Display', 'off');
% opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','off');
% laborInitial=[labor_1_SteadyState,labor_2_SteadyState];
% 
% tic
% 
% % Finding value and policy functions numerically
% while iteration <= maxIter && mDifference(iteration) > tolerance
%     expectedValue0 = mValue0 * mProb_a1a2';
% %     inputs.expectedValue0 = expectedValue0;
% 
%     for ia = 1:Na
% %         a = mGrid_a1a2(ia,:);
%         a_1 = mGrid_a1a2(ia,1);
%         a_2 = mGrid_a1a2(ia,2);
%         
%         for ik = 1:Nk
%             k = vGrid_k(ik);
%                         
%             [kPrime, vAux] = fminbnd(@(kPrime) ...
%                 -valueFunction(kPrime,k,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
%                 max(0.8*vGrid_k(ik),vGrid_k(1)),min(1.2*vGrid_k(ik),vGrid_k(end)),options);
%             mKPolicy(ik,ia) = kPrime;
%             mValue(ik,ia) = -vAux;            
%             
%  %% If we can retrieve other outputs of the objective function of fminbnd, it would be much easier and FASTER to get the labor and consumption policy function
%             if mDifference(iteration) < tolerance*10 % Store the policy functions in case it's the last iteration
%                 vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
%                 mLaborPolicy_1(ik,ia) = vLabor(1);
%                 mLaborPolicy_2(ik,ia) = vLabor(2);
%                 mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
%                 mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,vLabor(2));
%                 laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
% %                 inputs.laborInitial = laborInitial;
% 
%             end
%         end
%     end
%     
%     iteration = iteration + 1;
%     mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
%     mValue0         = mValue;
%     
%     if mod(iteration,10) == 2
%         fprintf(' Iteration: %2.0f, Sup diff: %2.6f\n', iteration-1, mDifference(iteration)); 
%     end
%     
% %     inputs.valueFunction0 = mValue0;
% %     inputs.laborFunction  = mLaborPolicy_1;
% %     policyFunction0       = mKPolicy;
% end
% 
% toc
% 
% %% figures for Value Function Iteration with Accelerator
% 
% figure;
% [kk,aa]=meshgrid(vGrid_k, mGrid_a1a2(:,1));
% mesh(kk, aa, mValue');
% 
% title('Value Under Different Shocks','interpreter','latex')
% xlabel('Capital Stock $k$','interpreter','latex')
% ylabel('shocks $z_1$ $z_2$','interpreter','latex')
% zlabel('Value','interpreter','latex')
% xlim([min(vGrid_k),max(vGrid_k)])
% ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
% savefig('q3_value_3D')
% 
% figure;
% mesh(kk, aa, mKPolicy');
% 
% title('Policy for Next Period Capital Under Different Shocks','interpreter','latex')
% xlabel('Capital Stock $k$','interpreter','latex')
% ylabel('shocks $z_1$ $z_2$','interpreter','latex')
% zlabel('Next Period Capital $k\prime$','interpreter','latex')
% xlim([min(vGrid_k),max(vGrid_k)])
% ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
% savefig('q3_kPolicy_3D')
% 
% figure;
% mesh(kk, aa, mLaborPolicy_1');
% 
% title('Policy for Good 1 Labor Under Different Shock','interpreter','latex')
% xlabel('Capital Stock $k$','interpreter','latex')
% ylabel('shocks $z_1$ $z_2$','interpreter','latex')
% zlabel('Good 1 Labor','interpreter','latex')
% xlim([min(vGrid_k),max(vGrid_k)])
% ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
% savefig('q3_laborPolicy_1_3D')
% 
% figure;
% mesh(kk, aa, mLaborPolicy_2');
% 
% title('Policy for Good 2 Labor Under Different Shock','interpreter','latex')
% xlabel('Capital Stock $k$','interpreter','latex')
% ylabel('shocks $z_1$ $z_2$','interpreter','latex')
% zlabel('Good 2 Labor','interpreter','latex')
% xlim([min(vGrid_k),max(vGrid_k)])
% ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
% savefig('q3_laborPolicy_2_3D')
% 
% figure;
% mesh(kk, aa, mConsumptionPolicy_1');
% 
% title('Policy for Good 1 Consumption Under Different Shock','interpreter','latex')
% xlabel('Capital Stock $k$','interpreter','latex')
% ylabel('shocks $z_1$ $z_2$','interpreter','latex')
% zlabel('Good 1 Consumption','interpreter','latex')
% xlim([min(vGrid_k),max(vGrid_k)])
% ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
% savefig('q3_consumptionPolicy_1_3D')
% 
% figure;
% mesh(kk, aa, mConsumptionPolicy_2');
% 
% title('Policy for Good 2 Consumption Under Different Shocks','interpreter','latex')
% xlabel('Capital Stock $k$','interpreter','latex')
% ylabel('shocks $z_1$ $z_2$','interpreter','latex')
% zlabel('Good 2 Consumption','interpreter','latex')
% xlim([min(vGrid_k),max(vGrid_k)])
% ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
% savefig('q3_consumptionPolicy_2_3D')
% 
% save ShashaWang_JFV_PS1_10_capital_grid_points_accelerator
% 
% %% For accuracy test, compute the euler equation error
% errorEulerEquationLinearInterpolation = eulerEquationErrorFunction(Nk,vGrid_k,mKPolicy,mLaborPolicy_1,mConsumptionPolicy_1,mConsumptionPolicy_2,Na,mGrid_a1a2,mProb_a1a2,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL);
% errorEulerEquationLinearInterpolationDecimalLog = log10( errorEulerEquationLinearInterpolation );
% 
% figure;
% mesh(kk, aa, errorEulerEquationLinearInterpolationDecimalLog');
% 
% title('Euler Equation Error $log_{10}$ Linear Interpolation','interpreter','latex')
% xlabel('Capital Stock $k$','interpreter','latex')
% ylabel('shocks $z_1$ $z_2$','interpreter','latex')
% zlabel('error','interpreter','latex')
% xlim([min(vGrid_k),max(vGrid_k)])
% ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
% savefig('q3_eulerEquationErrorLinearInterpolation_3D')
