% Homework 1 Econ714 JFV
% A two good production economy
% Written by Shasha Wang

main_setup % parameters, shocks

%% Value function iteration with a Fixed Grid

%% Step 1: Compute the Steady State

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

consumption_1_SteadyState = consumptionFunction1(1,kSteadyState,kSteadyState,labor_1_SteadyState,aalphaK,aalphaL,ddelta);
consumption_2_SteadyState = consumptionFunction2(1,labor_2_SteadyState);
utilitySteadyState = utilityFunction(consumption_1_SteadyState,consumption_2_SteadyState,labor_1_SteadyState,labor_2_SteadyState,mmu_1,mmu_2);

T = table(kSteadyState,labor_1_SteadyState,labor_2_SteadyState,consumption_1_SteadyState,consumption_2_SteadyState)

%% Step 2: set up k grid around the steady state
kSpread = 0.3;
kMax = kSteadyState * (1 + kSpread);
kMin = kSteadyState * (1 - kSpread);
Nk = 250;
vGrid_k = linspace(kMin,kMax,Nk)';
inputs.vGrid_k = vGrid_k;

%% Step 3. Pre-Iteration: Set labor to steady state

% Iterate on the Value function implied by the Social Planner¡¯s Problem using linear interpolation
% until the change in the sup norm between to iterations is less than 1e-6.
% Compute the Policy function.
% Describe the responses of the economy to a technology shock.

% We can see from the question that once given a1,a2,k,and chosen k',
% labor_1, labor_2 are both chosen by optimization, and consumption_1,
% consumption_2 are both chosen immediately.

% So we need to write a function laborFunction.m from a1,a2,k,k' to labor_1,labor_2
% And we need to write consumption functions consumptionFunction1.m,
% consumptionFunction2.m

% consumptionFunction1 = @(a_1,k,kPrime,labor_1,aalphaK,aalphaL,ddelta)a_1 * k^aalphaK * labor_1^aalphaL + (1-ddelta) * k - kPrime;
% consumptionFunction2 = @(a_2,labor_2)a_2 * labor_2;

%% Required matrices and vectors

mValue0                 = utilitySteadyState.*ones(Nk,Na); inputs.mValue0 = mValue0;
mValue                  = zeros(Nk,Na);
mKPolicy                = zeros(Nk,Na);
% mLaborPolicy_1          = zeros(Nk,Na);
% mLaborPolicy_2          = zeros(Nk,Na); 
mConsumptionPolicy_1    = zeros(Nk,Na);
mConsumptionPolicy_2    = zeros(Nk,Na);

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

    for ia = 1:Na
        a_1 = mGrid_a1a2(ia,1);
        a_2 = mGrid_a1a2(ia,2);
        
        for ik = 1:Nk
            k = vGrid_k(ik);
            
            if ik == 1
                        
                [kPrime, vAux] = fminbnd(@(kPrime) ...
                    -valueFunction_setLaborToSteadyState(kPrime,ik,k,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL,labor_1_SteadyState,labor_2_SteadyState),...
                    vGrid_k(1),min(1.2*vGrid_k(ik),vGrid_k(end)),options);
        
            else
                [kPrime, vAux] = fminbnd(@(kPrime) ...
                    -valueFunction_setLaborToSteadyState(kPrime,ik,k,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL,labor_1_SteadyState,labor_2_SteadyState),...
                    mKPolicy((ik-1),ia),min(1.2*vGrid_k(ik),vGrid_k(end)),options);
            end
            
            mKPolicy(ik,ia) = kPrime;
            mValue(ik,ia) = -vAux;       
            
 %% If we can retrieve other outputs of the objective function of fminbnd, it would be much easier and FASTER to get the labor and consumption policy function
            if mDifference(iteration) < tolerance*10 % Store the policy functions in case it's the last iteration
%                 vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
%                 mLaborPolicy_1(ik,ia) = vLabor(1);
%                 mLaborPolicy_2(ik,ia) = vLabor(2);
                mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,labor_1_SteadyState,aalphaK,aalphaL,ddelta);
                mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,labor_2_SteadyState);
%                 laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
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
save ShashaWang_JFV_PS1_250_capital_grid_points_valueFunctionIteration_setLaborToSteadyState

%  Iteration:  1, Sup diff: 0.005169
%  Iteration: 11, Sup diff: 0.000334
%  Iteration: 21, Sup diff: 0.000173
%  Iteration: 31, Sup diff: 0.000100
%  Iteration: 41, Sup diff: 0.000060
%  Iteration: 51, Sup diff: 0.000037
%  Iteration: 61, Sup diff: 0.000023
%  Iteration: 71, Sup diff: 0.000014
%  Iteration: 81, Sup diff: 0.000009
%  Iteration: 91, Sup diff: 0.000006
%  Iteration: 101, Sup diff: 0.000004
%  Iteration: 111, Sup diff: 0.000002
%  Iteration: 121, Sup diff: 0.000001
%  Iteration: 129, Sup diff: 0.000001
% ÀúÊ± 829.203325 Ãë¡£

%% Then do the regular Value Function Iteration using value function calculated above as the first guess

%% Required matrices and vectors

% mValue0                 = utilitySteadyState.*ones(Nk,Na); inputs.mValue0 = mValue0;
mValue                  = zeros(Nk,Na);
mKPolicy                = zeros(Nk,Na);
mLaborPolicy_1          = zeros(Nk,Na);
mLaborPolicy_2          = zeros(Nk,Na); 
mConsumptionPolicy_1    = zeros(Nk,Na);
mConsumptionPolicy_2    = zeros(Nk,Na);

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

    for ia = 1:Na
        a_1 = mGrid_a1a2(ia,1);
        a_2 = mGrid_a1a2(ia,2);
        
        for ik = 1:Nk
            k = vGrid_k(ik);
            
            if ik == 1
                [kPrime, vAux] = fminbnd(@(kPrime) ...
                    -valueFunction(kPrime,ik,k,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                    vGrid_k(1),min(1.2*vGrid_k(ik),vGrid_k(end)),options);
       
            else
                [kPrime, vAux] = fminbnd(@(kPrime) ...
                    -valueFunction(kPrime,ik,k,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                    mKPolicy((ik-1),ia),min(1.2*vGrid_k(ik),vGrid_k(end)),options);
            end
            
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
    
%     if mod(iteration,10) == 2
        fprintf(' Iteration: %2.0f, Sup diff: %2.6f\n', iteration-1, mDifference(iteration)); 
%     end
    
end

toc


%% For accuracy test, compute the euler equation error

% Let's use LINEAR INTERPOLATION for tomorrow's values

errorEulerEquationLinearInterpolation = eulerEquationErrorFunction(Nk,vGrid_k,mKPolicy,mLaborPolicy_1,mConsumptionPolicy_1,mConsumptionPolicy_2,Na,mGrid_a1a2,mProb_a1a2,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL);
errorEulerEquationLinearInterpolationDecimalLog = log10( errorEulerEquationLinearInterpolation );

[kk,aa]=meshgrid(vGrid_k, mGrid_a1a2(:,1));

figure
mesh(kk, aa, errorEulerEquationLinearInterpolationDecimalLog');

title('Euler Equation Error $log_{10}$ Linear Interpolation - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('error','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_eulerEquationErrorLinearInterpolation_fixed_grid')

save ShashaWang_JFV_PS1_250_capital_grid_points_valueFunctionIteration_setLaborToSteadyState_thenDoRealValueFunctionIteration

%% figures for Value Function Iteration with a Fixed Grid

figure
mesh(kk, aa, mValue');

title('Value - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Value','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_value_fixed_grid')

figure
mesh(kk, aa, mKPolicy');

title('Policy for Next Period Capital - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_kPolicy_fixed_grid')

figure
mesh(kk, aa, mLaborPolicy_1');

title('Policy for Good 1 Labor - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_1_fixed_grid')

figure
mesh(kk, aa, mLaborPolicy_2');

title('Policy for Good 2 Labor - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_2_fixed_grid')

figure
mesh(kk, aa, mConsumptionPolicy_1');

title('Policy for Good 1 Consumption - Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_1_fixed_grid')

figure
mesh(kk, aa, mConsumptionPolicy_2');

title('Policy for Good 2 Consumption- Fixed Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_2_fixed_grid')
