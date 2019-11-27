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

% Iterate on the Value function implied by the Social Planner’s Problem using linear interpolation
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
% 历时 829.203325 秒。

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
%  Iteration:  1, Sup diff: 0.001283
%  Iteration:  2, Sup diff: 0.001218
%  Iteration:  3, Sup diff: 0.001157
%  Iteration:  4, Sup diff: 0.000970
%  Iteration:  5, Sup diff: 0.000831
%  Iteration:  6, Sup diff: 0.000788
%  Iteration:  7, Sup diff: 0.000748
%  Iteration:  8, Sup diff: 0.000712
%  Iteration:  9, Sup diff: 0.000677
%  Iteration: 10, Sup diff: 0.000644
%  Iteration: 11, Sup diff: 0.000614
%  Iteration: 12, Sup diff: 0.000585
%  Iteration: 13, Sup diff: 0.000558
%  Iteration: 14, Sup diff: 0.000533
%  Iteration: 15, Sup diff: 0.000509
%  Iteration: 16, Sup diff: 0.000486
%  Iteration: 17, Sup diff: 0.000464
%  Iteration: 18, Sup diff: 0.000444
%  Iteration: 19, Sup diff: 0.000424
%  Iteration: 20, Sup diff: 0.000406
%  Iteration: 21, Sup diff: 0.000388
%  Iteration: 22, Sup diff: 0.000372
%  Iteration: 23, Sup diff: 0.000356
%  Iteration: 24, Sup diff: 0.000341
%  Iteration: 25, Sup diff: 0.000326
%  Iteration: 26, Sup diff: 0.000312
%  Iteration: 27, Sup diff: 0.000299
%  Iteration: 28, Sup diff: 0.000287
%  Iteration: 29, Sup diff: 0.000275
%  Iteration: 30, Sup diff: 0.000263
%  Iteration: 31, Sup diff: 0.000252
%  Iteration: 32, Sup diff: 0.000242
%  Iteration: 33, Sup diff: 0.000232
%  Iteration: 34, Sup diff: 0.000222
%  Iteration: 35, Sup diff: 0.000213
%  Iteration: 36, Sup diff: 0.000204
%  Iteration: 37, Sup diff: 0.000196
%  Iteration: 38, Sup diff: 0.000188
%  Iteration: 39, Sup diff: 0.000180
%  Iteration: 40, Sup diff: 0.000172
%  Iteration: 41, Sup diff: 0.000165
%  Iteration: 42, Sup diff: 0.000159
%  Iteration: 43, Sup diff: 0.000152
%  Iteration: 44, Sup diff: 0.000146
%  Iteration: 45, Sup diff: 0.000140
%  Iteration: 46, Sup diff: 0.000134
%  Iteration: 47, Sup diff: 0.000129
%  Iteration: 48, Sup diff: 0.000123
%  Iteration: 49, Sup diff: 0.000118
%  Iteration: 50, Sup diff: 0.000113
%  Iteration: 51, Sup diff: 0.000109
%  Iteration: 52, Sup diff: 0.000104
%  Iteration: 53, Sup diff: 0.000100
%  Iteration: 54, Sup diff: 0.000096
%  Iteration: 55, Sup diff: 0.000092
%  Iteration: 56, Sup diff: 0.000088
%  Iteration: 57, Sup diff: 0.000085
%  Iteration: 58, Sup diff: 0.000081
%  Iteration: 59, Sup diff: 0.000078
%  Iteration: 60, Sup diff: 0.000075
%  Iteration: 61, Sup diff: 0.000072
%  Iteration: 62, Sup diff: 0.000069
%  Iteration: 63, Sup diff: 0.000066
%  Iteration: 64, Sup diff: 0.000063
%  Iteration: 65, Sup diff: 0.000061
%  Iteration: 66, Sup diff: 0.000058
%  Iteration: 67, Sup diff: 0.000056
%  Iteration: 68, Sup diff: 0.000054
%  Iteration: 69, Sup diff: 0.000052
%  Iteration: 70, Sup diff: 0.000049
%  Iteration: 71, Sup diff: 0.000047
%  Iteration: 72, Sup diff: 0.000046
%  Iteration: 73, Sup diff: 0.000044
%  Iteration: 74, Sup diff: 0.000042
%  Iteration: 75, Sup diff: 0.000040
%  Iteration: 76, Sup diff: 0.000039
%  Iteration: 77, Sup diff: 0.000037
%  Iteration: 78, Sup diff: 0.000035
%  Iteration: 79, Sup diff: 0.000034
%  Iteration: 80, Sup diff: 0.000033
%  Iteration: 81, Sup diff: 0.000031
%  Iteration: 82, Sup diff: 0.000030
%  Iteration: 83, Sup diff: 0.000029
%  Iteration: 84, Sup diff: 0.000028
%  Iteration: 85, Sup diff: 0.000027
%  Iteration: 86, Sup diff: 0.000026
%  Iteration: 87, Sup diff: 0.000024
%  Iteration: 88, Sup diff: 0.000023
%  Iteration: 89, Sup diff: 0.000023
%  Iteration: 90, Sup diff: 0.000022
%  Iteration: 91, Sup diff: 0.000021
%  Iteration: 92, Sup diff: 0.000020
%  Iteration: 93, Sup diff: 0.000019
%  Iteration: 94, Sup diff: 0.000018
%  Iteration: 95, Sup diff: 0.000018
%  Iteration: 96, Sup diff: 0.000017
%  Iteration: 97, Sup diff: 0.000016
%  Iteration: 98, Sup diff: 0.000016
%  Iteration: 99, Sup diff: 0.000015
%  Iteration: 100, Sup diff: 0.000014
%  Iteration: 101, Sup diff: 0.000014
%  Iteration: 102, Sup diff: 0.000013
%  Iteration: 103, Sup diff: 0.000013
%  Iteration: 104, Sup diff: 0.000012
%  Iteration: 105, Sup diff: 0.000012
%  Iteration: 106, Sup diff: 0.000011
%  Iteration: 107, Sup diff: 0.000011
%  Iteration: 108, Sup diff: 0.000010
%  Iteration: 109, Sup diff: 0.000010
%  Iteration: 110, Sup diff: 0.000009
%  Iteration: 111, Sup diff: 0.000009
%  Iteration: 112, Sup diff: 0.000009
%  Iteration: 113, Sup diff: 0.000008
%  Iteration: 114, Sup diff: 0.000008
%  Iteration: 115, Sup diff: 0.000008
%  Iteration: 116, Sup diff: 0.000007
%  Iteration: 117, Sup diff: 0.000007
%  Iteration: 118, Sup diff: 0.000007
%  Iteration: 119, Sup diff: 0.000007
%  Iteration: 120, Sup diff: 0.000006
%  Iteration: 121, Sup diff: 0.000006
%  Iteration: 122, Sup diff: 0.000006
%  Iteration: 123, Sup diff: 0.000006
%  Iteration: 124, Sup diff: 0.000005
%  Iteration: 125, Sup diff: 0.000005
%  Iteration: 126, Sup diff: 0.000005
%  Iteration: 127, Sup diff: 0.000005
%  Iteration: 128, Sup diff: 0.000005
%  Iteration: 129, Sup diff: 0.000004
%  Iteration: 130, Sup diff: 0.000004
%  Iteration: 131, Sup diff: 0.000004
%  Iteration: 132, Sup diff: 0.000004
%  Iteration: 133, Sup diff: 0.000004
%  Iteration: 134, Sup diff: 0.000004
%  Iteration: 135, Sup diff: 0.000003
%  Iteration: 136, Sup diff: 0.000003
%  Iteration: 137, Sup diff: 0.000003
%  Iteration: 138, Sup diff: 0.000003
%  Iteration: 139, Sup diff: 0.000003
%  Iteration: 140, Sup diff: 0.000003
%  Iteration: 141, Sup diff: 0.000003
%  Iteration: 142, Sup diff: 0.000003
%  Iteration: 143, Sup diff: 0.000002
%  Iteration: 144, Sup diff: 0.000002
%  Iteration: 145, Sup diff: 0.000002
%  Iteration: 146, Sup diff: 0.000002
%  Iteration: 147, Sup diff: 0.000002
%  Iteration: 148, Sup diff: 0.000002
%  Iteration: 149, Sup diff: 0.000002
%  Iteration: 150, Sup diff: 0.000002
%  Iteration: 151, Sup diff: 0.000002
%  Iteration: 152, Sup diff: 0.000002
%  Iteration: 153, Sup diff: 0.000002
%  Iteration: 154, Sup diff: 0.000002
%  Iteration: 155, Sup diff: 0.000001
%  Iteration: 156, Sup diff: 0.000001
%  Iteration: 157, Sup diff: 0.000001
%  Iteration: 158, Sup diff: 0.000001
%  Iteration: 159, Sup diff: 0.000001
%  Iteration: 160, Sup diff: 0.000001
%  Iteration: 161, Sup diff: 0.000001
%  Iteration: 162, Sup diff: 0.000001
%  Iteration: 163, Sup diff: 0.000001
%  Iteration: 164, Sup diff: 0.000001
%  Iteration: 165, Sup diff: 0.000001
% 历时 34715.242191 秒。
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
