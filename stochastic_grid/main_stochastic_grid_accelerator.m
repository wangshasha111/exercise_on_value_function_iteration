% Homework 1 Econ714 JFV
% A two good production economy
% Written by Shasha Wang
% Reference: Rust 1997

%% Value Function Iteration with Stochastic grid and Accelerator

% Revised based on ..\main_accelerator.m

% Implement a stochastic grid scheme (Rust, 1997) for a Value function iteration,
% with 500 vertex points with a coverage of +-25% of kss (you can keep the grid of investment fixed).
% Compare accuracy and computing time between the simple grid scheme implemented in 2) 
% and the results from the multigrid scheme. 
% Present evidence to support your claims.

%% Algorithm explained
% The Random Bellman operator (RBO) is a mapping 
% where next period states are N states randomly chosen from all states.
% Keep the points chosen fixed for iterations till convergence.

% Note when computing the expected value, we need to use
% CONDITIONAL PROBABILITY - conditional on that the N states are chosen, i.e., we need
% to divide each transition probability by the sum of the transition probability of N chosen states

% ==============================================================

% The way to draw random numbers

% FIRST step: set up the seed

% stochastic_grid=rng;
% rng(stochastic_grid)

% SECOND step: use randperm to generate unique integers. Don't use randi,
% since numbers won't be unique.
% p = randperm(n,k) returns a ROW vector containing k unique integers selected randomly from 1 to n inclusive.
% X = randi(imax) returns a pseudorandom scalar integer between 1 and imax.
% X = randi(imax,sz1,...,szN) returns an sz1-by-...-by-szN array where sz1,...,szN indicates the size of each dimension. For example, randi(10,3,4) returns a 3-by-4 array of pseudorandom integers between 1 and 10.

% ==============================================================
% Grid Search VS fminsearch
% In the value function iteration, I used grid search instead of linear
% interpolation, because I can't randomize tomorrow's state using
% fminsearch/con/bnd. 
% Because we are doing grid search, CONCAVITY of value function and
% MONOTONICITY of policy function can be utilized to speed up the process.

main_setup % parameters, shocks

%% Value function iteration

% As in the fixed-grid problem, we also first set labor level to steady
% state level and compute a converged value function, then use it as the
% initial guess for the real and qualified value function iteration. 

%% Step 1: Compute the Steady State

% Use fsolve to solve the system of equations
    
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

kSpread = 0.25;
kMax = kSteadyState * (1 + kSpread);
kMin = kSteadyState * (1 - kSpread);
Nk = 500;
vGrid_k = linspace(kMin,kMax,Nk)';
inputs.vGrid_k = vGrid_k;

%% Step 3. Pre-Iteration: Set labor to steady state

%% Draw random grid points

% Most generally, we should draw random draws from all possible states, i.e., CAPITAL * SHOCKS.
% If we have 500 points for capital, 15 points for shocks, then altogether
% we have 7500 states, from which we draw.

% But we are told in the homework that the grid for capital should be fixed. 
% And we only have two dimensions of shocks. 
% That means the other dimension also needs to be fixed. 
% So we are going to take draws from different level of capital as well as from different shocks.

NkDraws = ceil(Nk/5);
seedKDraws=rng;
rng(seedKDraws);
ikStochasticGrid = randperm(Nk,NkDraws);%a ROW vector containing numberOfDraws unique integers selected randomly from 1 to Nk inclusive.
ikStochasticGrid = sort(ikStochasticGrid);

NaDraws = ceil(Na/2);
seedADraws=rng;
rng(seedADraws);
iaStochasticGrid = randperm(Na,NaDraws);%a ROW vector containing numberOfDraws unique integers selected randomly from 1 to Nk inclusive.
iaStochasticGrid = sort(iaStochasticGrid);





mProb_a1a2StochasticGrid = mProb_a1a2;
for iaPrime = iaStochasticGrid
     mProb_a1a2StochasticGrid(:,iaPrime) = zeros(Na,1);
end
mProb_a1a2StochasticGrid = mProb_a1a2StochasticGrid

vNormalizer = sum(mProb_a1a2(:,iaStochasticGrid),2);% 15 by 1 row vector. Use it to normalize expected value










%% Required matrices and vectors for Pre-Iteration

mValue0                 = utilitySteadyState.*ones(Nk,Na); inputs.mValue0 = mValue0;
mValue                  = zeros(Nk,Na);
mKPolicy                = zeros(Nk,Na);
% mLaborPolicy_1          = zeros(Nk,Na);
% mLaborPolicy_2          = zeros(Nk,Na); 
mConsumptionPolicy_1    = zeros(Nk,Na);
mConsumptionPolicy_2    = zeros(Nk,Na);

maxIter = 10000;
tolerance = 1e-6;

iteration = 1;
mDifference = zeros(maxIter,1);
mDifference(iteration) = 100;

options = optimset('Display', 'off');
% opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','off');
% laborInitial=[labor_1_SteadyState,labor_2_SteadyState];

tic

% Finding value and policy functions numerically
while iteration <= maxIter  ...% make sure the last iteration does the maximization
        && ( (mDifference(iteration) > tolerance)  | (mDifference(iteration) <= tolerance && mod(iteration,10)~=2)) 
        
    expectedValue0 = mValue0 * mProb_a1a2';% row: kPrime, column: a, i.e., expected value if today's shock is a and I choose k as tomorrow's capital stock

    if (mDifference(iteration) > tolerance)    
        for ia = 1:Na
            a_1 = mGrid_a1a2(ia,1);
            a_2 = mGrid_a1a2(ia,2);
            gridCapitalNextPeriod = 1;
            normalizedExpectedValue0 = expectedValue0(:,ia)/vNormalizer(ia); % Nk by 1 colum vector

            for ik = 1:Nk
                k = vGrid_k(ik);
  
                if mod(iteration,10) ==1 % do maximization
                    
                    valueHighSoFar = -1000.0;
                    kChoice  = vGrid_k(1);
            
                    for ikPrime = iaStochasticGrid(gridCapitalNextPeriod:end)
                        kPrime = vGrid_k(ikPrime);
                        outputValueFunction = valueFunction_stochasticGrid_setLaborToSteadyState(kPrime,ikPrime,ik,k,ia,a_1,a_2,normalizedExpectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL,labor_1_SteadyState,labor_2_SteadyState);
 
                        valueProvisional = outputValueFunction(1);
                        consumption_1 = outputValueFunction(2);
                        consumption_2 = outputValueFunction(3);

                        
                        
                        
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
                    mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,labor_1_SteadyState,aalphaK,aalphaL,ddelta);
                    mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,labor_2_SteadyState);

                else
                    currentUtility = utilityFunction(mConsumptionPolicy_1(ik,ia),mConsumptionPolicy_2(ik,ia),labor_1_SteadyState,labor_2_SteadyState,mmu_1,mmu_2);;
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

                mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,labor_1_SteadyState,aalphaK,aalphaL,ddelta);
                mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,labor_2_SteadyState);
            end
        end
        iteration = iteration + 1;
        mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
        mValue0         = mValue;

        fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 
        
        if mDifference(iteration) <= tolerance
            break
        end
    end

end

toc

fprintf(' Convergence achieved. Total Number of Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 

save ShashaWang_JFV_PS1_250_capital_grid_points_valueFunctionIteration_setLaborToSteadyState_accelerator
% lab computer
%  Iteration:  1, Sup diff: 0.00516861
%  Iteration: 11, Sup diff: 0.00067347
%  Iteration: 21, Sup diff: 0.00018309
%  Iteration: 31, Sup diff: 0.00011045
%  Iteration: 41, Sup diff: 0.00010236
%  Iteration: 51, Sup diff: 0.00006380
%  Iteration: 61, Sup diff: 0.00001971
%  Iteration: 71, Sup diff: 0.00001192
%  Iteration: 81, Sup diff: 0.00000730
%  Iteration: 91, Sup diff: 0.00000444
%  Iteration: 101, Sup diff: 0.00000273
%  Iteration: 111, Sup diff: 0.00000169
%  Iteration: 121, Sup diff: 0.00000105
%  Iteration: 123, Sup diff: 0.00000095
% Elapsed time is 84.668544 seconds.
%  Convergence achieved. Total Number of Iteration: 123, Sup diff: 0.00000095


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
while iteration <= maxIter  ...% make sure the last iteration does the maximization
        && ( (mDifference(iteration) > tolerance)  | (mDifference(iteration) <= tolerance && mod(iteration,10)~=2)) 
        
    expectedValue0 = mValue0 * mProb_a1a2';
    
    if (mDifference(iteration) > tolerance)    

        for ia = 1:Na
            a_1 = mGrid_a1a2(ia,1);
            a_2 = mGrid_a1a2(ia,2);

            for ik = 1:Nk
                k = vGrid_k(ik);
                
                if mod(iteration,10) == 1

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

                    vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
                    mLaborPolicy_1(ik,ia) = vLabor(1);
                    mLaborPolicy_2(ik,ia) = vLabor(2);
                    mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
                    mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,vLabor(2));
                    laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
                else
%                     currentUtility = interp1(vGrid_k,mCurrentUtilityFsolve(:,ia,ik),mKPolicy(ik,ia));
                    currentUtility = utilityFunction(mConsumptionPolicy_1(ik,ia),mConsumptionPolicy_2(ik,ia),mLaborPolicy_1(ik,ia),mLaborPolicy_2(ik,ia),mmu_1,mmu_2);;
                    expectedValue = interp1(vGrid_k,expectedValue0(:,ia),mKPolicy(ik,ia));
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

    else
        for ia = 1:Na
            a_1 = mGrid_a1a2(ia,1);
            a_2 = mGrid_a1a2(ia,2);

            for ik = 1:Nk
                k = vGrid_k(ik);
                
%                 if mod(iteration,10) == 1

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

                vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
                mLaborPolicy_1(ik,ia) = vLabor(1);
                mLaborPolicy_2(ik,ia) = vLabor(2);
                mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
                mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,vLabor(2));
                laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
%                 else
% %                     currentUtility = interp1(vGrid_k,mCurrentUtilityFsolve(:,ia,ik),mKPolicy(ik,ia));
%                     currentUtility = utilityFunction(mConsumptionPolicy_1(ik,ia),mConsumptionPolicy_2(ik,ia),mLaborPolicy_1(ik,ia),mLaborPolicy_2(ik,ia),mmu_1,mmu_2);;
%                     expectedValue = interp1(vGrid_k,expectedValue0(:,ia),mKPolicy(ik,ia));
%                     value = (1-bbeta)*currentUtility + bbeta * expectedValue;
%                     
%                     mValue(ik,ia) = value;
%                     
%                 end
            end
        end
        iteration = iteration + 1;
        mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
        mValue0         = mValue;

%         if mod(iteration,10) == 2
        fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 
%         end        

        if mDifference(iteration) <= tolerance
            break
        end
    end
%         iteration = iteration + 1;
%         mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
%         mValue0         = mValue;
% 
%         fprintf(' Iteration: %2.0f, Sup diff: %2.6f\n', iteration-1, mDifference(iteration)); 

end

toc

fprintf(' Convergence achieved. Total Number of Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 




%% For accuracy test, compute the euler equation error

% Let's use LINEAR INTERPOLATION for tomorrow's values

errorEulerEquationLinearInterpolation = eulerEquationErrorFunction(Nk,vGrid_k,mKPolicy,mLaborPolicy_1,mConsumptionPolicy_1,mConsumptionPolicy_2,Na,mGrid_a1a2,mProb_a1a2,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL);
errorEulerEquationLinearInterpolationDecimalLog = log10( errorEulerEquationLinearInterpolation );

[kk,aa]=meshgrid(vGrid_k, mGrid_a1a2(:,1));

figure
mesh(kk, aa, errorEulerEquationLinearInterpolationDecimalLog');

title('Euler Equation Error $log_{10}$ Linear Interpolation - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('error','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_eulerEquationErrorLinearInterpolation_accelerator')

save ShashaWang_JFV_PS1_250_capital_grid_points_valueFunctionIteration_setLaborToSteadyState_thenDoRealValueFunctionIteration_accelerator

%  Iteration:  1, Sup diff: 0.00128356
%  Iteration: 11, Sup diff: 0.00061426
%  Iteration: 21, Sup diff: 0.00038861
%  Iteration: 31, Sup diff: 0.00025225
%  Iteration: 41, Sup diff: 0.00016538
%  Iteration: 51, Sup diff: 0.00010888
%  Iteration: 61, Sup diff: 0.00007183
%  Iteration: 71, Sup diff: 0.00004744
%  Iteration: 81, Sup diff: 0.00003137
%  Iteration: 91, Sup diff: 0.00002075
%  Iteration: 101, Sup diff: 0.00001374
%  Iteration: 111, Sup diff: 0.00000910
%  Iteration: 121, Sup diff: 0.00000603
%  Iteration: 131, Sup diff: 0.00000400
%  Iteration: 141, Sup diff: 0.00000265
%  Iteration: 151, Sup diff: 0.00000176
%  Iteration: 161, Sup diff: 0.00000117
% Elapsed time is 2100.471542 seconds.
%  Convergence achieved. Total Number of Iteration: 166, Sup diff: 0.00000095

%% figures for Value Function Iteration with a Fixed Grid

figure
mesh(kk, aa, mValue');

title('Value - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Value','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_value_accelerator')

figure
mesh(kk, aa, mKPolicy');

title('Policy for Next Period Capital - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_kPolicy_accelerator')

figure
mesh(kk, aa, mLaborPolicy_1');

title('Policy for Good 1 Labor - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_1_accelerator')

figure
mesh(kk, aa, mLaborPolicy_2');

title('Policy for Good 2 Labor - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_2_accelerator')

figure
mesh(kk, aa, mConsumptionPolicy_1');

title('Policy for Good 1 Consumption - Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_1_accelerator')

figure
mesh(kk, aa, mConsumptionPolicy_2');

title('Policy for Good 2 Consumption- Accelerator','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_2_accelerator')
