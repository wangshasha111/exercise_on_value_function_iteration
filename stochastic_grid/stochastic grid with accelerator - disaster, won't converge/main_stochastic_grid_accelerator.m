% Homework 1 Econ714 JFV
% A two good production economy
% Written by Shasha Wang
% Stochastic grid scheme 
% Reference: Rust 1997

% Implement a stochastic grid scheme (Rust, 1997) for a Value function iteration,
% with 500 vertex points with a coverage of +-25% of kss (you can keep the grid of investment fixed).
% Compare accuracy and computing time between the simple grid scheme implemented in 2) 
% and the results from the multigrid scheme. 
% Present evidence to support your claims.

%% Value Function Iteration with Stochastic grid and Accelerator

% Revised based on ..\main_accelerator.m
% Even though through experiment I found out that Stochastic grid with
% accelerator won't converge, I keep the algorithm for future tries. 
% To NOT use the accelerator, 
%                 if mod(iteration,10) >=0 
% To use the accelerator,
%                 if mod(iteration,10) ==1 

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
% NOTE: we also have to deal with COMPLEX solutions for fsolve when we use
% grid search

% ==============================================================

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

NkDraws = ceil(Nk/100);
seedKDraws=rng;
rng(seedKDraws);
ikStochasticGrid = randperm(Nk,NkDraws);%a ROW vector containing numberOfDraws unique integers selected randomly from 1 to Nk inclusive.
ikStochasticGrid = sort(ikStochasticGrid);

NaDraws = ceil(Na/2);
seedADraws=rng;
rng(seedADraws);
iaStochasticGrid = randperm(Na,NaDraws);%a ROW vector containing numberOfDraws unique integers selected randomly from 1 to Nk inclusive.
iaStochasticGrid = sort(iaStochasticGrid);

% Get the probability matrix for the random draws
mProb_a1a2StochasticGrid = zeros(Na,Na);
for iaPrime = iaStochasticGrid
     mProb_a1a2StochasticGrid(:,iaPrime) = mProb_a1a2(:,iaPrime);
end
% Then normalize it to a proper probability matrix
mProb_a1a2StochasticGrid = mProb_a1a2StochasticGrid./sum(mProb_a1a2StochasticGrid,2);

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
        
    expectedValue0 = mValue0 * mProb_a1a2StochasticGrid';% row: kPrime, column: a, i.e., expected value if today's shock is a and I choose k as tomorrow's capital stock

    if (mDifference(iteration) > tolerance)    
        for ia = 1:Na
            a_1 = mGrid_a1a2(ia,1);
            a_2 = mGrid_a1a2(ia,2);
            
            iikPrime = 1;

            for ik = 1:Nk
                k = vGrid_k(ik);
  
%                 if mod(iteration,10) ==1 % do maximization
                if mod(iteration,10) >=0 % do maximization
                    
                    valueHighSoFar = -1000.0;
                    kChoice  = vGrid_k(1);
            
                    for ikPrime = ikStochasticGrid(iikPrime:end)
                        kPrime = vGrid_k(ikPrime);
                        [valueProvisional,consumption_1,consumption_2] = valueFunction_stochasticGrid_setLaborToSteadyState...
                            (kPrime,ikPrime,ik,k,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL,labor_1_SteadyState,labor_2_SteadyState);

                        if (valueProvisional>valueHighSoFar)
                            valueHighSoFar = valueProvisional;
                            kChoice = vGrid_k(ikPrime);
                            iikPrime = sum(ikPrime >= ikStochasticGrid);
                        else
                            break
                        end
                    end
               
                    mValue(ik,ia) = valueHighSoFar;
                    mKPolicy(ik,ia) = kChoice;
                    mConsumptionPolicy_1(ik,ia) = consumption_1;
                    mConsumptionPolicy_2(ik,ia) = consumption_2;

                else % accelerator
                    currentUtility = utilityFunction(mConsumptionPolicy_1(ik,ia),mConsumptionPolicy_2(ik,ia),labor_1_SteadyState,labor_2_SteadyState,mmu_1,mmu_2);
                    expectedValue = expectedValue0(sum(mKPolicy(ik,ia)>=vGrid_k),ia);
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
            iikPrime = 1;

            for ik = 1:Nk
                k = vGrid_k(ik);
                
                valueHighSoFar = -1000.0;
                kChoice  = vGrid_k(1);

                for ikPrime = ikStochasticGrid(iikPrime:end)
                    kPrime = vGrid_k(ikPrime);
                    [valueProvisional,consumption_1,consumption_2] = valueFunction_stochasticGrid_setLaborToSteadyState...
                        (kPrime,ikPrime,ik,k,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL,labor_1_SteadyState,labor_2_SteadyState);

                    if (valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional;
                        kChoice = vGrid_k(ikPrime);
                        iikPrime = sum(ikPrime >= ikStochasticGrid);
                    else
                        break
                    end
                end

                mValue(ik,ia) = valueHighSoFar;
                mKPolicy(ik,ia) = kChoice;
                mConsumptionPolicy_1(ik,ia) = consumption_1;
                mConsumptionPolicy_2(ik,ia) = consumption_2;                

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

save ShashaWang_JFV_PS1_500_stochastic_capital_grid_points_valueFunctionIteration_setLaborToSteadyState

%  My computer, without accelerator
%  Iteration:  1, Sup diff: 0.00473112
%  Iteration: 11, Sup diff: 0.00036681
%  Iteration: 21, Sup diff: 0.00022816
%  Iteration: 31, Sup diff: 0.00014992
%  Iteration: 41, Sup diff: 0.00009946
%  Iteration: 51, Sup diff: 0.00006610
%  Iteration: 61, Sup diff: 0.00004394
%  Iteration: 71, Sup diff: 0.00002921
%  Iteration: 81, Sup diff: 0.00001942
%  Iteration: 91, Sup diff: 0.00001291
%  Iteration: 101, Sup diff: 0.00000858
%  Iteration: 111, Sup diff: 0.00000571
%  Iteration: 121, Sup diff: 0.00000379
%  Iteration: 131, Sup diff: 0.00000252
%  Iteration: 141, Sup diff: 0.00000168
%  Iteration: 151, Sup diff: 0.00000111
%  Iteration: 155, Sup diff: 0.00000095
% ÀúÊ± 24.143055 Ãë¡£
%  Convergence achieved. Total Number of Iteration: 155, Sup diff: 0.00000095

% Lab computer thanksgiving day 11/28/2019 4:34pm
%  Iteration:  1, Sup diff: 0.00442837
%  Iteration: 11, Sup diff: 0.00046245
%  Iteration: 21, Sup diff: 0.00030112
%  Iteration: 31, Sup diff: 0.00019410
%  Iteration: 41, Sup diff: 0.00012641
%  Iteration: 51, Sup diff: 0.00008234
%  Iteration: 61, Sup diff: 0.00005365
%  Iteration: 71, Sup diff: 0.00003496
%  Iteration: 81, Sup diff: 0.00002279
%  Iteration: 91, Sup diff: 0.00001485
%  Iteration: 101, Sup diff: 0.00000969
%  Iteration: 111, Sup diff: 0.00000632
%  Iteration: 121, Sup diff: 0.00000412
%  Iteration: 131, Sup diff: 0.00000269
%  Iteration: 141, Sup diff: 0.00000175
%  Iteration: 151, Sup diff: 0.00000114
%  Iteration: 156, Sup diff: 0.00000092
% Elapsed time is 13.969091 seconds.
%  Convergence achieved. Total Number of Iteration: 156, Sup diff: 0.00000092
 
% 5 draws of capital
%  Iteration:  1, Sup diff: 0.00417316
%  Iteration: 11, Sup diff: 0.00032735
%  Iteration: 21, Sup diff: 0.00018324
%  Iteration: 31, Sup diff: 0.00011859
%  Iteration: 41, Sup diff: 0.00007828
%  Iteration: 51, Sup diff: 0.00005167
%  Iteration: 61, Sup diff: 0.00003409
%  Iteration: 71, Sup diff: 0.00002250
%  Iteration: 81, Sup diff: 0.00001484
%  Iteration: 91, Sup diff: 0.00000979
%  Iteration: 101, Sup diff: 0.00000646
%  Iteration: 111, Sup diff: 0.00000426
%  Iteration: 121, Sup diff: 0.00000281
%  Iteration: 131, Sup diff: 0.00000185
%  Iteration: 141, Sup diff: 0.00000122
%  Iteration: 147, Sup diff: 0.00000095
% Elapsed time is 12.511260 seconds.
%  Convergence achieved. Total Number of Iteration: 147, Sup diff: 0.00000095

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
        
    expectedValue0 = mValue0 * mProb_a1a2StochasticGrid';% row: kPrime, column: a, i.e., expected value if today's shock is a and I choose k as tomorrow's capital stock
    
    if (mDifference(iteration) > tolerance)    

        for ia = 1:Na
            a_1 = mGrid_a1a2(ia,1);
            a_2 = mGrid_a1a2(ia,2);

            iikPrime = 1;
            
            for ik = 1:Nk
                k = vGrid_k(ik);
                
%                 if mod(iteration,10) == 1 % do maximization
                if mod(iteration,10) >= 0 % do maximization
                    valueHighSoFar = -1000.0;
                    kChoice  = vGrid_k(1);
%                     laborInitial=[labor_1_SteadyState,labor_2_SteadyState];
                    
                    for ikPrime = ikStochasticGrid(iikPrime:end)
                        kPrime = vGrid_k(ikPrime);

                        [valueProvisional,labor_1,labor_2,consumption_1,consumption_2] = valueFunction_stochasticGrid...
                            (kPrime,ikPrime,ik,k,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL,laborInitial);

                        if (valueProvisional>valueHighSoFar)
                            valueHighSoFar = valueProvisional;
                            kChoice = vGrid_k(ikPrime);
                            iikPrime = sum(ikPrime >= ikStochasticGrid);
%                             laborInitial = [labor_1,labor_2];
                            mValue(ik,ia) = valueHighSoFar;
                            mKPolicy(ik,ia) = kChoice;
                            mLaborPolicy_1(ik,ia) = labor_1;
                            mLaborPolicy_2(ik,ia) = labor_2;
                            mConsumptionPolicy_1(ik,ia) = consumption_1;
                            mConsumptionPolicy_2(ik,ia) = consumption_2;
                        else
                            break
                        end
                    end
               

 
                else
                    currentUtility = utilityFunction(mConsumptionPolicy_1(ik,ia),mConsumptionPolicy_2(ik,ia),mLaborPolicy_1(ik,ia),mLaborPolicy_2(ik,ia),mmu_1,mmu_2);
                    expectedValue = expectedValue0(sum(mKPolicy(ik,ia)>=vGrid_k),ia);
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

            iikPrime = 1;
            
            for ik = 1:Nk
                k = vGrid_k(ik);
                
%                 if mod(iteration,10) == 1
                    valueHighSoFar = -1000.0;
                    kChoice  = vGrid_k(1);
%                     laborInitial=[labor_1_SteadyState,labor_2_SteadyState];
            
                    for ikPrime = ikStochasticGrid(iikPrime:end)
                        kPrime = vGrid_k(ikPrime);
                        

                        [valueProvisional,labor_1,labor_2,consumption_1,consumption_2] = valueFunction_stochasticGrid...
                            (kPrime,ikPrime,ik,k,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL,laborInitial);
                        

                        if (valueProvisional>valueHighSoFar)
                            valueHighSoFar = valueProvisional;
                            kChoice = vGrid_k(ikPrime);
                            iikPrime = sum(ikPrime >= ikStochasticGrid);
%                             laborInitial = [labor_1,labor_2];
                            mValue(ik,ia) = valueHighSoFar;
                            mKPolicy(ik,ia) = kChoice;
                            mLaborPolicy_1(ik,ia) = labor_1;
                            mLaborPolicy_2(ik,ia) = labor_2;
                            mConsumptionPolicy_1(ik,ia) = consumption_1;
                            mConsumptionPolicy_2(ik,ia) = consumption_2;
                        else
                            break
                        end
                    end
               
                    mValue(ik,ia) = valueHighSoFar;
                    mKPolicy(ik,ia) = kChoice;
                    mLaborPolicy_1(ik,ia) = labor_1;
                    mLaborPolicy_2(ik,ia) = labor_2;
                    mConsumptionPolicy_1(ik,ia) = consumption_1;
                    mConsumptionPolicy_2(ik,ia) = consumption_2;

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


end

toc

fprintf(' Convergence achieved. Total Number of Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 

save ShashaWang_JFV_PS1_500_stochastic_capital_grid_points_valueFunctionIteration_setLaborToSteadyState_thenDoRealValueFunctionIteration_accelerator


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


% Lab computer Thanksgiving day 11.28.2019 8:53pm
%  Iteration:  1, Sup diff: 0.00156838
%  Iteration:  2, Sup diff: 0.00119851
%  Iteration:  3, Sup diff: 0.00082321
%  Iteration:  4, Sup diff: 0.00078964
%  Iteration:  5, Sup diff: 0.00075744
%  Iteration:  6, Sup diff: 0.00072656
%  Iteration:  7, Sup diff: 0.00069694
%  Iteration:  8, Sup diff: 0.00066854
%  Iteration:  9, Sup diff: 0.00064129
%  Iteration: 10, Sup diff: 0.00061515
%  Iteration: 11, Sup diff: 0.00059008
%  Iteration: 12, Sup diff: 0.00056604
%  Iteration: 13, Sup diff: 0.00054297
%  Iteration: 14, Sup diff: 0.00052085
%  Iteration: 15, Sup diff: 0.00049962
%  Iteration: 16, Sup diff: 0.00047927
%  Iteration: 17, Sup diff: 0.00045974
%  Iteration: 18, Sup diff: 0.00044101
%  Iteration: 19, Sup diff: 0.00042304
%  Iteration: 20, Sup diff: 0.00040581
%  Iteration: 21, Sup diff: 0.00038928
%  Iteration: 22, Sup diff: 0.00037342
%  Iteration: 23, Sup diff: 0.00035821
%  Iteration: 24, Sup diff: 0.00034361
%  Iteration: 25, Sup diff: 0.00032962
%  Iteration: 26, Sup diff: 0.00031619
%  Iteration: 27, Sup diff: 0.00030331
%  Iteration: 28, Sup diff: 0.00029096
%  Iteration: 29, Sup diff: 0.00027911
%  Iteration: 30, Sup diff: 0.00026774
%  Iteration: 31, Sup diff: 0.00025684
%  Iteration: 32, Sup diff: 0.00024638
%  Iteration: 33, Sup diff: 0.00023635
%  Iteration: 34, Sup diff: 0.00022672
%  Iteration: 35, Sup diff: 0.00021749
%  Iteration: 36, Sup diff: 0.00020864
%  Iteration: 37, Sup diff: 0.00020014
%  Iteration: 38, Sup diff: 0.00019200
%  Iteration: 39, Sup diff: 0.00018418
%  Iteration: 40, Sup diff: 0.00017668
%  Iteration: 41, Sup diff: 0.00016949
%  Iteration: 42, Sup diff: 0.00016259
%  Iteration: 43, Sup diff: 0.00015598
%  Iteration: 44, Sup diff: 0.00014963
%  Iteration: 45, Sup diff: 0.00014354
%  Iteration: 46, Sup diff: 0.00013770
%  Iteration: 47, Sup diff: 0.00013210
%  Iteration: 48, Sup diff: 0.00012672
%  Iteration: 49, Sup diff: 0.00012157
%  Iteration: 50, Sup diff: 0.00011662
%  Iteration: 51, Sup diff: 0.00011188
%  Iteration: 52, Sup diff: 0.00010733
%  Iteration: 53, Sup diff: 0.00010296
%  Iteration: 54, Sup diff: 0.00009877
%  Iteration: 55, Sup diff: 0.00009476
%  Iteration: 56, Sup diff: 0.00009090
%  Iteration: 57, Sup diff: 0.00008721
%  Iteration: 58, Sup diff: 0.00008366
%  Iteration: 59, Sup diff: 0.00008026
%  Iteration: 60, Sup diff: 0.00007700
%  Iteration: 61, Sup diff: 0.00007387
%  Iteration: 62, Sup diff: 0.00007086
%  Iteration: 63, Sup diff: 0.00006798
%  Iteration: 64, Sup diff: 0.00006522
%  Iteration: 65, Sup diff: 0.00006257
%  Iteration: 66, Sup diff: 0.00006003
%  Iteration: 67, Sup diff: 0.00005759
%  Iteration: 68, Sup diff: 0.00005525
%  Iteration: 69, Sup diff: 0.00005300
%  Iteration: 70, Sup diff: 0.00005085
%  Iteration: 71, Sup diff: 0.00004879
%  Iteration: 72, Sup diff: 0.00004680
%  Iteration: 73, Sup diff: 0.00004490
%  Iteration: 74, Sup diff: 0.00004308
%  Iteration: 75, Sup diff: 0.00004133
%  Iteration: 76, Sup diff: 0.00003965
%  Iteration: 77, Sup diff: 0.00003804
%  Iteration: 78, Sup diff: 0.00003650
%  Iteration: 79, Sup diff: 0.00003502
%  Iteration: 80, Sup diff: 0.00003359
%  Iteration: 81, Sup diff: 0.00003223
%  Iteration: 82, Sup diff: 0.00003092
%  Iteration: 83, Sup diff: 0.00002967
%  Iteration: 84, Sup diff: 0.00002846
%  Iteration: 85, Sup diff: 0.00002731
%  Iteration: 86, Sup diff: 0.00002620
%  Iteration: 87, Sup diff: 0.00002514
%  Iteration: 88, Sup diff: 0.00002412
%  Iteration: 89, Sup diff: 0.00002314
%  Iteration: 90, Sup diff: 0.00002220
%  Iteration: 91, Sup diff: 0.00002130
%  Iteration: 92, Sup diff: 0.00002044
%  Iteration: 93, Sup diff: 0.00001961
%  Iteration: 94, Sup diff: 0.00001881
%  Iteration: 95, Sup diff: 0.00001805
%  Iteration: 96, Sup diff: 0.00001732
%  Iteration: 97, Sup diff: 0.00001662
%  Iteration: 98, Sup diff: 0.00001594
%  Iteration: 99, Sup diff: 0.00001530
%  Iteration: 100, Sup diff: 0.00001468
%  Iteration: 101, Sup diff: 0.00001408
%  Iteration: 102, Sup diff: 0.00001351
%  Iteration: 103, Sup diff: 0.00001296
%  Iteration: 104, Sup diff: 0.00001244
%  Iteration: 105, Sup diff: 0.00001193
%  Iteration: 106, Sup diff: 0.00001145
%  Iteration: 107, Sup diff: 0.00001099
%  Iteration: 108, Sup diff: 0.00001054
%  Iteration: 109, Sup diff: 0.00001011
%  Iteration: 110, Sup diff: 0.00000971
%  Iteration: 111, Sup diff: 0.00000931
%  Iteration: 112, Sup diff: 0.00000893
%  Iteration: 113, Sup diff: 0.00000857
%  Iteration: 114, Sup diff: 0.00000823
%  Iteration: 115, Sup diff: 0.00000789
%  Iteration: 116, Sup diff: 0.00000757
%  Iteration: 117, Sup diff: 0.00000727
%  Iteration: 118, Sup diff: 0.00000697
%  Iteration: 119, Sup diff: 0.00000669
%  Iteration: 120, Sup diff: 0.00000642
%  Iteration: 121, Sup diff: 0.00000616
%  Iteration: 122, Sup diff: 0.00000591
%  Iteration: 123, Sup diff: 0.00000567
%  Iteration: 124, Sup diff: 0.00000544
%  Iteration: 125, Sup diff: 0.00000522
%  Iteration: 126, Sup diff: 0.00000501
%  Iteration: 127, Sup diff: 0.00000481
%  Iteration: 128, Sup diff: 0.00000461
%  Iteration: 129, Sup diff: 0.00000443
%  Iteration: 130, Sup diff: 0.00000425
%  Iteration: 131, Sup diff: 0.00000408
%  Iteration: 132, Sup diff: 0.00000391
%  Iteration: 133, Sup diff: 0.00000375
%  Iteration: 134, Sup diff: 0.00000360
%  Iteration: 135, Sup diff: 0.00000346
%  Iteration: 136, Sup diff: 0.00000332
%  Iteration: 137, Sup diff: 0.00000318
%  Iteration: 138, Sup diff: 0.00000305
%  Iteration: 139, Sup diff: 0.00000293
%  Iteration: 140, Sup diff: 0.00000281
%  Iteration: 141, Sup diff: 0.00000270
%  Iteration: 142, Sup diff: 0.00000259
%  Iteration: 143, Sup diff: 0.00000248
%  Iteration: 144, Sup diff: 0.00000238
%  Iteration: 145, Sup diff: 0.00000229
%  Iteration: 146, Sup diff: 0.00000220
%  Iteration: 147, Sup diff: 0.00000211
%  Iteration: 148, Sup diff: 0.00000202
%  Iteration: 149, Sup diff: 0.00000194
%  Iteration: 150, Sup diff: 0.00000186
%  Iteration: 151, Sup diff: 0.00000179
%  Iteration: 152, Sup diff: 0.00000171
%  Iteration: 153, Sup diff: 0.00000164
%  Iteration: 154, Sup diff: 0.00000158
%  Iteration: 155, Sup diff: 0.00000151
%  Iteration: 156, Sup diff: 0.00000145
%  Iteration: 157, Sup diff: 0.00000139
%  Iteration: 158, Sup diff: 0.00000134
%  Iteration: 159, Sup diff: 0.00000128
%  Iteration: 160, Sup diff: 0.00000123
%  Iteration: 161, Sup diff: 0.00000118
%  Iteration: 162, Sup diff: 0.00000114
%  Iteration: 163, Sup diff: 0.00000109
%  Iteration: 164, Sup diff: 0.00000105
%  Iteration: 165, Sup diff: 0.00000100
%  Iteration: 166, Sup diff: 0.00000096
%  Iteration: 167, Sup diff: 0.00000092
% Elapsed time is 10393.859303 seconds.
%  Convergence achieved. Total Number of Iteration: 167, Sup diff: 0.00000092

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
