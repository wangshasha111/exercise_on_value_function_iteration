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

%% Value Function Iteration with Stochastic grid

% Revised based on ..\main_accelerator.m
% Even though through experiment I found out that Stochastic grid with
% accelerator won't converge, I keep the algorithm for future tries. 
% To NOT use the accelerator, 
%                 if mod(iteration,10) >=0 
% To use the accelerator,
%                 if mod(iteration,10) ==1 

%% Algorithm explained
% The Random Bellman Operator (RBO) is a mapping 
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

% SECOND step: use randperm to generate unique integers. Don't use randi since numbers won't be unique.
% p = randperm(n,k) returns a ROW vector containing k unique integers selected randomly from 1 to n inclusive.

% Side note: how to use randi, though we don't use it here
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

%% Efficiency Matrices to save time
% Because fsolve takes a long time in iteration to solve for labor_1,
% labor_2, I create the matrices for them as well as for consumption_1 and
% consumption_2 to retrieve from later.

% Note we have to put "kPrime" and "k" on the first two dimensions in order not
% to cause any problem of dimensionality during interpolation

mLabor_1Fsolve = NaN(Nk,Nk,Na); % kPrime,k,[a_1,a_2]
mLabor_2Fsolve = NaN(Nk,Nk,Na);
mConsumption_1Fsolve = NaN(Nk,Nk,Na);
mConsumption_2Fsolve = NaN(Nk,Nk,Na);
mCurrentUtilityFsolve = NaN(Nk,Nk,Na);
% mMarginalUtilityTodayFsolve = zeros(Nk,Na,Nk);
% mMarginalUtilityTomorrowFsolve = ;

% laborInitial=[labor_1_SteadyState,labor_2_SteadyState];
opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','off');
options = optimset('Display', 'off');

tic
for ia = 1:Na
    a_1 = mGrid_a1a2(ia,1);
    a_2 = mGrid_a1a2(ia,2);
    
    
    parfor ik = 1:Nk
        laborInitial=[labor_1_SteadyState,labor_2_SteadyState];
        k = vGrid_k(ik);
        
        for ikPrime = 1:Nk
            kPrime = vGrid_k(ikPrime);
            
            vLaborFsolve = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
%                 [labor, fval] = fminbnd(@(labor) ...
%                     -laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta),...
%                     labor_1_SteadyState,labor_2_SteadyState,options)
            if isreal(vLaborFsolve)==0
                break
            end
            mLabor_1Fsolve(ikPrime,ik,ia) = vLaborFsolve(1);
            mLabor_2Fsolve(ikPrime,ik,ia) = vLaborFsolve(2);
            mConsumption_1Fsolve(ikPrime,ik,ia) = consumptionFunction1(a_1,k,kPrime,vLaborFsolve(1),aalphaK,aalphaL,ddelta);
            mConsumption_2Fsolve(ikPrime,ik,ia) = consumptionFunction2(a_2,vLaborFsolve(2));
            mCurrentUtilityFsolve(ikPrime,ik,ia) = utilityFunction(mConsumption_1Fsolve(ikPrime,ik,ia),mConsumption_2Fsolve(ikPrime,ik,ia),vLaborFsolve(1),vLaborFsolve(2),mmu_1,mmu_2);
%             mMarginalUtilityTodayFsolve(ikPrime,ia,ik) = mmu_1 * (mConsumption_2Fsolve(ikPrime,ia,ik))^mmu_2 * (mConsumption_1Fsolve(ikPrime,ia,ik))^(mmu_1-1);
            laborInitial=[vLaborFsolve(1),vLaborFsolve(2)];
        end
        
    end
    
    fprintf('Progress = %d%% \n', round(ia/Na*100)); 

end
toc
elapsedTimeMinutes=toc/60;


save('efficiencyMatricesNk500StartingWithNaNMatrix','mLabor_1Fsolve','mLabor_2Fsolve','mConsumption_1Fsolve','mConsumption_2Fsolve','mCurrentUtilityFsolve','elapsedTimeMinutes')
% 历时 2333.831306s 秒。% my computer
% 历时 2152.400222 秒。 % my computer parfor 2019-11-30 22:15:53

% Elapsed time is 385.744013 seconds. % Lab computer
% Elapsed time is 1631.222194 seconds. %Lab computer
% Elapsed time is 1543.995917 seconds.%Lab computer

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

%% Draw random grid points

% Most generally, we should draw random draws from all possible states, i.e., CAPITAL * SHOCKS.
% If we have 500 points for capital, 15 points for shocks, then altogether
% we have 7500 states, from which we draw.

% But we are told in the homework that the grid for capital should be fixed. 
% And we only have two dimensions of shocks. 
% That means the other dimension also needs to be fixed. 
% So we are going to take draws from different level of capital as well as from different shocks.

% Note that we want to keep the first and the last grid point of capital
NkDraws = ceil(Nk/2);
NkStochasticGrid = NkDraws + 2;
seedKDraws=rng;
rng(seedKDraws);
ikStochasticGrid = randperm(Nk-2,NkDraws)+1;%a ROW vector containing numberOfDraws unique integers selected randomly from 2 to Nk-1 inclusive.
ikStochasticGrid = sort(ikStochasticGrid);
ikStochasticGrid = [1,ikStochasticGrid,Nk];

NaDraws = floor(Na/1);
seedADraws=rng;
rng(seedADraws);
iaStochasticGrid = randperm(Na,NaDraws);%a ROW vector containing numberOfDraws unique integers selected randomly from 1 to Nk inclusive.
iaStochasticGrid = sort(iaStochasticGrid);

% Get the probability matrix for the random draws
mProb_a1a2StochasticGrid = zeros(NaDraws,NaDraws);
for iiaPrime = 1:length(iaStochasticGrid)
    iaPrime = iaStochasticGrid(iiaPrime);
     mProb_a1a2StochasticGrid(:,iiaPrime) = mProb_a1a2(iaStochasticGrid,iaPrime);
end
% Then normalize it to a proper probability matrix
mProb_a1a2StochasticGrid = mProb_a1a2StochasticGrid./sum(mProb_a1a2StochasticGrid,2);

%% Then do the regular Value Function Iteration using value function calculated above as the first guess

%% Required matrices and vectors

mValue0                 = utilitySteadyState.*ones(NkStochasticGrid,NaDraws); 
mValue                  = zeros(NkStochasticGrid,NaDraws);
mKPolicy                = zeros(NkStochasticGrid,NaDraws);
mLaborPolicy_1          = zeros(NkStochasticGrid,NaDraws);
mLaborPolicy_2          = zeros(NkStochasticGrid,NaDraws); 
mConsumptionPolicy_1    = zeros(NkStochasticGrid,NaDraws);
mConsumptionPolicy_2    = zeros(NkStochasticGrid,NaDraws);

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
while iteration <= maxIter  &&  (mDifference(iteration) > tolerance)
        
    expectedValue0 = mValue0 * mProb_a1a2StochasticGrid';% row: kPrime, column: a, i.e., expected value if today's shock is a and I choose k as tomorrow's capital stock
    

    for ia = iaStochasticGrid(1:end)
        a_1 = mGrid_a1a2(ia,1);
        a_2 = mGrid_a1a2(ia,2);
        iia = sum(ia>=iaStochasticGrid);

        iikPrime = 1;

        for ik = ikStochasticGrid(1:end)

            k = vGrid_k(ik);
            iik = sum(ik>=ikStochasticGrid);

%             valueHighSoFar = -1000.0;
%             kChoice  = vGrid_k(1);
% %                     laborInitial=[labor_1_SteadyState,labor_2_SteadyState];
% 
%             for ikPrime = ikStochasticGrid(iikPrime:end)
%                 kPrime = vGrid_k(ikPrime);
% 
%                 [valueProvisional,labor_1,labor_2,consumption_1,consumption_2] = valueFunction_stochasticGrid...
%                     (kPrime,ikPrime,ik,k,ia,a_1,a_2,iik,iia,iikPrime,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL,laborInitial);
% 
%                 if (valueProvisional>valueHighSoFar)
%                     valueHighSoFar = valueProvisional;
%                     kChoice = vGrid_k(ikPrime);
%                     iikPrime = sum(ikPrime >= ikStochasticGrid);
% %                             laborInitial = [labor_1,labor_2];
% 
%                     mValue(iik,iia) = valueHighSoFar;
%                     mKPolicy(iik,iia) = kChoice;
%                     mLaborPolicy_1(iik,iia) = labor_1;
%                     mLaborPolicy_2(iik,iia) = labor_2;
%                     mConsumptionPolicy_1(iik,iia) = consumption_1;
%                     mConsumptionPolicy_2(iik,iia) = consumption_2;
%                 else
%                     break
%                 end
% 
%             end

           mCurrentUtility = mCurrentUtilityFsolve(ikStochasticGrid,ik,ia);
            
            [value,iikPrime] = max((1-bbeta)*mCurrentUtility + bbeta * expectedValue0(:,iia));
            ikPrime = ikStochasticGrid(iikPrime);
            mKPolicyIndex(iik,iia) = ikPrime;
            mKPolicy(iik,iia) = vGrid_k(ikPrime);
            mValue(iik,iia) = value;     
            
 %% If we can retrieve other outputs of the objective function of fminbnd, it would be much easier and FASTER to get the labor and consumption policy function
            if mDifference(iteration) < tolerance*10 % Store the policy functions in case it's the last iteration
%                 vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
%                 mLaborPolicy_1(ik,ia) = vLabor(1);
%                 mLaborPolicy_2(ik,ia) = vLabor(2);
%                 mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
%                 mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,vLabor(2));
%                 laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
%                 inputs.laborInitial = laborInitial;
                mLaborPolicy_1(iik,iia) = mLabor_1Fsolve(ikPrime,ik,ia);
                mLaborPolicy_2(iik,iia) = mLabor_2Fsolve(ikPrime,ik,ia);
                mConsumptionPolicy_1(iik,iia) = mConsumption_1Fsolve(ikPrime,ik,ia);
                mConsumptionPolicy_2(iik,iia) = mConsumption_2Fsolve(ikPrime,ik,ia);;
            end
            
%             mValue(iik,iia) = valueHighSoFar;
%             mKPolicy(iik,iia) = kChoice;
%             mLaborPolicy_1(iik,iia) = labor_1;
%             mLaborPolicy_2(iik,iia) = labor_2;
%             mConsumptionPolicy_1(iik,iia) = consumption_1;
%             mConsumptionPolicy_2(iik,iia) = consumption_2;

        end
    end
    iteration = iteration + 1;
    mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
    mValue0         = mValue;

    if mod(iteration,10) == 2
        fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 
    end

end

toc

fprintf(' Convergence achieved. Total Number of Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 


%% For accuracy test, compute the euler equation error

% Let's use LINEAR INTERPOLATION for tomorrow's values

errorEulerEquationLinearInterpolation = eulerEquationErrorFunction(NkStochasticGrid,vGrid_k(ikStochasticGrid),mKPolicy,mLaborPolicy_1,mConsumptionPolicy_1,mConsumptionPolicy_2,NaDraws,mGrid_a1a2,mProb_a1a2,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL);
errorEulerEquationLinearInterpolationDecimalLog = log10( errorEulerEquationLinearInterpolation );

[kk,aa]=meshgrid(vGrid_k(ikStochasticGrid), mGrid_a1a2(iaStochasticGrid,1));

figure
mesh(kk, aa, errorEulerEquationLinearInterpolationDecimalLog');

title('Euler Equation Error $log_{10}$ Linear Interpolation - Stochastic Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('error','interpreter','latex')
xlim([min(vGrid_k(ikStochasticGrid)),max(vGrid_k(ikStochasticGrid))])
ylim([min(mGrid_a1a2((iaStochasticGrid),1)),max(mGrid_a1a2((iaStochasticGrid),1))])
savefig('q3_eulerEquationErrorLinearInterpolation_stochastic_grid')

save ShashaWang_JFV_PS1_500_stochastic_capital_grid_points_valueFunctionIteration

% Lab computer 12/2/2019 01:56
%  Iteration:  1, Sup diff: 0.00875014
%  Iteration: 11, Sup diff: 0.00067335
%  Iteration: 21, Sup diff: 0.00039615
%  Iteration: 31, Sup diff: 0.00023352
%  Iteration: 41, Sup diff: 0.00013843
%  Iteration: 51, Sup diff: 0.00008258
%  Iteration: 61, Sup diff: 0.00004958
%  Iteration: 71, Sup diff: 0.00002996
%  Iteration: 81, Sup diff: 0.00001822
%  Iteration: 91, Sup diff: 0.00001115
%  Iteration: 101, Sup diff: 0.00000686
%  Iteration: 111, Sup diff: 0.00000425
%  Iteration: 121, Sup diff: 0.00000265
%  Iteration: 131, Sup diff: 0.00000169
%  Iteration: 141, Sup diff: 0.00000112
% Elapsed time is 1.447077 seconds.
%  Convergence achieved. Total Number of Iteration: 144, Sup diff: 0.00000099

%% figures for Value Function Iteration with a Fixed Grid

figure
mesh(kk, aa, mValue');

title('Value - Stochastic Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Value','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_value_stochastic_grid')

figure
mesh(kk, aa, mKPolicy');

title('Policy for Next Period Capital - Stochastic Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_kPolicy_stochastic_grid')

figure
mesh(kk, aa, mLaborPolicy_1');

title('Policy for Good 1 Labor - Stochastic Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_1_stochastic_grid')

figure
mesh(kk, aa, mLaborPolicy_2');

title('Policy for Good 2 Labor - Stochastic Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_2_stochastic_grid')

figure
mesh(kk, aa, mConsumptionPolicy_1');

title('Policy for Good 1 Consumption - Stochastic Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_1_stochastic_grid')

figure
mesh(kk, aa, mConsumptionPolicy_2');

title('Policy for Good 2 Consumption - Stochastic Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_2_stochastic_grid')
