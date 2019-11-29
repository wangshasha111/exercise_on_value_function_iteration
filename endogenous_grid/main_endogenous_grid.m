% Homework 1 Econ714 JFV
% A two good production economy
% Written by Shasha Wang
% Endogenous grid scheme 
% Reference: JMMT 2016

%% 4. Value Function Iteration with an Endogenous Grid
% Repeat previous exercise with an endogenous grid.

%% Algorithm explained
% Step 1: set labor to the steady state level and use endogenous grid method
% to iteration the value function till convergence. Specifically,
% % (This should be easily matrixized)
% % Guess mValue - NkPrime by Na, i.e., tomorrow's expected value matrix and compute its derivative
% % Given mValue(ikPrime,:), by the first order condition, use the derivative of V(kPrime,a) for all a's (denoted by V(kPrime,:)) to get consumptions while keeping labor at steady state level.
% % Given consumption, we know current utility u(kPrime,:). Recall we know already mValue(ikPrime,:) by guessing.
% % Sum up V(kPrime,:) and current utility u(kPrime,:), we get V(output,a) for all a's; remember output = consumption + kPrime
% % Compute tomorrow's expected value of V(output,:), i.e., bbeta * V(output,:) * mProb_a1a2', which is the updated V(kPrime,:)
% % Update mValue using it, and iterate until convergence

% Step 2: use the converged value function as the initial guess for the real
% qualified value function iteration, as in the fixed grid case.

main_setup % parameters, shocks

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

kSpread = 0.3;
kMax = kSteadyState * (1 + kSpread);
kMin = kSteadyState * (1 - kSpread);
Nk = 250;
vGrid_k = linspace(kMin,kMax,Nk)';
inputs.vGrid_k = vGrid_k;

%% Step 3. Endogenous Grid: Set labor to steady state

