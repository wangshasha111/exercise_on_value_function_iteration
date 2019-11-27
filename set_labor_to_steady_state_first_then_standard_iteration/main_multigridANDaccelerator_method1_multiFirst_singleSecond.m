% Homework 1 Econ714 JFV
% A two good production economy
% Written by Shasha Wang

main_setup % parameters, shocks

%% Value function iteration with MULTIGRID & ACCELERATOR

% There are TWO ways to do the Multigrid, since we have two steps in the value function iteration:
% Step 1: Run a Semi-value function iteration, i.e., set the labor level
% equal to the steady state and then get the converged value matrix.
% Step 2: Use the matrix as the first guess of the subsequent qualified VFI 
% where we choose all the choice variables.

% So we have two ways of doing the multigrid: (we can put MULTIGRID algorithm in either of the two steps)
% We could first run multigrid in the semi-value function iteration and
% then use the converged value matrix as the first guess for 
% the real qualified value function iteration for the maximum grid
% points.

% Or we can first run single grid (minimum grid points) for the semi-value function iteration and
% then use the converged value matrix as the first guess for 
% the real qualified value function iteration which is done multigrid.
 
% Intuitively, second method is more reasonable since we are doing "multigrid" on the VFI that we actually want to do. 

% But let's try both. Method 2 is in another file, in which I generalized
% the algorithm to the most extent.

%% Method 1

%% Multigrid on PreIteration and SingleGrid on RealIteration

%% Step 1: Compute the Steady State
    
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
% Nk = 250;
% vGrid_k = linspace(kMin,kMax,Nk)';

%% Step 3. Pre-Iteration: Set labor to steady state

% Iterate on the Value function implied by the Social Planner¡¯s Problem using linear interpolation
% until the change in the sup norm between to iterations is less than 1e-6.
% Compute the Policy function.
% Describe the responses of the economy to a technology shock.

kGridLength = [100,500,5000]; % number of points in grid for capital 
Nk = max(kGridLength);

%% Create cells to restore results of Pre-Iteration
% For consistency I use cells, even though for iteration and mDifference I
% could use matrices.
% Now I realized I could've used cell to combine the two exogenous shocks

cValueMultigridPreIteration = cell(length(kGridLength),1);
cKPolicyMultigridPreIteration = cell(length(kGridLength),1);
% cLaborPolicy_1MultigridPreIteration = cell(length(kGridLength),1);
% cLaborPolicy_2MultigridPreIteration = cell(length(kGridLength),1);
cConsumptionPolicy_1MultigridPreIteration = cell(length(kGridLength),1);
cConsumptionPolicy_2MultigridPreIteration = cell(length(kGridLength),1);
cIterationMultigridPreIteration  = cell(length(kGridLength),1);
cDifferenceMultigridPreIteration  = cell(length(kGridLength),1);

%% Required matrices and vectors

mValue0        = utilitySteadyState.*ones(kGridLength(1),Na);
inputs.mValue0 = mValue0;
mValue         = zeros(kGridLength(1),Na);
mKPolicy        = zeros(kGridLength(1),Na);
% mLaborPolicy_1 = zeros(kGridLength(1),Na);
% mLaborPolicy_2 = zeros(kGridLength(1),Na); 
mConsumptionPolicy_1   = zeros(kGridLength(1),Na);
mConsumptionPolicy_2   = zeros(kGridLength(1),Na);

tic 

for i=1:length(kGridLength)
    vGrid_k           = linspace(kMin,kMax,kGridLength(i))';
    inputs.vGrid_k = vGrid_k;

%% Main iteration

    maxIter = 10000;
    tolerance     = 1e-6;

    iteration       = 1;
    mDifference = zeros(maxIter,1);
    mDifference(iteration) = 100;

    options = optimset('Display', 'off');
% opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','off');
% laborInitial=[labor_1_SteadyState,labor_2_SteadyState];


% Finding value and policy functions numerically
    while iteration <= maxIter  ...% make sure the last iteration does the maximization
             && ( (mDifference(iteration) > tolerance)  | (mDifference(iteration) <= tolerance && mod(iteration,10)~=2)) 
        
        expectedValue0 = mValue0 * mProb_a1a2';% row: kPrime, column: a, i.e., expected value if today's shock is a and I choose k as tomorrow's capital stock

        if (mDifference(iteration) > tolerance)    
            for ia = 1:Na
                a_1 = mGrid_a1a2(ia,1);
                a_2 = mGrid_a1a2(ia,2);

                for ik = 1:kGridLength(i)
                    k = vGrid_k(ik);

                    if mod(iteration,10) == 1

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

                for ik = 1:kGridLength(i)
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
    fprintf(' Convergence Achieved. \n Number of Grid Points: %2.0f\n Iteration: %2.0f, Sup difference: %2.8f\n ', kGridLength(i),iteration-1, mDifference(iteration)); 

%% Store result for this grid point
    cValueMultigridPreIteration{i}= mValue;
    cKPolicyMultigridPreIteration{i} = mKPolicy;
%     cLaborPolicy_1Multigrid{i} = mLaborPolicy_1;
%     cLaborPolicy_2Multigrid{i} = mLaborPolicy_2;
    cConsumptionPolicy_1MultigridPreIteration{i} = mConsumptionPolicy_1;
    cConsumptionPolicy_2MultigridPreIteration{i} = mConsumptionPolicy_2;
    cIterationMultigridPreIteration{i} = iteration - 1;
    cDifferenceMultigridPreIteration{i} = mDifference(2:iteration);

    if i ~= length(kGridLength)
        mValue0 = interp1(vGrid_k, mValue,linspace(kMin, kMax, kGridLength(i+1))');
        mValue  = mValue0; % Just for consistency of dimensionality. mValue can be anything you want, as long as the dimensionality is right
    %         kPolicy = interp1(grid_k,kPolicy,linspace(kMin, kMax, kGridLength(i+1)));
        mKPolicy         = zeros(kGridLength(i+1),Na);
%         mLaborPolicy_1         = zeros(kGridLength(i+1),Na);
%         mLaborPolicy_2         = zeros(kGridLength(i+1),Na);
        mConsumptionPolicy_1         = zeros(kGridLength(i+1),Na);
        mConsumptionPolicy_2         = zeros(kGridLength(i+1),Na);
        
        
    end

end

toc

fprintf(' Convergence Achieved. \n Number of Grid Points: %2.0f\n Iteration: %2.0f, Sup difference: %2.8f\n ', kGridLength(i),iteration-1, mDifference(iteration)); 

save ShashaWang_JFV_PS1_250_capital_grid_points_valueFunctionIteration_setLaborToSteadyState_multigridANDaccelerator_method1

%  Iteration:  1, Sup diff: 0.00530891
%  Iteration: 11, Sup diff: 0.00073227
%  Iteration: 21, Sup diff: 0.00019339
%  Iteration: 31, Sup diff: 0.00010574
%  Iteration: 41, Sup diff: 0.00011728
%  Iteration: 51, Sup diff: 0.00003955
%  Iteration: 61, Sup diff: 0.00002003
%  Iteration: 71, Sup diff: 0.00001179
%  Iteration: 81, Sup diff: 0.00000718
%  Iteration: 91, Sup diff: 0.00000440
%  Iteration: 101, Sup diff: 0.00000271
%  Iteration: 111, Sup diff: 0.00000167
%  Iteration: 121, Sup diff: 0.00000104
%  Iteration: 123, Sup diff: 0.00000094
%  Convergence Achieved. 
%  Number of Grid Points: 100
%  Iteration: 123, Sup difference: 0.00000094

%  Iteration:  1, Sup diff: 0.00000091
%  Convergence Achieved. 
%  Number of Grid Points: 500
%  Iteration:  1, Sup difference: 0.00000091
 
%  Iteration:  1, Sup diff: 0.00046054
%  Iteration: 11, Sup diff: 0.00029520
%  Iteration: 21, Sup diff: 0.00022697
%  Iteration: 31, Sup diff: 0.00016110
%  Iteration: 41, Sup diff: 0.00025614
%  Iteration: 51, Sup diff: 0.00038761
%  Iteration: 61, Sup diff: 0.00014773
%  Iteration: 71, Sup diff: 0.00011937
%  Iteration: 81, Sup diff: 0.00015672
%  Iteration: 91, Sup diff: 0.00014987
%  Iteration: 93, Sup diff: 0.00000093
%  Convergence Achieved. 
%  Number of Grid Points: 5000
%  Iteration: 93, Sup difference: 0.00000093

%  Elapsed time is 1358.047397 seconds.
%  Convergence achieved. Total Number of Iteration: 93, Sup diff: 0.00000093
%  Convergence Achieved. 
%  Number of Grid Points: 5000
%  Iteration: 93, Sup difference: 0.00000093


%% Then do the real Value Function Iteration using value function calculated above as the first guess

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

title('Euler Equation Error $log_{10}$ Linear Interpolation - Multigrid & Accelerator - Multi First Single Second','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('error','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_eulerEquationErrorLinearInterpolation_multigridANDaccelerator_method1_multiFirst_singleSecond')

save ShashaWang_JFV_PS1_250_capital_grid_points_valueFunctionIteration_setLaborToSteadyState_thenDoRealValueFunctionIteration_multigridANDaccelerator_method1
%  Iteration:  1, Sup diff: 0.00205348
%  Iteration:  2, Sup diff: 0.00195390
%  Iteration:  3, Sup diff: 0.00186060
%  Iteration:  4, Sup diff: 0.00108450
%  Iteration:  5, Sup diff: 0.00084524
%  Iteration:  6, Sup diff: 0.00080195
%  Iteration:  7, Sup diff: 0.00076164
%  Iteration:  8, Sup diff: 0.00072402
%  Iteration:  9, Sup diff: 0.00068884
%  Iteration: 10, Sup diff: 0.00065587
%  Iteration: 11, Sup diff: 0.00062520
%  Iteration: 12, Sup diff: 0.00059589
%  Iteration: 13, Sup diff: 0.00056839
%  Iteration: 14, Sup diff: 0.00054250
%  Iteration: 15, Sup diff: 0.00051804
%  Iteration: 16, Sup diff: 0.00049490
%  Iteration: 17, Sup diff: 0.00047298
%  Iteration: 18, Sup diff: 0.00045219
%  Iteration: 19, Sup diff: 0.00043246
%  Iteration: 20, Sup diff: 0.00041371
%  Iteration: 21, Sup diff: 0.00039588
%  Iteration: 22, Sup diff: 0.00037891
%  Iteration: 23, Sup diff: 0.00036274
%  Iteration: 24, Sup diff: 0.00034734
%  Iteration: 25, Sup diff: 0.00033264
%  Iteration: 26, Sup diff: 0.00031863
%  Iteration: 27, Sup diff: 0.00030524
%  Iteration: 28, Sup diff: 0.00029246
%  Iteration: 29, Sup diff: 0.00028025
%  Iteration: 30, Sup diff: 0.00026858
%  Iteration: 31, Sup diff: 0.00025743
%  Iteration: 32, Sup diff: 0.00024675
%  Iteration: 33, Sup diff: 0.00023654
%  Iteration: 34, Sup diff: 0.00022677
%  Iteration: 35, Sup diff: 0.00021742
%  Iteration: 36, Sup diff: 0.00020847
%  Iteration: 37, Sup diff: 0.00019990
%  Iteration: 38, Sup diff: 0.00019169
%  Iteration: 39, Sup diff: 0.00018383
%  Iteration: 40, Sup diff: 0.00017629
%  Iteration: 41, Sup diff: 0.00016908
%  Iteration: 42, Sup diff: 0.00016216
%  Iteration: 43, Sup diff: 0.00015554
%  Iteration: 44, Sup diff: 0.00014919
%  Iteration: 45, Sup diff: 0.00014310
%  Iteration: 46, Sup diff: 0.00013726
%  Iteration: 47, Sup diff: 0.00013167
%  Iteration: 48, Sup diff: 0.00012631
%  Iteration: 49, Sup diff: 0.00012116
%  Iteration: 50, Sup diff: 0.00011623
%  Iteration: 51, Sup diff: 0.00011151
%  Iteration: 52, Sup diff: 0.00010697
%  Iteration: 53, Sup diff: 0.00010263
%  Iteration: 54, Sup diff: 0.00009846
%  Iteration: 55, Sup diff: 0.00009446
%  Iteration: 56, Sup diff: 0.00009063
%  Iteration: 57, Sup diff: 0.00008695
%  Iteration: 58, Sup diff: 0.00008342
%  Iteration: 59, Sup diff: 0.00008004
%  Iteration: 60, Sup diff: 0.00007679
%  Iteration: 61, Sup diff: 0.00007368
%  Iteration: 62, Sup diff: 0.00007069
%  Iteration: 63, Sup diff: 0.00006783
%  Iteration: 64, Sup diff: 0.00006508
%  Iteration: 65, Sup diff: 0.00006245
%  Iteration: 66, Sup diff: 0.00005992
%  Iteration: 67, Sup diff: 0.00005749
%  Iteration: 68, Sup diff: 0.00005517
%  Iteration: 69, Sup diff: 0.00005293
%  Iteration: 70, Sup diff: 0.00005079
%  Iteration: 71, Sup diff: 0.00004874
%  Iteration: 72, Sup diff: 0.00004677
%  Iteration: 73, Sup diff: 0.00004488
%  Iteration: 74, Sup diff: 0.00004306
%  Iteration: 75, Sup diff: 0.00004132
%  Iteration: 76, Sup diff: 0.00003965
%  Iteration: 77, Sup diff: 0.00003805
%  Iteration: 78, Sup diff: 0.00003651
%  Iteration: 79, Sup diff: 0.00003504
%  Iteration: 80, Sup diff: 0.00003362
%  Iteration: 81, Sup diff: 0.00003226
%  Iteration: 82, Sup diff: 0.00003096
%  Iteration: 83, Sup diff: 0.00002971
%  Iteration: 84, Sup diff: 0.00002851
%  Iteration: 85, Sup diff: 0.00002736
%  Iteration: 86, Sup diff: 0.00002626
%  Iteration: 87, Sup diff: 0.00002520
%  Iteration: 88, Sup diff: 0.00002418
%  Iteration: 89, Sup diff: 0.00002320
%  Iteration: 90, Sup diff: 0.00002227
%  Iteration: 91, Sup diff: 0.00002137
%  Iteration: 92, Sup diff: 0.00002051
%  Iteration: 93, Sup diff: 0.00001968
%  Iteration: 94, Sup diff: 0.00001889
%  Iteration: 95, Sup diff: 0.00001812
%  Iteration: 96, Sup diff: 0.00001739
%  Iteration: 97, Sup diff: 0.00001669
%  Iteration: 98, Sup diff: 0.00001602
%  Iteration: 99, Sup diff: 0.00001537
%  Iteration: 100, Sup diff: 0.00001475
%  Iteration: 101, Sup diff: 0.00001416
%  Iteration: 102, Sup diff: 0.00001359
%  Iteration: 103, Sup diff: 0.00001304
%  Iteration: 104, Sup diff: 0.00001252
%  Iteration: 105, Sup diff: 0.00001201
%  Iteration: 106, Sup diff: 0.00001153
%  Iteration: 107, Sup diff: 0.00001106
%  Iteration: 108, Sup diff: 0.00001062
%  Iteration: 109, Sup diff: 0.00001019
%  Iteration: 110, Sup diff: 0.00000978
%  Iteration: 111, Sup diff: 0.00000939
%  Iteration: 112, Sup diff: 0.00000901
%  Iteration: 113, Sup diff: 0.00000865
%  Iteration: 114, Sup diff: 0.00000830
%  Iteration: 115, Sup diff: 0.00000796
%  Iteration: 116, Sup diff: 0.00000764
%  Iteration: 117, Sup diff: 0.00000734
%  Iteration: 118, Sup diff: 0.00000704
%  Iteration: 119, Sup diff: 0.00000676
%  Iteration: 120, Sup diff: 0.00000649
%  Iteration: 121, Sup diff: 0.00000623
%  Iteration: 122, Sup diff: 0.00000597
%  Iteration: 123, Sup diff: 0.00000573
%  Iteration: 124, Sup diff: 0.00000550
%  Iteration: 125, Sup diff: 0.00000528
%  Iteration: 126, Sup diff: 0.00000507
%  Iteration: 127, Sup diff: 0.00000487
%  Iteration: 128, Sup diff: 0.00000467
%  Iteration: 129, Sup diff: 0.00000448
%  Iteration: 130, Sup diff: 0.00000430
%  Iteration: 131, Sup diff: 0.00000413
%  Iteration: 132, Sup diff: 0.00000396
%  Iteration: 133, Sup diff: 0.00000380
%  Iteration: 134, Sup diff: 0.00000365
%  Iteration: 135, Sup diff: 0.00000350
%  Iteration: 136, Sup diff: 0.00000336
%  Iteration: 137, Sup diff: 0.00000323
%  Iteration: 138, Sup diff: 0.00000310
%  Iteration: 139, Sup diff: 0.00000297
%  Iteration: 140, Sup diff: 0.00000285
%  Iteration: 141, Sup diff: 0.00000274
%  Iteration: 142, Sup diff: 0.00000263
%  Iteration: 143, Sup diff: 0.00000252
%  Iteration: 144, Sup diff: 0.00000242
%  Iteration: 145, Sup diff: 0.00000233
%  Iteration: 146, Sup diff: 0.00000223
%  Iteration: 147, Sup diff: 0.00000214
%  Iteration: 148, Sup diff: 0.00000206
%  Iteration: 149, Sup diff: 0.00000197
%  Iteration: 150, Sup diff: 0.00000189
%  Iteration: 151, Sup diff: 0.00000182
%  Iteration: 152, Sup diff: 0.00000175
%  Iteration: 153, Sup diff: 0.00000168
%  Iteration: 154, Sup diff: 0.00000161
%  Iteration: 155, Sup diff: 0.00000154
%  Iteration: 156, Sup diff: 0.00000148
%  Iteration: 157, Sup diff: 0.00000142
%  Iteration: 158, Sup diff: 0.00000136
%  Iteration: 159, Sup diff: 0.00000131
%  Iteration: 160, Sup diff: 0.00000126
%  Iteration: 161, Sup diff: 0.00000121
%  Iteration: 162, Sup diff: 0.00000116
%  Iteration: 163, Sup diff: 0.00000111
%  Iteration: 164, Sup diff: 0.00000107
%  Iteration: 165, Sup diff: 0.00000102
%  Iteration: 166, Sup diff: 0.00000098
%  Iteration: 167, Sup diff: 0.00000094
% Elapsed time is 13981.929367 seconds. 3.88 hours
%  Convergence achieved. Total Number of Iteration: 167, Sup diff: 0.00000094

%% figures for Value Function Iteration with a Fixed Grid

figure
mesh(kk, aa, mValue');

title('Value - Multigrid & Accelerator - Multi First Single Second','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Value','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_value_multigridANDaccelerator_method1_multiFirst_singleSecond')

figure
mesh(kk, aa, mKPolicy');

title('Policy for Next Period Capital - Multigrid & Accelerator - Multi First Single Second','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_kPolicy_multigridANDaccelerator_method1_multiFirst_singleSecond')

figure
mesh(kk, aa, mLaborPolicy_1');

title('Policy for Good 1 Labor - Multigrid & Accelerator - Multi First Single Second','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_1_multigridANDaccelerator_method1_multiFirst_singleSecond')

figure
mesh(kk, aa, mLaborPolicy_2');

title('Policy for Good 2 Labor - Multigrid & Accelerator - Multi First Single Second','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_2_multigridANDaccelerator_method1_multiFirst_singleSecond')

figure
mesh(kk, aa, mConsumptionPolicy_1');

title('Policy for Good 1 Consumption - Multigrid & Accelerator - Multi First Single Second','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_1_multigridANDaccelerator_method1_multiFirst_singleSecond')

figure
mesh(kk, aa, mConsumptionPolicy_2');

title('Policy for Good 2 Consumption- Multigrid & Accelerator - Multi First Single Second','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_2_multigridANDaccelerator_method1_multiFirst_singleSecond')

