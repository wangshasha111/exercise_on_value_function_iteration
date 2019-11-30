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
% qualified value function iteration, as in the fixed grid case. Remember
% to use accelerator

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
                
            mLabor_1Fsolve(ikPrime,ik,ia) = vLaborFsolve(1);
            mLabor_2Fsolve(ikPrime,ik,ia) = vLaborFsolve(2);
            mConsumption_1Fsolve(ikPrime,ik,ia) = consumptionFunction1(a_1,k,kPrime,vLaborFsolve(1),aalphaK,aalphaL,ddelta);
            mConsumption_2Fsolve(ikPrime,ik,ia) = consumptionFunction2(a_2,vLaborFsolve(2));
            mCurrentUtilityFsolve(ikPrime,ik,ia) = utilityFunction(mConsumption_1Fsolve(ikPrime,ik,ia),mConsumption_2Fsolve(ikPrime,ik,ia),vLaborFsolve(1),vLaborFsolve(2),mmu_1,mmu_2);
%             mMarginalUtilityTodayFsolve(ikPrime,ia,ik) = mmu_1 * (mConsumption_2Fsolve(ikPrime,ia,ik))^mmu_2 * (mConsumption_1Fsolve(ikPrime,ia,ik))^(mmu_1-1);
            laborInitial=[vLaborFsolve(1),vLaborFsolve(2)];
        end
        
    end
    
    fprintf('Progress calculating period utility = %d%% \n', ia/Na*100); 

end
toc
elapsedTimeMinutes=toc/60;

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

save('efficiencyMatricesNk250','mLabor_1Fsolve','mLabor_2Fsolve','mConsumption_1Fsolve','mConsumption_2Fsolve','mCurrentUtilityFsolve','elapsedTimeMinutes')
% 历时 2333.831306s 秒。% my computer
% 历时 2152.400222 秒。 % my computer parfor

%% Step 3. Value Function Iteration without interpolation 
% to get a good initial guess for endogenous grid

mValue0        = utilitySteadyState.*ones(Nk,Na);
inputs.mValue0 = mValue0;
mValue         = zeros(Nk,Na);
mKPolicy        = zeros(Nk,Na);
mKPolicyIndex= zeros(Nk,Na);
mLaborPolicy_1 = zeros(Nk,Na);
mLaborPolicy_2 = zeros(Nk,Na); 
mConsumptionPolicy_1   = zeros(Nk,Na);
mConsumptionPolicy_2   = zeros(Nk,Na);

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
            
%             if ik ==1
%                         
%                 [kPrime, fval] = fminbnd(@(kPrime) ...
%                     -valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
%                     vGrid_k(1),min(1.2*vGrid_k(ik),vGrid_k(end)),options);       
%             else
%                 [kPrime, fval] = fminbnd(@(kPrime) ...
%                     -valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
%                     mKPolicy((ik-1),ia),min(1.2*vGrid_k(ik),vGrid_k(end)),options);
%             end
            
%             for ikPrime = 1:NkPrime
            mCurrentUtility = mCurrentUtilityFsolve(:,ik,ia);
            
            [value,ikPrime] = max((1-bbeta)*mCurrentUtility + bbeta * expectedValue0(:,ia));
            mKPolicyIndex(ik,ia) = ikPrime;
            mKPolicy(ik,ia) = vGrid_k(ikPrime);
            mValue(ik,ia) = value;     
            
 %% If we can retrieve other outputs of the objective function of fminbnd, it would be much easier and FASTER to get the labor and consumption policy function
            if mDifference(iteration) < tolerance*10 % Store the policy functions in case it's the last iteration
%                 vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
%                 mLaborPolicy_1(ik,ia) = vLabor(1);
%                 mLaborPolicy_2(ik,ia) = vLabor(2);
%                 mConsumptionPolicy_1(ik,ia) = consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
%                 mConsumptionPolicy_2(ik,ia) = consumptionFunction2(a_2,vLabor(2));
%                 laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
%                 inputs.laborInitial = laborInitial;
                mLaborPolicy_1(ik,ia) = mLabor_1Fsolve(ikPrime,ik,ia);
                mLaborPolicy_2(ik,ia) = mLabor_2Fsolve(ikPrime,ik,ia);
                mConsumptionPolicy_1(ik,ia) = mConsumption_1Fsolve(ikPrime,ik,ia);
                mConsumptionPolicy_2(ik,ia) = mConsumption_2Fsolve(ikPrime,ik,ia);;
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
elapsed_seconds_VFI_no_interpolation=toc;
fprintf('Convergence Achieved. Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 

%% Step 4. Endogenous Grid: Set labor to steady state

%Setting up Vtilde at n=0
vGrid_kPrime = vGrid_k;
NkPrime = Nk;
mValueTildeGuess = bbeta * mValue * mProb_a1a2';
mValueTildeOld = mValueTildeGuess;

% Preallocate
% mEndogenousValueFunction   = zeros(Nk,Na);
% mEndogenousMarketResources = zeros(Nk,Na);
mValueFunctionNew          = zeros(Nk,Na);
% mValueFunctionTildeNew     = zeros(Nk,Na);
% mEndogenousConsumption     = zeros(Nk,Na);

% Setting the Market Resources
mGridOutput = zeros(Nk,Na);
for ia = 1:Na
        a_1 = mGrid_a1a2(ia,1);
        mGridOutput(:,ia) = a_1* (vGrid_k.^mmu_1) .* (labor_1_SteadyState.^mmu_2)+(1-ddelta)*vGrid_k;
end


maxIter = 10000;
tolerance     = 1e-6;

iteration       = 1;
mDifference = zeros(maxIter,1);
mDifference(iteration) = 100;

tic
while iteration<=maxIter && mDifference(iteration) > tolerance
    
    % % Guess mValue - NkPrime by Na, i.e., tomorrow's expected value matrix and compute its derivative
    mValueTildeDerivative = getDerivative(mValueTildeOld,NkPrime,Na,vGrid_kPrime);

    % % Given mValue(ikPrime,:), by the first order condition, use the derivative of V(kPrime,a) for all a's (denoted by V(kPrime,:)) to get consumptions while keeping labor at steady state level.
    mConsumption_2 = mGrid_a1a2(:,2)' * labor_2_SteadyState;
    vDenominator = mmu_1*(mConsumption_2).^mmu_2 * (1-bbeta);
    mConsumption_1 = (mValueTildeDerivative./vDenominator).^(1/(mmu_1-1));
    
    mEndogenousOutput = mConsumption_1 + vGrid_kPrime;

    % % Given consumption, we know current utility u(kPrime,:). Recall we know already mValue(ikPrime,:) by guessing.
    % mUtility = mConsumption.^mmu_1 .* (mGrid_a1a2(:,2)' * labor_2_SteadyState).^mmu_2 - (labor_1_SteadyState+labor_2_SteadyState)^2/2
    mCurrentUtility = utilityFunction(mConsumption_1,mConsumption_2,labor_1_SteadyState,labor_2_SteadyState,mmu_1,mmu_2);

    % % Sum up V(kPrime,:) and current utility u(kPrime,:), we get V(output,a) for all a's; remember output = consumption + kPrime
    mEndogenousValueFunction = (1-bbeta) * mCurrentUtility + mValueTildeOld;
    
    for ia = 1:Na
        for ikPrime = 1:NkPrime

            if isreal(mEndogenousOutput(ikPrime,ia)) == 0  ||  isreal(mEndogenousValueFunction(ikPrime,ia)) ==0
                break;
            else
%                 for iOutput = 1:Nk                    
%                     mValueFunctionNew(iOutput,ia) = interp1(mEndogenousOutput(:,ia),mEndogenousValueFunction(:,ia),mGridOutput(iOutput,ia));       

                    % Note that mere interp1 will cause error since it will interpolate outside of the range and return NaN. 
                    % But if you add 'linear','extrap', you will be able to exterpolate, but the value function will not converge.
                    % So be sure to add 'spline','extrap'
                    mValueFunctionNew(:,ia) = interp1(mEndogenousOutput(:,ia),mEndogenousValueFunction(:,ia),mGridOutput(:,ia),'spline','extrap'); 
%                 end           
            end
            
        end % end kPrime
    end % end a
                   
    % % Compute tomorrow's expected value of V(output,:), i.e., bbeta * V(output,:) * mProb_a1a2', which is the updated V(kPrime,:)
    mValueTildeNew = bbeta * mValueFunctionNew * mProb_a1a2';

    % % Update mValue using it, and iterate until convergence

    iteration = iteration + 1;
    mDifference(iteration) = max(abs(mValueTildeNew - mValueTildeOld),[],'all');
    mValueTildeOld = mValueTildeNew;
    
    if mod(iteration,10) == 2
        fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 
    end
end
toc
%  Iteration:  1, Sup diff: 0.00055379
%  Iteration: 11, Sup diff: 0.00032796
%  Iteration: 21, Sup diff: 0.00021114
%  Iteration: 31, Sup diff: 0.00013817
%  Iteration: 41, Sup diff: 0.00009100
%  Iteration: 51, Sup diff: 0.00006010
%  Iteration: 61, Sup diff: 0.00003975
%  Iteration: 71, Sup diff: 0.00002631
%  Iteration: 81, Sup diff: 0.00001742
%  Iteration: 91, Sup diff: 0.00001154
%  Iteration: 101, Sup diff: 0.00000765
%  Iteration: 111, Sup diff: 0.00000507
%  Iteration: 121, Sup diff: 0.00000336
%  Iteration: 131, Sup diff: 0.00000223
%  Iteration: 141, Sup diff: 0.00000148
%  Iteration: 151, Sup diff: 0.00000098
% 历时 73.815803 秒。

figure;
[kk,aa]=meshgrid(vGrid_k, mGrid_a1a2(:,1));
mesh(kk, aa, mValueTildeNew');

title('Value - Endogenous Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Value','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_value_endogenous_grid')

%% regular Value Function Iteration 

%% Required matrices and vectors

mValue0                 = mValueTildeNew;
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

save ShashaWang_JFV_PS1_endogenous_250_capital_grid_points_valueFunctionIteration_accelerator

%  Iteration:  1, Sup diff: 0.00043476
%  Iteration:  2, Sup diff: 0.00040420
%  Iteration:  3, Sup diff: 0.00034748
%  Iteration:  4, Sup diff: 0.00032525
%  Iteration:  5, Sup diff: 0.00030496
%  Iteration:  6, Sup diff: 0.00028638
%  Iteration:  7, Sup diff: 0.00026934
%  Iteration:  8, Sup diff: 0.00025368
%  Iteration:  9, Sup diff: 0.00023924
%  Iteration: 10, Sup diff: 0.00022592



%% figures for Value Function Iteration with a Fixed Grid

figure;
[kk,aa]=meshgrid(vGrid_k, mGrid_a1a2(:,1));
mesh(kk, aa, mValue');

title('Value - Endogenous Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Value','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_value_endogenous_grid')

figure;
mesh(kk, aa, mKPolicy');

title('Policy for Next Period Capital - Endogenous Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Next Period Capital $k\prime$','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_kPolicy_endogenous_grid')

figure;
mesh(kk, aa, mLaborPolicy_1');

title('Policy for Good 1 Labor - Endogenous Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_1_endogenous_grid')

figure;
mesh(kk, aa, mLaborPolicy_2');

title('Policy for Good 2 Labor - Endogenous Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Labor','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_laborPolicy_2_endogenous_grid')

figure;
mesh(kk, aa, mConsumptionPolicy_1');

title('Policy for Good 1 Consumption - Endogenous Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 1 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_1_endogenous_grid')

figure;
mesh(kk, aa, mConsumptionPolicy_2');

title('Policy for Good 2 Consumption - Endogenous Grid','interpreter','latex')
xlabel('Capital Stock $k$','interpreter','latex')
ylabel('Shocks $z_1$ $z_2$','interpreter','latex')
zlabel('Good 2 Consumption','interpreter','latex')
xlim([min(vGrid_k),max(vGrid_k)])
ylim([min(mGrid_a1a2(:,1)),max(mGrid_a1a2(:,1))])
savefig('q3_consumptionPolicy_2_endogenous_grid')