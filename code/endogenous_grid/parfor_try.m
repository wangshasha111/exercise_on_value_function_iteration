%% Step 5: Regular Value Function Iteration 

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

        if mod(iteration,10) == 1
        
            parfor ia = 1:Na
                a_1 = mGrid_a1a2(ia,1);
                a_2 = mGrid_a1a2(ia,2);

                % for parfor
                tempValue = zeros(Nk,1);
                tempKPolicy = zeros(Nk,1);
                tempLaborPolicy_1 = zeros(Nk,1);
                tempLaborPolicy_2 = zeros(Nk,1);
                tempConsumptionPolicy_1 = zeros(Nk,1);
                tempConsumptionPolicy_2 = zeros(Nk,1);

                for ik = 1:Nk
                    k = vGrid_k(ik);
                    laborInitial=[labor_1_SteadyState,labor_2_SteadyState];
                
%                 if mod(iteration,10) == 1

                    if ik == 1
                        [kPrime, vAux] = fminbnd(@(kPrime) ...
                            -valueFunction(kPrime,ik,k,ia,a_1,a_2,vGrid_k,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                            vGrid_k(1),min(1.2*vGrid_k(ik),vGrid_k(end)),options);

                    else
                        [kPrime, vAux] = fminbnd(@(kPrime) ...
                            -valueFunction(kPrime,ik,k,ia,a_1,a_2,vGrid_k,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                            tempKPolicy(ik-1),min(1.2*vGrid_k(ik),vGrid_k(end)),options);
                    end
                    
                    tempValue(ik) = -vAux
                    tempKPolicy(ik)= kPrime;
                    
                    vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
                    tempLaborPolicy_1(ik)= vLabor(1);
                    tempLaborPolicy_2(ik)= vLabor(2);
                    tempConsumptionPolicy_1(ik)= consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
                    tempConsumptionPolicy_2(ik)= consumptionFunction2(a_2,vLabor(2));
                    laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process

                end %k
                    mKPolicy(:,ia) = tempKPolicy;
                    mValue(:,ia) = tempValue;           

%                     vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
                    mLaborPolicy_1(:,ia) = tempLaborPolicy_1;
                    mLaborPolicy_2(:,ia) = tempLaborPolicy_2;
                    mConsumptionPolicy_1(:,ia) = tempConsumptionPolicy_1;
                    mConsumptionPolicy_2(:,ia) = tempConsumptionPolicy_2;

            end % a
        else
            for ia = 1:Na
                for ik = 1:Nk
%                     currentUtility = interp1(vGrid_k,mCurrentUtilityFsolve(:,ia,ik),mKPolicy(ik,ia));
                    currentUtility = utilityFunction(mConsumptionPolicy_1(ik,ia),mConsumptionPolicy_2(ik,ia),mLaborPolicy_1(ik,ia),mLaborPolicy_2(ik,ia),mmu_1,mmu_2);;
                    expectedValue = interp1(vGrid_k,expectedValue0(:,ia),mKPolicy(ik,ia));
                    value = (1-bbeta)*currentUtility + bbeta * expectedValue;
                    
                    mValue(ik,ia) = value;
                end %k
            end %a
                    
        end
        iteration = iteration + 1;
        mDifference(iteration) = max(abs(mValue - mValue0),[],'all');
        mValue0         = mValue;

%         if mod(iteration,10) == 2
        fprintf(' Iteration: %2.0f, Sup diff: %2.8f\n', iteration-1, mDifference(iteration)); 
%         end

    else
        parfor ia = 1:Na
            a_1 = mGrid_a1a2(ia,1);
            a_2 = mGrid_a1a2(ia,2);

            % for parfor
            tempValue = zeros(Nk,1);
            tempKPolicy = zeros(Nk,1);
            tempLaborPolicy_1 = zeros(Nk,1);
            tempLaborPolicy_2 = zeros(Nk,1);
            tempConsumptionPolicy_1 = zeros(Nk,1);
            tempConsumptionPolicy_2 = zeros(Nk,1);
                
            for ik = 1:Nk
                k = vGrid_k(ik);
                laborInitial=[labor_1_SteadyState,labor_2_SteadyState];
                
%                 if mod(iteration,10) == 1

                if ik == 1
                    [kPrime, vAux] = fminbnd(@(kPrime) ...
                        -valueFunction(kPrime,ik,k,ia,a_1,a_2,vGrid_k,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                        vGrid_k(1),min(1.2*vGrid_k(ik),vGrid_k(end)),options);

                else
                    [kPrime, vAux] = fminbnd(@(kPrime) ...
                        -valueFunction(kPrime,ik,k,ia,a_1,a_2,vGrid_k,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL),...
                        tempKPolicy(ik-1),min(1.2*vGrid_k(ik),vGrid_k(end)),options);
                end

                tempValue(ik) = -vAux
                tempKPolicy(ik)= kPrime;

                vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
                tempLaborPolicy_1(ik)= vLabor(1);
                tempLaborPolicy_2(ik)= vLabor(2);
                tempConsumptionPolicy_1(ik)= consumptionFunction1(a_1,k,kPrime,vLabor(1),aalphaK,aalphaL,ddelta);
                tempConsumptionPolicy_2(ik)= consumptionFunction2(a_2,vLabor(2));
                laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process

            end %k
            mKPolicy(:,ia) = tempKPolicy;
            mValue(:,ia) = tempValue;           

%                     vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
            mLaborPolicy_1(:,ia) = tempLaborPolicy_1;
            mLaborPolicy_2(:,ia) = tempLaborPolicy_2;
            mConsumptionPolicy_1(:,ia) = tempConsumptionPolicy_1;
            mConsumptionPolicy_2(:,ia) = tempConsumptionPolicy_2;

        end % a
            
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
