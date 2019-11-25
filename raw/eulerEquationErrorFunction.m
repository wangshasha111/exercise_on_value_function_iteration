function errorEulerEquation = eulerEquationErrorFunction(Nk,vGrid_k,mKPolicy,mLaborPolicy_1,mConsumptionPolicy_1,mConsumptionPolicy_2,Na,mGrid_a1a2,mProb_a1a2,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL)
    % Then let's use LINEAR INTERPOLATION for tomorrow's values
    global inputs;
    mMarginalUtilityTodayFsolve = inputs.mMarginalUtilityTodayFsolve;
    
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
    %         [v,ikPrime]=min(abs(kPrime - vGrid_k));

            consumption_1 = mConsumptionPolicy_1(ik,ia);
            consumption_2 = mConsumptionPolicy_2(ik,ia);
    %         labor_1 = mLaborPolicy_1(ik,ia);
    %         labor_2 = mLaborPolicy_2(ik,ia);
    %         laborInitial=[vLabor(1),vLabor(2)]; % update the initial guess for labor policy to speed up the process
    %         leftHandSideOfEulerEquation(ik,ia) = mmu_1 * (a_2 * labor_2)^mmu_2 * ...
    %             (a_1 * k^aalphaK * labor_1^aalphaL + (1-ddelta) * k - kPrime)^(mmu_1-1);

%             marginalUtilityToday = mmu_1 * (consumption_2)^mmu_2 * (consumption_1)^(mmu_1-1);
            marginalUtilityToday = interp1(vGrid_k, mMarginalUtilityTodayFsolve(ik,ia,:),kPrime);
            
            leftHandSideOfEulerEquation(ik,ia) = marginalUtilityToday;

    % Compute expected marginal utility of tomorrow
            expected = zeros(Na,1);

            for iaPrime = 1:Na
%                 aPrime = mGrid_a1a2(iaPrime,:);
                aPrime_1 = mGrid_a1a2(iaPrime,1);
                aPrime_2 = mGrid_a1a2(iaPrime,2); 

    %             kPrimePrime = mKPolicy(ikPrime,iaPrime);
    %             consumptionPrime_1 = mConsumptionPolicy_1(ikPrime,iaPrime);
    %             consumptionPrime_2 = mConsumptionPolicy_2(ikPrime,iaPrime);
    %             laborPrime_1 = mLaborPolicy_1(ikPrime,iaPrime);
%                 consumptionPrime_1 = interp1(vGrid_k,mConsumptionPolicy_1(:,iaPrime),kPrime);
%                 consumptionPrime_2 = interp1(vGrid_k,mConsumptionPolicy_2(:,iaPrime),kPrime);
%                 laborPrime_1 = interp1(vGrid_k,mLaborPolicy_1(:,iaPrime),kPrime);

                marginalUtilityTomorrow = interp1(vGrid_k, mMarginalUtilityTodayFsolve(ik,ia,:),kPrime);
%                 marginalUtilityTomorrow = mmu_1 * (consumptionPrime_2)^mmu_2 * (consumptionPrime_1)^(mmu_1-1);
                returnOnCapital = 1 - ddelta + aPrime_1 * laborPrime_1^aalphaL * aalphaK * (kPrime)^(aalphaK-1);

                unexpected = bbeta * marginalUtilityTomorrow * returnOnCapital;
                expected(iaPrime) = unexpected * mProb_a1a2(ia,iaPrime);
            end
            rightHandSideOfEulerEquation(ik,ia) = sum(expected);
        end
    end

    errorEulerEquation = abs(leftHandSideOfEulerEquation - rightHandSideOfEulerEquation);

end