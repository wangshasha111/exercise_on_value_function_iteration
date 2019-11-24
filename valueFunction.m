function [value] = valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL)
% function [value,labor_1,labor_2,consumption_1,consumption_2] = valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL)

    global inputs;
    vGrid_a1 = inputs.vGrid_a1;
    vGrid_a2 = inputs.vGrid_a2;
    mProb_a1 = inputs.mProb_a1;
    mProb_a2 = inputs.mProb_a2;
    vGrid_k = inputs.vGrid_k;

    mLabor_1Fsolve = inputs.mLabor_1Fsolve; % k,[a_1,a_2],kPrime
    mLabor_2Fsolve = inputs.mLabor_2Fsolve;
    mConsumption_1Fsolve = inputs.mConsumption_1Fsolve;
    mConsumption_2Fsolve = inputs.mConsumption_2Fsolve;
    mCurrentUtilityFsolve = inputs.mCurrentUtilityFsolve;

% laborInitial=[0.2,0.2];
% opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','off');
% vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
% labor_1 = vLabor(1);
% labor_2 = vLabor(2);
% consumption_1 = consumptionFunction1(a_1,k,kPrime,labor_1,aalphaK,aalphaL,ddelta);
% consumption_2 = consumptionFunction2(a_2,labor_2);
% currentUtility = utilityFunction(consumption_1,consumption_2,labor_1,labor_2,mmu_1,mmu_2);

% find the index of kPrime, a1, a2, k
% ikPrimeLow  = max(sum(kPrime > vGrid_k),1);
% ikPrimeHigh = ikPrimeLow + 1;
% ia_1 = sum(a_1 >= vGrid_a1);
% ia_2 = sum(a_2 >= vGrid_a2);
% ik = sum(k >= vGrid_k);
    labor_1 = interp1(vGrid_k,mLabor_1Fsolve(:,ia,ik),kPrime);
    labor_2 = interp1(vGrid_k,mLabor_2Fsolve(:,ia,ik),kPrime);
    consumption_1 = interp1(vGrid_k,mConsumption_1Fsolve(:,ia,ik),kPrime);
    consumption_2 = interp1(vGrid_k,mConsumption_2Fsolve(:,ia,ik),kPrime);

    currentUtility = interp1(vGrid_k,mCurrentUtilityFsolve(:,ia,ik),kPrime);

    if (consumption_1 >= 0) && (labor_1 >= 0) ...
            && (consumption_2 >= 0) && (labor_2 >= 0) 

    %     expectedValue  = expectedValue0(ikPrimeLow,ia) + ...
    %         (expectedValue0(ikPrimeHigh,ia) - expectedValue0(ikPrimeLow,ia))./...
    %         (vGrid_k(ikPrimeHigh) - vGrid_k(ikPrimeLow)).*(kPrime - vGrid_k(ikPrimeLow));

        expectedValue = interp1(vGrid_k,expectedValue0(:,ia),kPrime);

        value = (1-bbeta)*currentUtility + bbeta * expectedValue;

    else
        value = -1e10;
    end
end
