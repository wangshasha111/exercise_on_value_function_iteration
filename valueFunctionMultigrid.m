% This function is only different from valueFunction.m in using interp2 instead of
% interp1
function [value] = valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,vGrid_kMultigrid,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL)
% function [value,labor_1,labor_2,consumption_1,consumption_2] = valueFunction(kPrime,k,ik,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL)
% I want to do single dimension maximization, but I also want to retrieve
% other outputs under optimization fmincon or fminbnd or fminsearch. How to
% do that?
    global inputs;
    vGrid_a1 = inputs.vGrid_a1;
    vGrid_a2 = inputs.vGrid_a2;
    mProb_a1 = inputs.mProb_a1;
    mProb_a2 = inputs.mProb_a2;
    vGrid_k = inputs.vGrid_k;

    mLabor_1Fsolve = inputs.mLabor_1Fsolve;  % kPrime,k,[a_1,a_2]
    mLabor_2Fsolve = inputs.mLabor_2Fsolve;% kPrime,k,[a_1,a_2]
    mConsumption_1Fsolve = inputs.mConsumption_1Fsolve;% kPrime,k,[a_1,a_2]
    mConsumption_2Fsolve = inputs.mConsumption_2Fsolve;% kPrime,k,[a_1,a_2]
    mCurrentUtilityFsolve = inputs.mCurrentUtilityFsolve;% kPrime,k,[a_1,a_2]

    % do matrix transformations for linear interpolation
%     mLabor_1Fsolve = permute(mLabor_1Fsolve,[1,3,2]);% kPrime,k,[a_1,a_2]
%     mLabor_2Fsolve = permute(mLabor_2Fsolve,[1,3,2]);% kPrime,k,[a_1,a_2]
%     mConsumption_1Fsolve = permute(mConsumption_1Fsolve,[1,3,2]);% kPrime,k,[a_1,a_2]
%     mConsumption_2Fsolve = permute(mConsumption_2Fsolve,[1,3,2]);% kPrime,k,[a_1,a_2]
%     mCurrentUtilityFsolve = permute(mCurrentUtilityFsolve,[1,3,2]);% kPrime,k,[a_1,a_2]
    
    % use linear interpolation to find the optimized value
    labor_1 = interp2(vGrid_k,vGrid_k,mLabor_1Fsolve(:,:,ia),kPrime,k);
    labor_2 = interp2(vGrid_k,vGrid_k,mLabor_2Fsolve(:,:,ia),kPrime,k);
    consumption_1 = interp2(vGrid_k,vGrid_k,mConsumption_1Fsolve(:,:,ia),kPrime,k);
    consumption_2 = interp2(vGrid_k,vGrid_k,mConsumption_2Fsolve(:,:,ia),kPrime,k);
    currentUtility = interp2(vGrid_k,vGrid_k,mCurrentUtilityFsolve(:,:,ia),kPrime,k);

    if (consumption_1 >= 0) && (labor_1 >= 0) ...
            && (consumption_2 >= 0) && (labor_2 >= 0) 

    %     expectedValue  = expectedValue0(ikPrimeLow,ia) + ...
    %         (expectedValue0(ikPrimeHigh,ia) - expectedValue0(ikPrimeLow,ia))./...
    %         (vGrid_k(ikPrimeHigh) - vGrid_k(ikPrimeLow)).*(kPrime - vGrid_k(ikPrimeLow));

        expectedValue = interp1(vGrid_kMultigrid,expectedValue0(:,ia),kPrime);

        value = (1-bbeta)*currentUtility + bbeta * expectedValue;

    else
        value = -1e10;
    end
end
