function [value,labor_1,labor_2,consumption_1,consumption_2] = valueFunction_stochasticGrid(kPrime,ikPrime,ik,k,ia,a_1,a_2,iik,iia,iikPrime,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL,laborInitial)

    global inputs;
    vGrid_a1 = inputs.vGrid_a1;
    vGrid_a2 = inputs.vGrid_a2;
    mProb_a1 = inputs.mProb_a1;
    mProb_a2 = inputs.mProb_a2;
    vGrid_k = inputs.vGrid_k;
    
    mLabor_1Fsolve = inputs.mLabor_1Fsolve;
    mLabor_2Fsolve = inputs.mLabor_2Fsolve;
    mConsumption_1Fsolve = inputs.mConsumption_1Fsolve;
    mConsumption_2Fsolve = inputs.mConsumption_2Fsolve;
    mCurrentUtilityFsolve = inputs.mCurrentUtilityFsolve;	
    % laborFunction = inputs.laborFunction;

%     laborInitial=[0.2,0.2];
%    opts1 = optimoptions('fsolve','Tolx',1e-6, 'Display','off');
%    vLabor = fsolve(@(labor) laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta), laborInitial,opts1);
%    labor_1 = vLabor(1);
%    labor_2 = vLabor(2);
    


%    	consumption_1 = consumptionFunction1(a_1,k,kPrime,labor_1,aalphaK,aalphaL,ddelta);
%    	consumption_2 = consumptionFunction2(a_2,labor_2);
%    	currentUtility = utilityFunction(consumption_1,consumption_2,labor_1,labor_2,mmu_1,mmu_2);
	currentUtility = mCurrentUtilityFsolve(ikPrime,ik,ia);
    labor_1 = mLabor_1Fsolve(ikPrime,ik,ia);
    labor_2 = mLabor_2Fsolve(ikPrime,ik,ia);    
   	consumption_1 = mConsumption_1Fsolve(ikPrime,ik,ia);
   	consumption_2 = mConsumption_2Fsolve(ikPrime,ik,ia);

        expectedValue = expectedValue0(iikPrime,iia);

        value = (1-bbeta)*currentUtility + bbeta * expectedValue;

end
