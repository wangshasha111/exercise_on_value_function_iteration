classdef valueFunctionClass < handle
    properties
        value = 0;
        labor_1 = 0;
        labor_2 = 0;
        consumption_1 = 0;
        consumption_2 = 0;
    end
   
    methods       
        function obj = valueFunction_setLaborToSteadyState(obj,kPrime,ik,k,ia,a_1,a_2,expectedValue0,bbeta,mmu_1,mmu_2,ddelta,aalphaK,aalphaL,labor_1_SteadyState,labor_2_SteadyState,value,labor_1,labor_2,consumption_1,consumption_2)
            
            global inputs;
            vGrid_a1 = inputs.vGrid_a1;
            vGrid_a2 = inputs.vGrid_a2;
            mProb_a1 = inputs.mProb_a1;
            mProb_a2 = inputs.mProb_a2;
            vGrid_k = inputs.vGrid_k;

            obj.labor_1 = obj.labor_1 + labor_1_SteadyState;
            obj.labor_2 = obj.labor_2 + labor_2_SteadyState;
            obj.consumption_1 = obj.consumption_1 + consumptionFunction1(a_1,k,kPrime,obj.labor_1,aalphaK,aalphaL,ddelta);
            obj.consumption_2 = obj.consumption_2 + consumptionFunction2(a_2,obj.labor_2);

            currentUtility = utilityFunction(obj.consumption_1,obj.consumption_2,obj.labor_1,obj.labor_2,mmu_1,mmu_2);

            if (obj.consumption_1 >= 0) && (obj.consumption_2 >= 0)  
                expectedValue = interp1(vGrid_k,expectedValue0(:,ia),kPrime);
                obj.value = (1-bbeta)*currentUtility + bbeta * expectedValue;
            else
                obj.value = -1e10;
            end

        end
      
        function
          
        end
    end
   
end

% Then in your current code create an object and the function handle, to the function on your object
% myObj = MyClass( );
% f = @myObj.myFunc