function [consumption_1] = consumptionFunction1(a_1,k,kPrime,labor_1,aalphaK,aalphaL,ddelta)
consumption_1 = a_1 * k.^aalphaK.* labor_1.^aalphaL + (1-ddelta) * k - kPrime;
end