function [utility] = utilityFunction(consumption_1,consumption_2,labor_1,labor_2,mmu_1,mmu_2)
utility = consumption_1.^mmu_1.* consumption_2.^mmu_2 -((labor_1+labor_2).^2)/2;
end