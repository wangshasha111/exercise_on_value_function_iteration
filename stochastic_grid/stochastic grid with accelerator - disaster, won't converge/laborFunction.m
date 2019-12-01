function [y] = laborFunction(labor,a_1,a_2,k,kPrime,mmu_1,mmu_2,aalphaK,aalphaL,ddelta)
% We can see from the question that once given a1,a2,k,and chosen k',
% labor_1, labor_2 are both chosen by optimization, and consumption_1,
% consumption_2 are both chosen immediately.
% labor has two entries.
% y also has two entries.
y(1) = labor(1) + labor(2) - (a_2*labor(2))^mmu_2*mmu_1*(a_1*k^aalphaK*labor(1)^aalphaL + (1-ddelta)*k-kPrime)^(mmu_1-1)*a_1*k^aalphaK*aalphaL*labor(1)^(aalphaL-1);
y(2) = labor(1) + labor(2) - (a_1*k^aalphaK*labor(1)^aalphaL + (1-ddelta)*k-kPrime)^(mmu_1)*a_2^mmu_2*mmu_2*labor(2)^(mmu_2-1);
end

    