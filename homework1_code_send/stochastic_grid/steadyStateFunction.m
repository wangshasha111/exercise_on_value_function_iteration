function [y] = steady_state_function(input_ss,bbeta,ddelta,aalphaK,aalphaL,mmu_1,mmu_2) % input is k,labor_1,labor_2
%     y(1) = input_ss(1)/input_ss(2) - ((1/bbeta-0.9)/aalphaK)^(-1/aalphaL);
    y(1) = 1-bbeta*(aalphaK * input_ss(1)^(aalphaK-1) * input_ss(2)^aalphaL + 1 - ddelta);
%     y(2) = aalphaL*input_ss(3) - input_ss(2)^aalphaL * (1-0.1/((1/bbeta-0.9)/aalphaK));
    y(2) = aalphaL*input_ss(3)* input_ss(2)^(aalphaL-1) - input_ss(2)^ aalphaL + 0.1 * input_ss(1)^(1-aalphaK); 
%     y(3) = input_ss(2) + input_ss(3) - mmu_2 * (input_ss(2)/input_ss(3))^0.5 * (((1/bbeta-0.9)/aalphaK)^(-aalphaK/aalphaL)-0.1*((1/bbeta-0.9)/aalphaK)^(-1/aalphaL))^0.5;
    y(3) = input_ss(2) + input_ss(3) - mmu_2 * input_ss(3)^(mmu_2-1)*(input_ss(1)^aalphaK * input_ss(2)^aalphaL-ddelta * input_ss(1))^mmu_1;
end