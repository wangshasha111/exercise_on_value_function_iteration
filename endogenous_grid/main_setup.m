% Homework 1 Econ714 JFV
% A two good production economy
% Written by Shasha Wang

close all;
clear;

% cd 'E:\Dropbox\fall 19-20\jesus\homework_1\set_labor_to_steady_state_first_then_standard_iteration';

%% Calibration
bbeta = 0.96;

mmu_1 = 0.5; % utility weight for good 1
mmu_2 = 0.5; % utility weight for good 2

ddelta = 0.1; % depreciation for capital used to produce good 1

aalphaK = 0.33;
aalphaL = 0.67;

%% Shocks
vGrid_a1 = exp([-0.0673; -0.0336; 0; 0.0336; 0.0673]);
Na_1 = length(vGrid_a1);
mProb_a1 = [0.9727 0.0273 0 0 0;...
            0.0041 0.9806 0.0153 0 0;...
            0 0.0082 0.9836 0.0082 0;...
            0 0 0.0153 0.9806 0.0041;...
            0 0 0 0.0273 0.9727];
        
vGrid_a2 = [0.9;1;1.1];
Na_2 = length(vGrid_a2);
mProb_a2 = [0.9 0.1 0;...
            0.05 0.9 0.05;...
            0 0.1 0.9];

% Combine the two shocks into one shock
mProb_a1a2 = kron(mProb_a1,mProb_a2);
[A1,A2] = meshgrid(vGrid_a2,vGrid_a1);
temp=cat(2,A2',A1');
mGrid_a1a2=reshape(temp,[],2);
Na = Na_1 * Na_2;

global inputs;
inputs.vGrid_a1 = vGrid_a1;
inputs.vGrid_a2 = vGrid_a2;
inputs.mProb_a1 = mProb_a1;
inputs.mProb_a2 = mProb_a2;
inputs.mGrid_a1a2 = mGrid_a1a2;
inputs.mProb_a1a2 = mProb_a1a2;

