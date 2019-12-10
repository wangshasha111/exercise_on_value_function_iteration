% Homework 2
% ECON 714 - Prof. Jesus Fernandez-Villaverde
% Projection
% Sara Casella
% December 2018

clear
%clc
close all

%% Calibration

% Parameters
theta = 0.5;
gamma = 10;
beta  = 0.99;
delta = 0.1;

% save parameters in an array
parameters.theta = theta;
parameters.gamma = gamma;
parameters.beta  = beta;
parameters.delta = delta;

% Shocks processes

vProductivity = exp([-0.0673 -0.0336 0 0.0336 0.0673]);
nGridProductivity = length(vProductivity);

mTransitionProductivity  = [0.9727 0.0273 0      0      0;
                            0.0041 0.9806 0.0153 0      0; 
                            0      0.0082 0.9836 0.0082 0;	
                            0      0      0.0153 0.9806 0.0041;
                            0      0      0      0.0273 0.9727];
                             
vCapitalShare = [0.25 0.3 0.35];
nGridCapitalShare= length(vCapitalShare);

mTransitionCapitalShare =  [0.9  0.07 0.03;
                            0.05 0.9 0.05;
                            0.03 0.07 0.9];
                        
% Collapse the two shocks into one made of all the possible combinations
% 5*3 = 15

vShocks = zeros(2, length(vProductivity)*length(vCapitalShare));
nGridShocks = length(vShocks);

for i=1:length(vProductivity)
    for j=1:length(vCapitalShare)
        
      vShocks(:,j+(i-1)*length(vCapitalShare))=...
          [vProductivity(i), vCapitalShare(j)]';
                                                   
    end 
end

mTransition = kron(mTransitionProductivity,mTransitionCapitalShare); 

%% Deterministic Steady State

z_ss = 1;        % productivity
alpha_ss = 0.3;  % share of capital

l_ss = 100;       % labor normalized to 10
k_ss = (beta*alpha_ss*z_ss*l_ss^(1-alpha_ss)/(1+beta*(delta-1)))^(1/(1-alpha_ss));
c_ss = -delta*k_ss + z_ss*k_ss^alpha_ss*l_ss^(1-alpha_ss);
y_ss = z_ss*k_ss^alpha_ss*l_ss^(1-alpha_ss);

% Normalization in the utility function
eta = (1-alpha_ss)*z_ss*k_ss^alpha_ss*l_ss^(-1-alpha_ss)/c_ss;

% Steady state value. Note we scale up utility by (1-beta)
v_ss = (log(c_ss)-eta*l_ss^2/2);

parameters.eta = eta;

%% Min and Max values for capital and labor

% Capital
minCapital = 0.7*k_ss;
maxCapital = 1.5*k_ss;

% Labor
getLabor = @(z,k,kprime,alpha,lowerBound,upperBound) fminbnd(@(l)...
     ((1-alpha)*z*k^alpha*l.^(-(1+alpha))-eta*z*k^alpha*l.^(1-alpha)-eta*(1-delta)*k+eta*kprime).^2,...
     lowerBound, upperBound); % from mkt clearing

maxLabor = zeros(1,nGridShocks);
minLabor = zeros(1,nGridShocks);

for iShocks = 1:nGridShocks
    productivity = vShocks(1,iShocks);
    alpha        = vShocks(2,iShocks);
    minLabor(iShocks) = getLabor(productivity,maxCapital,minCapital,alpha,0.0001,5*l_ss);
    maxLabor(iShocks) = getLabor(productivity,minCapital,maxCapital,alpha,0.0001,5*l_ss);
end

%% Declare size vectors

% Productivity and capital
shock_num = nGridShocks;           % number of nodes for technology process Z
nodeSizeList = [3,8];              % number of polynomials for capital
listSize  = length(nodeSizeList);  % steps in multinode iteration for capital

% Euler errors
grid_num  = 3000;                  % # of grid points for  capital (to compute euler errors)


%% Computation of the projection coefficients

tic
for multinode_step=1:listSize
    
    node_num = nodeSizeList(multinode_step);
    M = node_num*shock_num;
    
    % Find Zeros of the Chebychev Polynomial on order M 
    chebZeros=cos((2*(1:node_num)'-1)*pi/(2*node_num));

    % Define Chebychev polynomials recursively
    chebPols=ones(node_num,node_num);
    chebPols(:,2)=chebZeros;

    for i=3:node_num
        chebPols(:,i)=2*chebZeros.*chebPols(:,i-1)-chebPols(:,i-2);
    end

    % Project collocation points in the capital space
    capitalGridColloc=((chebZeros+1)*(maxCapital-minCapital))/2+minCapital;

    % Initial Guess for Chebyshev coefficients
    coefGuess = zeros(2*M,1);
    if(multinode_step == 1)
        for shockIndex = 1:shock_num
            coefGuess((shockIndex-1)*node_num+1) = v_ss;
            coefGuess((shockIndex-1)*node_num+M+1) = l_ss;
        end
    else
        for shockIndex = 1:2*shock_num
            coefGuess((shockIndex-1)*node_num+1: (shockIndex-1)*node_num+node_num_old)...
            = coefGuess_old((shockIndex-1)*node_num_old+1: shockIndex*node_num_old);
        end
    end
    
    % Solve for Chebyshev coefficients
    projectionCoefficients = getChebychevCoefs(parameters,minCapital,maxCapital,...
        minLabor,maxLabor,coefGuess,capitalGridColloc,chebPols,vShocks,mTransition,...
        node_num,shock_num,M);
    coefGuess_old = projectionCoefficients;
    node_num_old = node_num;
end
toc
computationTimeCheb=toc;

%% Results

nGridComplete=grid_num;
GridStep = (maxCapital - minCapital)/(nGridComplete-1);
kGridComplete = (minCapital: GridStep : maxCapital)';

[mPolicyCapital,mPolicyConsumption,mPolicyLabor,mValueFunction,EulerErrorCheb,MaxErrorCheb]= ...
                    getEulerErrors(parameters,...
                    projectionCoefficients,vShocks,mTransition,...
                    minCapital,maxCapital,kGridComplete, nGridShocks,...
                    node_num,grid_num,M);
                
meanErrorCheb=mean(mean(EulerErrorCheb));
maxErrorCheb=max(max(EulerErrorCheb));              

save ChebyshevProjection.mat;

%% Plot 
% All together
figure
subplot(2,1,1)
plot(kGridComplete,mValueFunction)
title('Value function Chebyshev')
xlabel('k')
ylabel('V')
xlim([minCapital-10, maxCapital+10])
subplot(2,1,2)
h1 = plot(kGridComplete,mPolicyCapital);
hold on
h2 = plot(kGridComplete,kGridComplete);
title('Capital Decision Rule Chebyshev')
xlabel('k')
ylabel('k prime')
legend(h2,{'45° line'},'location','northwest')
xlim([minCapital-10, maxCapital+10])

figure
plot(kGridComplete,EulerErrorCheb)
xla=xlabel('Capital');
yla=ylabel('Log10(Euler Errors)');
tit=title('Location of the Errors on the Grid Chebyshev');

avgLogEulerErrorsCheb=mean(EulerErrorCheb,2);

totAvgLogEulerErrorsCheb=mean(avgLogEulerErrorsCheb);

maxEulerErrorsCheb=max(max(EulerErrorCheb));


