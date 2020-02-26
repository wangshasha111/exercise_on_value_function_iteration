%% EGM 2


function [mValueFunctionEGM, mEndogenousCapital, mEndogenousLabour] = EGM(parameters,vGridCapital,vShocks,...
    mTransition,c_ss,k_ss,mPolicyCapital,mPolicyLabor,mValueFunctionGuess,tolerance,options)

gamma = parameters.gamma;
theta = parameters.theta;
eta = parameters.eta;
beta = parameters.beta;
delta = parameters.delta;

nGridCapital = length(vGridCapital);
nGridShocks = length(vShocks);

% Step 3.a in paper: Use an interpolation method and solve for all 
% k_endogenous such that k'(k_endogenous,z,alpha)=Gk_t+1, we use the policy
% function for capital from VFI.

mEndogenousCapital = zeros(nGridCapital,nGridShocks);

for iShocks = 1:nGridShocks

            for iCapital = 1:nGridCapital
                capitalLow = max(sum(vGridCapital(iCapital)>mPolicyCapital(:,iShocks)),1);
                 if capitalLow == nGridCapital
                    capitalLow = nGridCapital-1;
                 end
                 capitalHigh = capitalLow + 1;
                 
                 while (mPolicyCapital(capitalHigh,iShocks)-mPolicyCapital(capitalLow,iShocks) == 0)
                    if capitalLow ~= nGridCapital-1
                    capitalHigh = capitalHigh+1;
                    else
                    capitalLow = capitalLow-1;
                    end
                 end
                 
                 
                 mEndogenousCapital(iCapital,iShocks) = vGridCapital(capitalLow)+...
                     (vGridCapital(iCapital)-mPolicyCapital(capitalLow,iShocks))*...
                     (vGridCapital(capitalHigh)-vGridCapital(capitalLow))/(mPolicyCapital(capitalHigh,iShocks)...
                      -mPolicyCapital(capitalLow,iShocks));
            end

end


% Step 3.b: Use an interpolation to define the Labour Policy obtained in
% the VFI on the K_endogenous grid defined in Step 3.a

mEndogenousLabour = zeros(nGridCapital,nGridShocks);

for iShocks = 1:nGridShocks
    for iCapital = 1:nGridCapital

               capitalLow = max(sum(mEndogenousCapital(iCapital,iShocks)>vGridCapital),1);
               
          
                if capitalLow == nGridCapital
                capitalLow = nGridCapital-1;
                end
                
              capitalHigh = capitalLow+1;
              
              mEndogenousLabour(iCapital,iShocks) = mPolicyLabor(capitalLow,iShocks)...
                  +(mEndogenousCapital(iCapital,iShocks)-vGridCapital(capitalLow))*...
                   (mPolicyLabor(capitalHigh,iShocks)-mPolicyLabor(capitalLow,iShocks))...
                   /(vGridCapital(capitalHigh)-vGridCapital(capitalLow));
    end
    
end


%save('mEndogenousLabour.mat','mEndogenousLabour');

% Step 3.C: Do EGM

% Initial guess for V tilde.
mValueFunctionTilde2 = beta*(mValueFunctionGuess.^(1-gamma)*mTransition').^((1-theta)/(1-gamma));

maxDifference3 = 1;
iteration3 = 0;
%tolerance3 = 0.1*tolerance;
tolerance3=0.1*tolerance;
maxIterations = 1000;

mEndogenousValueFunction2   = zeros(nGridCapital,nGridShocks);
mEndogenousCapital2         = zeros(nGridCapital,nGridShocks);
mValueFunctionEGM      = zeros(nGridCapital,nGridShocks);
%mValueFunctionTildeUpdated2 = zeros(nGridCapital,nGridShocks);

while (maxDifference3>tolerance3)  && iteration3<maxIterations
    for iShocks = 1:nGridShocks
        productivity = vShocks(1,iShocks);
        alpha        = vShocks(2,iShocks);
  
            for iCapital = 1:nGridCapital  
                
% Step 3.d: Do the derivative and consumption calculations   

                % Taking the derivative of ValueFunctionTilde
                if iCapital == nGridCapital
                    
                    derivativeValueFunctionTilde2 = (mValueFunctionTilde2(iCapital,iShocks)-...
                    mValueFunctionTilde2(iCapital-1,iShocks))/(vGridCapital(iCapital)-vGridCapital(iCapital-1));
                
                elseif iCapital == 1
                    
                    derivativeValueFunctionTilde2 = (mValueFunctionTilde2(iCapital+1,iShocks)-...
                    mValueFunctionTilde2(iCapital,iShocks))/(vGridCapital(iCapital+1)-vGridCapital(iCapital));
                                                   
                else
                    
                    derivativeValueFunctionTilde2 = ((mValueFunctionTilde2(iCapital,iShocks)-mValueFunctionTilde2(iCapital-1,iShocks))/...
                                                    (vGridCapital(iCapital)-vGridCapital(iCapital-1))+...
                                                    (mValueFunctionTilde2(iCapital+1,iShocks)-mValueFunctionTilde2(iCapital,iShocks))/...
                                                    (vGridCapital(iCapital+1)-vGridCapital(iCapital)))/2;
                end
                
                

                % Consumption
                ConsumptionFunction = @(x) (x-theta*(1-beta)*(log(x)-eta*(mEndogenousLabour(iCapital,iShocks)^2/2))^(theta-1)...
                    *derivativeValueFunctionTilde2^(-1))^2;
                      
               %Consumption = fzero(ConsumptionFunction, [exp(eta*(mEndogenousLabour(iCapital,iShocks)^2/2))+1, 1000]); 
               %Consumption = fsolve(ConsumptionFunction,exp(eta*(mEndogenousLabour(iCapital,iShocks)^2/2))+0.1,options);

%                Consumption = fsolve(ConsumptionFunction,2*c_ss,options);
                Consumption = fminbnd(ConsumptionFunction,exp(eta*(mEndogenousLabour(iCapital,iShocks)^2/2))+0.1,500);
                
                % Calculate the implied capital
                
                Labour = mEndogenousLabour(iCapital,iShocks);
    
                CapitalFunction =  @(x) (Consumption + vGridCapital(iCapital)...
                   -productivity*x^alpha*Labour^(1-alpha)-(1-delta)*x)^2;
%                Capital = fsolve(CapitalFunction,2*k_ss,options);
%                Capital = fsolve(CapitalFunction,k_ss,options);
                Capital = fminbnd(CapitalFunction,0.5*k_ss,2*k_ss);
                mEndogenousCapital2(iCapital,iShocks)=Capital;
                
% Step 3.e Update Value Function

                if (Consumption < 0 || isreal(((1-beta)*(log(Consumption) - eta*(Labour^2)/2)^(1-theta) + mValueFunctionTilde2(iCapital,iShocks))^(1/(1-theta)))==0)
                    disp('warning consumption or utility < 0')

                    break;
                else

                mEndogenousValueFunction2(iCapital,iShocks) =...
                 ((1-beta)*(log(Consumption) - eta*(Labour^2)/2)^(theta) ...
                 + mValueFunctionTilde2(iCapital,iShocks))^(1/theta);
                end
            end
            
            
                % Redefine it on the grid for capital
                
                for nCapital2 = 1:nGridCapital

                    CapitalLow = max(sum(vGridCapital(nCapital2)>mEndogenousCapital2(:,iShocks)));
                    if CapitalLow == 0
                        CapitalLow =1;
                    elseif CapitalLow == nGridCapital
                        CapitalLow = nGridCapital-1;
                    end
                        
                    CapitalHigh = CapitalLow + 1;
                    
                    mValueFunctionEGM(nCapital2,iShocks) = mEndogenousValueFunction2(CapitalLow,iShocks)...
                        +(vGridCapital(nCapital2)-mEndogenousCapital2(CapitalLow,iShocks))...
                        *(mEndogenousValueFunction2(CapitalHigh,iShocks)-mEndogenousValueFunction2(CapitalLow,iShocks))...
                        /(mEndogenousCapital2(CapitalHigh,iShocks)-mEndogenousCapital2(CapitalLow,iShocks));
                end
                
% Step 3.f Update Policy Function for Labour

              for nCapital3 = 1:nGridCapital

                CapitalLow = max(sum(mEndogenousCapital2(nCapital3,iShocks)>vGridCapital));
                if CapitalLow == 0
                    CapitalLow = 1;
                elseif CapitalLow == nGridCapital
                    CapitalLow = nGridCapital-1;
                end
                    

                CapitalHigh = CapitalLow+1;
                mEndogenousLabour(nCapital3,iShocks) = mPolicyLabor(CapitalLow,iShocks)+...
                (mEndogenousCapital2(nCapital3,iShocks)-vGridCapital(CapitalLow))*...
                (mPolicyLabor(CapitalHigh,iShocks)-mPolicyLabor(CapitalLow,iShocks))...
                /(vGridCapital(CapitalHigh)-vGridCapital(CapitalLow)); 
              end
    end
    
    if (isreal(beta*((mValueFunctionEGM.^(1-gamma))* mTransition').^(theta/(1-gamma))) == 0)
        
        mValueFunctionTildeUpdated2 = mValueFunctionTilde2;
        iteration3 = iteration3+1;
    else
    
    % Calculate Vtilde in t+1
    mValueFunctionTildeUpdated2 = beta*(mValueFunctionEGM.^(1-gamma)*mTransition').^(theta/(1-gamma));
    
    % Criteria for convergence
    maxDifference3 = max(abs(mValueFunctionTildeUpdated2(:)-mValueFunctionTilde2(:)));
    mValueFunctionTilde2 = mValueFunctionTildeUpdated2;
    
    iteration3 = iteration3+1;
    if (mod(iteration3,10)==0 || iteration3 ==1)
        fprintf(' IterationEGM = %d, Sup Diff = %2.8f\n', iteration3, maxDifference3); 
    end
    end
end
