function [mValueFunction, mPolicyCapital, mPolicyLabor, mPolicyConsumption, error] = VFIinterpolation(...
mValueFunctionGuess,vGridCapital,vShocks,mTransition,minCapital,maxCapital,laborMin,...
laborMax,functions,parameters, tolerance, maxIterations)

gamma = parameters.gamma;

nGridCapital = length(vGridCapital);
nGridShocks = length(vShocks);

%Initial guess for value function: previously found guess
mValueFunction = mValueFunctionGuess; 

% Preallocation for required matrices
mValueFunctionNew     = zeros(nGridCapital,nGridShocks);
mPolicyCapital        = zeros(nGridCapital,nGridShocks);
mPolicyLabor          = zeros(nGridCapital,nGridShocks);
mPolicyConsumption    = zeros(nGridCapital,nGridShocks);

% Tolerance and max iterations
error         = 1;
iteration     = 1;


while error > tolerance && iteration <= maxIterations
    
    % Continuation Value with Epstein Zin Preferences
     expectedValueFunction = mValueFunction.^(1-gamma)*mTransition';

    for iShocks = 1:nGridShocks
        
        productivity = vShocks(1,iShocks);
        alpha        = vShocks(2,iShocks);
                
         for iCapital = 1:nGridCapital
             
             capitalToday = vGridCapital(iCapital);

             interk = vGridCapital;
             interv = expectedValueFunction(:,iShocks)';

             capitalChoice = fminbnd(@(kprime) (-getValue(productivity,...
                alpha,capitalToday,kprime,parameters,functions,...
                laborMin(iShocks),laborMax(iShocks),...
                interk,interv)),minCapital,maxCapital);
             [value, laborChoice, consumptionChoice] = getValue(productivity,...
                alpha,capitalToday,capitalChoice,...
                parameters,functions,...
                laborMin(iShocks),laborMax(iShocks),...
                interk,interv);
                
             mValueFunctionNew(iCapital,iShocks) = value;
             mPolicyCapital(iCapital,iShocks)    = capitalChoice;
             mPolicyLabor(iCapital,iShocks)      = laborChoice;
             mPolicyConsumption(iCapital,iShocks) = consumptionChoice;
        
         end   
        
    end
    
    error = max(abs(mValueFunctionNew(:)-mValueFunction(:)));
    mValueFunction = mValueFunctionNew;
    
%    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration VFI = %d, Sup Diff = %2.8f\n', iteration, error); 
%    end
    
    iteration = iteration+1;
 

end
