function [policyFunction,consumptionFunction,laborFunction ,valueFunction,eulerError,maxError]= ...
         eulerr_grid(theta,inputs)

bbeta  = inputs.params(1);
ddelta = inputs.params(2);
ppsi   = inputs.params(3);
rrho   = inputs.params(4);
eeta   = inputs.params(5);

gridProd    = inputs.gridProd;
gridAlpha   = inputs.gridAlpha;
nodeNum     = inputs.simulParams(1);
sizeProd    = inputs.simulParams(2);
sizeAlpha   = inputs.simulParams(3);
capitalMin  = inputs.simulParams(6);
capitalMax  = inputs.simulParams(7);
M           = nodeNum*sizeProd*sizeAlpha;
sizeShocks  = sizeProd*sizeAlpha;

transitionMatrix = inputs.transitionMatrix;

chebyPoly   = inputs.chebyPoly;
     
eulerSize       = inputs.eulerSize;
interval        = inputs.interval;
gridEuler = zeros(eulerSize,1);

for i = 1:eulerSize
    gridEuler(i) = capitalMin+(i-1)*interval/(eulerSize-1);
end
    
% This function computes the Euler Errors on the capital and exogenous shock grid

gridEulerScaled = (2*gridEuler-(capitalMin+capitalMax))/(capitalMax-capitalMin);

chebyPolyEuler      = ones(eulerSize,nodeNum);
chebyPolyEuler(:,2) = gridEulerScaled;

for i1 = 3:nodeNum
    chebyPolyEuler(:,i1) = 2*gridEulerScaled.*chebyPolyEuler(:,i1-1)-chebyPolyEuler(:,i1-2);
end     

eulerError          = zeros(eulerSize,sizeShocks);
valueFunction       = zeros(eulerSize,sizeShocks);
laborFunction       = zeros(eulerSize,sizeShocks);
consumptionFunction = zeros(eulerSize,sizeShocks);
policyFunction      = zeros(eulerSize,sizeShocks);

for indexShocks = 1:sizeShocks
    
    indexProd  = ceil(indexShocks/sizeAlpha);
    indexAlpha = mod(indexShocks-1,sizeAlpha)+1;
    aalpha     = gridAlpha(indexAlpha);
    z          = gridProd(indexProd);

    for indexCapital = 1:eulerSize % Loop 1 over collocation point on k
        
        capital    = gridEuler(indexCapital);

        thetaVFunctionShock = theta(((indexShocks-1)*nodeNum+1):indexShocks*nodeNum,1);
        thetaLaborShock     = theta(((indexShocks-1)*nodeNum+1):indexShocks*nodeNum,2);
        valueFunction(indexCapital,indexShocks) = dot(thetaVFunctionShock,chebyPolyEuler(indexCapital,:)); % Value function at each collocation points
        laborFunction(indexCapital,indexShocks) = dot(thetaLaborShock,chebyPolyEuler(indexCapital,:));       % Labor at each collocation points

        labor      = laborFunction(indexCapital,indexShocks);
        
        output = exp(z)*capital^aalpha*labor^(1-aalpha);
        consumptionFunction(indexCapital,indexShocks) = (1-aalpha)*exp(z)*capital^aalpha*labor^(-aalpha)/(eeta*labor);            
        policyFunction(indexCapital,indexShocks) = output+(1-ddelta)*capital-consumptionFunction(indexCapital,indexShocks);            

    end % Loop 1 over collocation point on k ends

    % Scale k prime from [k_min,k_max] to [-1,1]
    policyFunctionScaledDown = (2*policyFunction(:,indexShocks)-(capitalMin+capitalMax))/(capitalMax-capitalMin);

    % value of polynomials at each scaled k prime
    chebyPolyEulerPrime = ones(eulerSize,nodeNum);
    chebyPolyEulerPrime(:,2) = policyFunctionScaledDown;

    for i1 = 3:nodeNum
        chebyPolyEulerPrime(:,i1) = 2*policyFunctionScaledDown.*chebyPolyEulerPrime(:,i1-1)-chebyPolyEulerPrime(:,i1-2);
    end     

    % Calculate residual        
    for indexCapital = 1:eulerSize % Loop 2 over collocation point on k  
        
        capitalPrime      = policyFunction(indexCapital,indexShocks);

        continuationValue = zeros(sizeShocks,1);
        temp              = zeros(sizeShocks,1);

        for indexShocksPrime = 1:sizeShocks
            
            indexProdPrime  = ceil(indexShocksPrime/sizeAlpha);
            indexAlphaPrime = mod(indexShocksPrime-1,sizeAlpha)+1;
            aalphaPrime     = gridAlpha(indexAlphaPrime);
            zPrime          = gridProd(indexProdPrime);

            thetaVFunctionShock = theta(((indexShocksPrime-1)*nodeNum+1):indexShocksPrime*nodeNum,1);
            thetaLaborShock = theta(((indexShocksPrime-1)*nodeNum+1):indexShocksPrime*nodeNum,2);
            continuationValue(indexShocksPrime) = dot(thetaVFunctionShock,chebyPolyEulerPrime(indexCapital,:));
            laborPrime = dot(thetaLaborShock,chebyPolyEulerPrime(indexCapital,:));

            consumptionPrime   = ((1 - aalphaPrime) * exp(zPrime)*capitalPrime^aalphaPrime*laborPrime^(-aalphaPrime))/(eeta*laborPrime);

            consPrimeMarginalUtil  = (log(consumptionPrime) - eeta * laborPrime^2/2)^(rrho-1)/consumptionPrime;
            capPrimeMarginalProd   = aalphaPrime*exp(zPrime)*capitalPrime^(aalphaPrime-1)*laborPrime^(1-aalphaPrime);
            temp(indexShocksPrime) = continuationValue(indexShocksPrime)^(ppsi-rrho)*consPrimeMarginalUtil*(capPrimeMarginalProd+(1-ddelta));

        end

        eulerRHS = bbeta*dot(transitionMatrix(indexShocks,:),continuationValue.^ppsi)^(rrho/ppsi-1)*dot(transitionMatrix(indexShocks,:),temp);

        labor       = laborFunction(indexCapital,indexShocks);
        consumption = consumptionFunction(indexCapital,indexShocks);

        eulerLHS   = (log(consumption) - eeta * labor^2/2)^(rrho-1)/consumption;

        eulerError(indexCapital,indexShocks) = eulerLHS - eulerRHS;
    end % Loop 2 over k ends

end % Loop over z ends

eulerError = log10(abs(eulerError));
maxError   = max(eulerError,[],1);