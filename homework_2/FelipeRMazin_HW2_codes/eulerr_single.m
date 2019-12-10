function [value,capitalPrime,consumption,labor,output,rb,rk,rkCond,eulerError] = ... 
        eulerr_single(theta,inputs,indexShocks,capital)

    bbeta       = inputs.params(1);
    ddelta      = inputs.params(2);
    ppsi        = inputs.params(3);
    rrho        = inputs.params(4);
    eeta        = inputs.params(5);
    
    gridCapital = inputs.gridCapital;
    gridProd    = inputs.gridProd;
    gridAlpha   = inputs.gridAlpha;
    nodeNum     = inputs.simulParams(1);
    sizeProd    = inputs.simulParams(2);
    sizeAlpha   = inputs.simulParams(3);
    capitalMin  = inputs.simulParams(6);
    capitalMax  = inputs.simulParams(7);
    sizeShocks  = sizeProd*sizeAlpha;
    
    transitionMatrix = inputs.transitionMatrix;

    indexProd  = ceil(indexShocks/sizeAlpha);
    indexAlpha = mod(indexShocks-1,sizeAlpha)+1;
    
    z = gridProd(indexProd);
    aalpha = gridAlpha(indexAlpha);

    capitalScaled = 2*(capital-capitalMin)/(capitalMax - capitalMin) - 1;

    chebyPoly = zeros(nodeNum,1);
    for i = 1:nodeNum % ith polynomial
        chebyPoly(i) = cos(real(i-1)*acos(capitalScaled));
    end
    
    thetaVFunctionShock = theta(((indexShocks-1)*nodeNum+1):indexShocks*nodeNum,1); % coefficient for value function
    thetaLaborShock     = theta(((indexShocks-1)*nodeNum+1):indexShocks*nodeNum,2); % coefficient for labor 

    value = dot(thetaVFunctionShock,chebyPoly);   % Value function at each collocation points    
    labor = dot(thetaLaborShock,chebyPoly);       % Labor at each collocation points

    output       = exp(z)*capital^aalpha*labor^(1-aalpha);
    consumption  = (1-aalpha)*exp(z)*capital^aalpha*labor^(-aalpha)/(eeta*labor);            
    capitalPrime = output+(1-ddelta)*capital-consumption;

    capitalPrimeScaled = 2*(capitalPrime-capitalMin)/(capitalMax - capitalMin) - 1;  
    chebyPolyPrime     = zeros(nodeNum,1);

    for i = 1:nodeNum % ith polynomial
        chebyPolyPrime(i) = cos(real(i-1)*acos(capitalPrimeScaled));
    end

    % calculate residual  
    continuationValue = zeros(sizeShocks,1);
    for indexShocksPrime = 1:sizeShocks
        thetaVFunctionShock = theta(((indexShocksPrime-1)*nodeNum+1):indexShocksPrime*nodeNum,1);
        continuationValue(indexShocksPrime) = dot(thetaVFunctionShock,chebyPolyPrime);
    end
    
    consMarginalUtil = (log(consumption) - eeta * labor^2/2)^(rrho-1)/consumption;

    capPrimeMarginalProd   = zeros(sizeShocks,1);
    temp   = zeros(sizeShocks,1);
    mPrime = zeros(sizeShocks,1);

    for indexShocksPrime = 1:sizeShocks
        
        indexProdPrime   = ceil(indexShocksPrime/sizeAlpha);
        indexAlphaPrime  = mod(indexShocksPrime-1,sizeAlpha)+1;
        
        zPrime           = gridProd(indexProdPrime);
        aalphaPrime      = gridAlpha(indexAlphaPrime);

        thetaLaborShock  = theta(((indexShocksPrime-1)*nodeNum+1):indexShocksPrime*nodeNum,2);
        laborPrime       = dot(thetaLaborShock,chebyPolyPrime);

        outputPrime      = exp(zPrime)*capitalPrime^aalphaPrime*laborPrime^(1-aalphaPrime);
        consumptionPrime = (1-aalphaPrime)*exp(zPrime)*capitalPrime^aalphaPrime*laborPrime^(-aalphaPrime)/(eeta*laborPrime);
        
        consPrimeMarginalUtil                  = (log(consumptionPrime) - eeta * laborPrime^2/2)^(rrho-1)/consumptionPrime;
        capPrimeMarginalProd(indexShocksPrime) = aalphaPrime*exp(zPrime)*capitalPrime^(aalphaPrime-1)*laborPrime^(1-aalphaPrime);
        temp(indexShocksPrime)                 = continuationValue(indexShocksPrime)^(ppsi-rrho)*consPrimeMarginalUtil*...
            (capPrimeMarginalProd(indexShocksPrime)+(1-ddelta));
        mPrime(indexShocksPrime)               = bbeta*consPrimeMarginalUtil/consMarginalUtil*continuationValue(indexShocksPrime)^(ppsi-rrho);
    
    end    

    eulerRHS = bbeta*dot(transitionMatrix(indexShocks,:),continuationValue.^ppsi)^(rrho/ppsi-1)*dot(transitionMatrix(indexShocks,:),temp);

    eulerLHS   = (log(consumption) - eeta * labor^2/2)^(rrho-1)/consumption;

    eulerError = eulerLHS - eulerRHS;

    eulerError = log10(abs(eulerError));

    rb = (dot(transitionMatrix(indexShocks,:),continuationValue.^ppsi)^(rrho/ppsi-1)*dot(transitionMatrix(indexShocks,:),mPrime))^(-1)-1;    % bond return between t and t+1
    rk = aalpha*output/capital-ddelta;                   % risky asset return between t-1 and t
    rkCond = dot(transitionMatrix(indexShocks,:),capPrimeMarginalProd);        % conditional rk next period

    if( eulerError < -17 )
        eulerError = -17;
    end