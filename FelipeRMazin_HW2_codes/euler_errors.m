function e = euler_errors(valueFunction,policyFunction,laborFunction,consumptionFunction)

global inputs

bbeta  = inputs.params(1);
ddelta = inputs.params(2);
ppsi   = inputs.params(3);
rrho   = inputs.params(4);
eeta   = inputs.params(5);

gridCapital = inputs.gridCapital;
sizeCapital = length(gridCapital);
gridProd    = inputs.gridProd;
gridAlpha   = inputs.gridAlpha;
sizeProd    = inputs.simulParams(2);
sizeAlpha   = inputs.simulParams(3);
sizeShocks  = sizeProd*sizeAlpha;

transitionMatrix = inputs.transitionMatrix;

expectedValue = (valueFunction.^ppsi) * transitionMatrix';
for indexShocks = 1:sizeShocks
    indexProd  = ceil(indexShocks/sizeAlpha);
    indexAlpha = mod(indexShocks-1,sizeAlpha)+1;
    
    aalpha     = gridAlpha(indexAlpha);
    z          = gridProd(indexProd);
    for indexCapital = 1:sizeCapital

        capital      = gridCapital(indexCapital);
        capitalPrime = policyFunction(indexCapital,indexShocks);
        labor        = laborFunction(indexCapital,indexShocks);
        consumption  = consumptionFunction(indexCapital,indexShocks);

        continuationValue = zeros(sizeShocks,1);
        temp              = zeros(sizeShocks,1);

        for indexShocksPrime = 1:sizeShocks

            indexProdPrime  = ceil(indexShocksPrime/sizeAlpha);
            indexAlphaPrime = mod(indexShocksPrime-1,sizeAlpha)+1;
            aalphaPrime     = gridAlpha(indexAlphaPrime);
            zPrime          = gridProd(indexProdPrime);
            
            indexCapitalLow  = max(sum(capitalPrime > gridCapital),1);
            indexCapitalHigh = indexCapitalLow + 1;

            continuationValue(indexShocksPrime) = valueFunction(indexCapitalLow,indexShocks) + (valueFunction(indexCapitalHigh,indexShocks) - ...
                valueFunction(indexCapitalLow,indexShocks))./(gridCapital(indexCapitalHigh) - gridCapital(indexCapitalLow)).*(capitalPrime - gridCapital(indexCapitalLow));
            laborPrime        = laborFunction(indexCapitalLow,indexShocks) + (laborFunction(indexCapitalHigh,indexShocks) - ...
                laborFunction(indexCapitalLow,indexShocks))./(gridCapital(indexCapitalHigh) - gridCapital(indexCapitalLow)).*(capitalPrime - gridCapital(indexCapitalLow));

            consumptionPrime   = ((1 - aalphaPrime) * exp(zPrime)*capitalPrime^aalphaPrime*laborPrime^(-aalphaPrime))/(eeta*laborPrime);

            consPrimeMarginalUtil  = (log(consumptionPrime) - eeta * laborPrime^2/2)^(rrho-1)/consumptionPrime;
            capPrimeMarginalProd   = aalphaPrime*exp(zPrime)*capitalPrime^(aalphaPrime-1)*laborPrime^(1-aalphaPrime);
            temp(indexShocksPrime) = continuationValue(indexShocksPrime)^(ppsi-rrho)*consPrimeMarginalUtil*(capPrimeMarginalProd+(1-ddelta));

        end

        eulerRHS = bbeta*dot(transitionMatrix(indexShocks,:),continuationValue.^ppsi)^(rrho/ppsi-1)*dot(transitionMatrix(indexShocks,:),temp);


        eulerLHS   = (log(consumption) - eeta * labor^2/2)^(rrho-1)/consumption;

        e(indexCapital,indexShocks) = eulerLHS - eulerRHS;
    
    
    end

end

e = log10(abs(e));