function [f, labor, consumption] = value_function(capitalPrime,indexCapital,indexShocks)

global inputs

bbeta  = inputs.params(1);
ddelta = inputs.params(2);
ppsi   = inputs.params(3);
rrho   = inputs.params(4);
eeta   = inputs.params(5);

gridCapital = inputs.gridCapital;
gridProd    = inputs.gridProd;
gridAlpha   = inputs.gridAlpha;
sizeAlpha   = inputs.simulParams(3);

indexProd  = ceil(indexShocks/sizeAlpha);
indexAlpha = mod(indexShocks-1,sizeAlpha)+1;

capital    = gridCapital(indexCapital);
aalpha     = gridAlpha(indexAlpha);
z          = gridProd(indexProd);

expectedValue0 = inputs.expectedValue0;
laborFunction  = inputs.laborFunction;

options = optimset('Display','off');
l = @(labor) (1-aalpha) * exp(z) * capital^aalpha * labor^(-1-aalpha) - (exp(z)*capital^aalpha * ...
    labor^(1-aalpha) + (1 - ddelta) * capital - capitalPrime) * eeta;

labor = fzero(@(labor) l(labor),laborFunction(indexCapital,indexShocks),options);

consumption = exp(z)*capital^aalpha*labor^(1-aalpha) + (1 - ddelta)*capital - capitalPrime;

indexCapitalLow  = max(sum(capitalPrime > gridCapital),1);
indexCapitalHigh = indexCapitalLow + 1;

currentUtility = log(consumption) - eeta * labor^2/2;

if (consumption >= 0) && (labor >= 0) && (currentUtility >= 0)

    expectedValue  = expectedValue0(indexCapitalLow,indexShocks) + (expectedValue0(indexCapitalHigh,indexShocks) - ...
        expectedValue0(indexCapitalLow,indexShocks))./(gridCapital(indexCapitalHigh) - gridCapital(indexCapitalLow)).*(capitalPrime - gridCapital(indexCapitalLow));

    f = ((1-bbeta)*currentUtility^rrho + bbeta * expectedValue^(rrho/ppsi))^(1/rrho);
    
else
    f = -1e10;
end
