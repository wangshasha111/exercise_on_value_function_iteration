function theta = residual_function(thetaGuess,inputs)

bbeta  = inputs.params(1);
ddelta = inputs.params(2);
ppsi   = inputs.params(3);
rrho   = inputs.params(4);
eeta   = inputs.params(5);

gridCapital = inputs.gridCapital;
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

% Solves for the coefficients associated to the Chebychev polynomials. 

options = optimset('Display','Iter','TolFun',10^(-15),'TolX',10^(-15));
theta = fsolve(@notime_iter,thetaGuess,options);

%% Residual Function
    function residual = notime_iter(theta)
        
        residualSection = zeros(nodeNum,2);
        residual = zeros(M,2);

        for indexShocks = 1:sizeShocks
            
            indexProd  = ceil(indexShocks/sizeAlpha);
            indexAlpha = mod(indexShocks-1,sizeAlpha)+1;
            aalpha     = gridAlpha(indexAlpha);
            z          = gridProd(indexProd);
    
            value = zeros(nodeNum,1);
            policyLabor       = zeros(nodeNum,1);
            policyCapital     = zeros(nodeNum,1);
            policyConsumption = zeros(nodeNum,1);
    
            thetaVFunctionShock = theta(((indexShocks-1)*nodeNum+1):indexShocks*nodeNum,1);
            thetaLaborShock     = theta(((indexShocks-1)*nodeNum+1):indexShocks*nodeNum,2);
    
            for indexCapital = 1:nodeNum % Loop 1 over collocation points on capital
        
                value(indexCapital) = dot(thetaVFunctionShock,chebyPoly(indexCapital,:));  % value function at each collocation points
                labor = dot(thetaLaborShock,chebyPoly(indexCapital,:));               % labor at each collocation points
                capital = gridCapital(indexCapital);

                if labor < 2
                    labor = 2;
                    disp('labor breaks lower bound')
                elseif labor > 120
                    labor = 120;
                    disp('labor breaks upper bound')
                end
                
                policyLabor(indexCapital) = labor;

                output      = exp(z)*capital^aalpha*labor^(1-aalpha);
                consumption = (1-aalpha)*exp(z)*capital^aalpha*labor^(-aalpha)/(eeta*labor);

                capitalPrime = output+(1-ddelta)*capital-consumption;
                if capitalPrime < capitalMin
                    capitalPrime = capitalMin+0.01;
                    %disp('capital prime breaks lower bound')
                elseif capitalPrime > capitalMax
                    capitalPrime = min(output+(1-ddelta)*capital-consumption,capitalMax) - 0.01;
                    %disp('capital prime breaks upper bound')
                end
            
                policyCapital(indexCapital) = capitalPrime;

                consumption = output+(1-ddelta)*capital-capitalPrime;
            
                if consumption < 0
                    disp('warning: consumption < 0')
                    consumption = 10^(-10);
                end
            
                policyConsumption(indexCapital) = consumption;

            end % Loop 1 over collocation point on k ends

            % scale capital prime from [capitalMin,capitalMax] to [-1,1]
            policyCapitalScaledDown = (2*policyCapital-(capitalMin+capitalMax))/(capitalMax-capitalMin);
    
            % value of polynomials at each scaled k prime
            chebyPolyPrime      = ones(nodeNum,nodeNum);
            chebyPolyPrime(:,2) = policyCapitalScaledDown;
            
            for i1=3:nodeNum
                chebyPolyPrime(:,i1) = 2*policyCapitalScaledDown.*chebyPolyPrime(:,i1-1)-chebyPolyPrime(:,i1-2);
            end
    
            % Calculate residual
            for indexCapital = 1:nodeNum % Loop 2 over collocation points on capital
        
                continuationValue = zeros(sizeShocks,1);
                temp              = zeros(sizeShocks,1);
                capitalPrime      = policyCapital(indexCapital);
            
                for indexShocksPrime = 1:sizeShocks
                    
                    indexProdPrime  = ceil(indexShocksPrime/sizeAlpha);
                    indexAlphaPrime = mod(indexShocksPrime-1,sizeAlpha)+1;
                    aalphaPrime     = gridAlpha(indexAlphaPrime);
                    zPrime          = gridProd(indexProdPrime);

                    thetaVFunctionShockPrime             = theta(((indexShocksPrime-1)*nodeNum+1):indexShocksPrime*nodeNum,1);
                    thetaLaborShockPrime                 = theta(((indexShocksPrime-1)*nodeNum+1):indexShocksPrime*nodeNum,2);
                    continuationValue(indexShocksPrime)  = dot(thetaVFunctionShockPrime,chebyPolyPrime(indexCapital,:));
                    laborPrime                           = dot(thetaLaborShockPrime,chebyPolyPrime(indexCapital,:));     

                    if laborPrime < 2
                        laborPrime = 2;
                        disp('labor tomorrow breaks lower bound')
                    elseif laborPrime > 120
                        laborPrime = 120;
                        disp('labor tomorrow breaks upper bound')
                    end

                    outputPrime        = exp(zPrime)*capitalPrime^aalphaPrime*laborPrime^(1-aalphaPrime);
                    consumptionPrime   = ((1 - aalphaPrime) * exp(zPrime)*capitalPrime^aalphaPrime*laborPrime^(-aalphaPrime))/(eeta*laborPrime);
                    capitalDoublePrime = outputPrime+(1-ddelta)*capitalPrime-consumptionPrime;

                    if capitalDoublePrime < capitalMin
                        capitalDoublePrime = capitalMin+0.01;
                        %disp('capital the day after tomorrow breaks lower bound')
                    elseif capitalDoublePrime > capitalMax
                        capitalDoublePrime = min(outputPrime+(1-ddelta)*capitalPrime-consumptionPrime,capitalMax) - 0.01;
                        %disp('capital the day after tomorrow breaks upper bound')
                    end

                    consumptionPrime = outputPrime+(1-ddelta)*capitalPrime-capitalDoublePrime;
                    if consumptionPrime < 0
                        disp('warning: consumptionPrime < 0')
                        consumptionPrime = 10^(-10);
                    end             

                    consPrimeMarginalUtil  = (log(consumptionPrime) - eeta * laborPrime^2/2)^(rrho-1)/consumptionPrime;
                    capPrimeMarginalProd   = aalphaPrime*exp(zPrime)*capitalPrime^(aalphaPrime-1)*laborPrime^(1-aalphaPrime);
                    temp(indexShocksPrime) = continuationValue(indexShocksPrime)^(ppsi-rrho)*consPrimeMarginalUtil*(capPrimeMarginalProd+(1-ddelta));
            
                end

                eulerRHS   = bbeta*dot(transitionMatrix(indexShocks,:),continuationValue.^ppsi)^(rrho/ppsi-1)*dot(transitionMatrix(indexShocks,:),temp);

                labor       = policyLabor(indexCapital);
                consumption = policyConsumption(indexCapital);
                
                eulerLHS   = (log(consumption) - eeta * labor^2/2)^(rrho-1)/consumption;
                bellmanRHS = ((1-bbeta)*(log(consumption) - eeta * labor^2/2)^rrho...
                              +bbeta*dot(transitionMatrix(indexShocks,:),continuationValue.^ppsi)^(rrho/ppsi))^(1/rrho);
                bellmanLHS = value(indexCapital);

                residualSection(indexCapital,1) = eulerRHS - eulerLHS;
                residualSection(indexCapital,2) = bellmanRHS - bellmanLHS;

            end % Loop 2 over k ends
    
            residual(((indexShocks-1)*nodeNum+1):indexShocks*nodeNum,:) = residualSection;

         end

      end

  end