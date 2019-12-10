function [value,labor,consumption] = getValue(productivity,alpha,capitalToday,...
    capitalTomorrow,parameters,functions,laborMin,laborMax,interk,interv)
                

theta = parameters.theta;
gamma = parameters.gamma;
beta = parameters.beta;

getLabor = functions.getLabor;
getConsumption = functions.getConsumption;
getPeriodUtility = functions.getPeriodUtility;

labor = getLabor(productivity,capitalToday,capitalTomorrow,alpha,laborMin,laborMax);
consumption = getConsumption(productivity,capitalToday,labor,alpha);
periodUtility = getPeriodUtility(consumption,labor);

 if periodUtility < 0
     periodUtility = 0;
 end

value = ((1-beta)*periodUtility^theta+...
                    beta*interp1(interk,interv,capitalTomorrow)^(theta/(1-gamma)))^(1/theta);