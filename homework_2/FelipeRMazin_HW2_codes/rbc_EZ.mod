// Basic RBC Model with EZ preferences
// 
// Jesus Fernandez-Villaverde
//
// Bala Cynwyd, August 16, 2010

//----------------------------------------------------------------
// 0. Housekeeping
//----------------------------------------------------------------

close all

//----------------------------------------------------------------
// 1. Endogenous variables
//----------------------------------------------------------------

var 

// Utility variables
utility valueFunction expectedValue kernel

// Allocation variables 
output consumption labor capital investment

// Input prices
interestRate wage

// Productivity
z

// Capital Share
aalpha;

//----------------------------------------------------------------
// 2. Exogenous variables
//----------------------------------------------------------------

varexo epsZ epsAlpha;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

// Utility function
bbeta eeta ppsi rrho

// Technology 
ddelta ggamma pphi0 pphi1 sigmaEpsilonZ sigmaEpsilonAlpha;

//----------------------------------------------------------------
// 4. Calibration
//----------------------------------------------------------------

// Preferences 
bbeta   = 0.99;     // Discount factor
rrho    = 0.5;      // Elasticity of intertemporal substitution factor: 1/(1-rrho)
ppsi    = -9;       // Risk aversion factor

// Technology
ddelta            = 0.1;  // Depreciation rate
ggamma            = 0.95;
pphi0             = 0.03;
pphi1             = 0.9;
sigmaEpsilonZ     = 0.005;
sigmaEpsilonAlpha = 0.01;

//----------------------------------------------------------------
// 5. Steady State
//----------------------------------------------------------------


laborSS         = 50;
aalphaSS        = 0.3;
zSS             = 0;
capitalSS       = (bbeta*aalphaSS*laborSS^(1-aalphaSS)/(1+bbeta*(ddelta-1)))^(1/(1-aalphaSS));
outputSS        = capitalSS^aalphaSS*laborSS^(1-aalphaSS);
consumptionSS   = outputSS-ddelta*capitalSS;
eeta            = (1-aalphaSS)*capitalSS^aalphaSS*laborSS^(-1-aalphaSS)/consumptionSS;
utilitySS       = (log(consumptionSS)-eeta*laborSS^2/2);
valueFunctionSS = (log(consumptionSS)-eeta*laborSS^2/2);
expectedValueSS = valueFunctionSS^ppsi;
interestRateSS  = aalphaSS*capitalSS^(aalphaSS-1)*laborSS^(1-aalphaSS);
kernelSS        = bbeta;

//----------------------------------------------------------------
// 6. Model
//----------------------------------------------------------------

model; 
  
  // 1. Utility
  utility = (log(consumption) - eeta * labor^2/2);
  
  // 2. Value function
  valueFunction = ((1-bbeta)*utility^rrho +bbeta*expectedValue^(rrho/ppsi))^(1/rrho);
  
  // 3. Expected value function
  expectedValue = valueFunction(+1)^ppsi;

  // 4. Static leisure-consumption
  wage = consumption * eeta * labor;
  
  // 5. Pricing Kernel
  kernel = bbeta*expectedValue(-1)^(rrho/ppsi-1)*(valueFunction^(ppsi-rrho)*(utility/utility(-1))^(rrho-1)*
  consumption(-1)/consumption);

  // 6. Euler equation
  kernel(+1)* (interestRate(+1) + (1-ddelta)) = 1; 
  
  // 7. Production function
  output = exp(z)*(capital(-1)^aalpha)*labor^(1-aalpha);
  
  // 8. Interest rate
  interestRate = aalpha*output/(capital(-1));

  // 9. Wage
  wage = (1-aalpha)*output/labor;
  
  // 10. Resource constraint
  consumption+investment = output; 

  // 11. Law of motion for capital
  capital = (1-ddelta)*capital(-1)+investment;  

  // 12. Law of motion for productivity
  z = ggamma * z(-1) + sigmaEpsilonZ * epsZ;
  
  // 13. Law of motion for capital share
  aalpha = pphi0 + pphi1 * aalpha(-1) + sigmaEpsilonAlpha * epsAlpha;

end;

//----------------------------------------------------------------
// 7. Computation
//----------------------------------------------------------------

initval;
  labor         = laborSS;
  capital       = capitalSS;
  output        = outputSS;
  consumption   = consumptionSS;
  investment    = ddelta*capital;
  interestRate  = interestRateSS;
  wage          = (1-aalpha)*outputSS/laborSS;
  kernel        = kernelSS;
  utility       = log(consumptionSS)-eeta*laborSS^2/2;
  valueFunction = utilitySS;
  expectedValue = utilitySS^ppsi;
  z             = zSS;
  aalpha        = aalphaSS;
  epsZ          = 0;
  epsAlpha      = 0;
end;
    
shocks;
  var epsZ = 1;
  var epsAlpha = 1;
end;

steady (solve_algo = 3);  

check;

stoch_simul(irf = 100, order = 3) capital labor consumption z aalpha;