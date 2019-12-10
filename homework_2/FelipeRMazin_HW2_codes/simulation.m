function [vSeriesDrop,kSeriesDrop,cSeriesDrop,lSeriesDrop,ySeriesDrop,...
    eeSeriesDrop,rbSeriesDrop,rkSeriesDrop,rkCondSeriesDrop] = ...
    simulation(theta,inputs)

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

T     = inputs.simulParams(4);
dropT = inputs.simulParams(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the economy for T periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    valueSeries       = zeros(T,1); 
    capitalSeries     = zeros(T,1);
    consumptionSeries = zeros(T,1);
    laborSeries       = zeros(T,1);
    outputSeries      = zeros(T,1);
    eulerErrorSeries  = zeros(T,1);
    rbSeries     = zeros(T,1);
    rkSeries     = zeros(T,1);
    rkCondSeries = zeros(T,1);

    % Generate Markov Chain
    rand('state',0);
    mu = zeros(1,sizeShocks);        % initial distribution
    mu(1,round(sizeShocks/2)) = 1;   % Always start from Z=0
    indexShocks = zeros(1,T+1);         % Stores indexes of the realization of "z" over the simulation
    indexShocks(1) = rando(mu);         % generate first ind value (time 0, not time 1)

    for i=1:T
        indexShocks(i+1) = rando(transitionMatrix(indexShocks(i),:));
    end
    
    capitalSeries(1) = inputs.steadyState;
    
    for indexTime = 1:T
        [valueSeries(indexTime),capitalSeries(indexTime+1),consumptionSeries(indexTime),laborSeries(indexTime),...
            outputSeries(indexTime),rbSeries(indexTime),...
            rkSeries(indexTime),rkCondSeries(indexTime),eulerErrorSeries(indexTime)] = ...
            eulerr_single(theta,inputs,indexShocks(indexTime),capitalSeries(indexTime)); 
    end

    vSeriesDrop      = valueSeries(dropT+1:T); 
    kSeriesDrop      = capitalSeries(dropT+1:T);    
    cSeriesDrop      = consumptionSeries(dropT+1:T);
    lSeriesDrop      = laborSeries(dropT+1:T);
    ySeriesDrop      = outputSeries(dropT+1:T);
    eeSeriesDrop     = eulerErrorSeries(dropT+1:T);
    rbSeriesDrop     = rbSeries(dropT+1:T);
    rkSeriesDrop     = rkSeries(dropT+1:T);
    rkCondSeriesDrop = rkCondSeries(dropT+1:T);