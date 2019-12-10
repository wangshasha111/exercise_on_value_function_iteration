function [residual, g1, g2, g3] = perturbation_hw2_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 10, 1);

%
% Model equations
%

T30 = (1-params(1))*y(5)^params(2)+params(1)*y(7)^(params(2)/(1-params(3)));
T42 = y(3)^y(10);
T45 = y(2)^(1+y(10));
T46 = params(4)*T45;
T55 = params(1)*y(7)^((params(2)+params(3)-1)/(1-params(3)));
T56 = y(6)^(params(2)+params(3)-1);
T64 = y(3)^(y(10)-1);
T65 = y(10)*exp(y(9))*T64;
T66 = y(2)^(1-y(10));
lhs =y(5);
rhs =log(y(1))-params(4)*y(2)^2/2;
residual(1)= lhs-rhs;
lhs =y(6);
rhs =T30^(1/params(2));
residual(2)= lhs-rhs;
lhs =y(7);
rhs =y(6)^(1-params(3));
residual(3)= lhs-rhs;
lhs =y(1);
rhs =(1-y(10))*exp(y(9))*T42/T46;
residual(4)= lhs-rhs;
lhs =y(8);
rhs =T55/T56;
residual(5)= lhs-rhs;
lhs =1;
rhs =y(8)*y(4);
residual(6)= lhs-rhs;
lhs =y(4);
rhs =1+T65*T66-params(5);
residual(7)= lhs-rhs;
lhs =y(1)+y(3);
rhs =y(3)*(1-params(5))+T66*exp(y(9))*T42;
residual(8)= lhs-rhs;
lhs =y(9);
rhs =y(9)*params(6)+params(7)*x(1);
residual(9)= lhs-rhs;
lhs =y(10);
rhs =params(8)+y(10)*params(9)+params(10)*x(2);
residual(10)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(10, 10);

  %
  % Jacobian matrix
  %

T108 = getPowerDeriv(y(2),1-y(10),1);
T113 = getPowerDeriv(y(3),y(10),1);
T128 = getPowerDeriv(T30,1/params(2),1);
  g1(1,1)=(-(1/y(1)));
  g1(1,2)=params(4)*2*y(2)/2;
  g1(1,5)=1;
  g1(2,5)=(-((1-params(1))*getPowerDeriv(y(5),params(2),1)*T128));
  g1(2,6)=1;
  g1(2,7)=(-(T128*params(1)*getPowerDeriv(y(7),params(2)/(1-params(3)),1)));
  g1(3,6)=(-(getPowerDeriv(y(6),1-params(3),1)));
  g1(3,7)=1;
  g1(4,1)=1;
  g1(4,2)=(-((-((1-y(10))*exp(y(9))*T42*params(4)*getPowerDeriv(y(2),1+y(10),1)))/(T46*T46)));
  g1(4,3)=(-((1-y(10))*exp(y(9))*T113/T46));
  g1(4,9)=(-((1-y(10))*exp(y(9))*T42/T46));
  g1(4,10)=(-((T46*(T42*(-exp(y(9)))+(1-y(10))*exp(y(9))*T42*log(y(3)))-(1-y(10))*exp(y(9))*T42*params(4)*T45*log(y(2)))/(T46*T46)));
  g1(5,6)=(-((-(T55*getPowerDeriv(y(6),params(2)+params(3)-1,1)))/(T56*T56)));
  g1(5,7)=(-(params(1)*getPowerDeriv(y(7),(params(2)+params(3)-1)/(1-params(3)),1)/T56));
  g1(5,8)=1;
  g1(6,4)=(-y(8));
  g1(6,8)=(-y(4));
  g1(7,2)=(-(T65*T108));
  g1(7,3)=(-(T66*y(10)*exp(y(9))*getPowerDeriv(y(3),y(10)-1,1)));
  g1(7,4)=1;
  g1(7,9)=(-(T65*T66));
  g1(7,10)=(-(T66*(exp(y(9))*T64+y(10)*exp(y(9))*T64*log(y(3)))+T65*T66*(-log(y(2)))));
  g1(8,1)=1;
  g1(8,2)=(-(exp(y(9))*T42*T108));
  g1(8,3)=1-(1-params(5)+T66*exp(y(9))*T113);
  g1(8,9)=(-(T66*exp(y(9))*T42));
  g1(8,10)=(-(exp(y(9))*T42*T66*(-log(y(2)))+T66*exp(y(9))*T42*log(y(3))));
  g1(9,9)=1-params(6);
  g1(10,10)=1-params(9);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,100);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],10,1000);
end
end
end
end
