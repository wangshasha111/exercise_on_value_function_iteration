// ECON 714
// Homework 2
// Question 8 - perturbation

//----------------------------------------------------------------
// 0. Housekeeping
//----------------------------------------------------------------

close all

//----------------------------------------------------------------
// 1. Endogenous variables (10=4+4+2)
//----------------------------------------------------------------

var 

// Allocation (4)
c l k r

// Utility (4)
u v ev m

// Shocks (2)
z alpha;

//----------------------------------------------------------------
// 2. Exogenous variables (2)
//----------------------------------------------------------------
 
varexo 

ez ealpha;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

// Preferences
beta theta gamma eta

// Technology
delta

// Stochastic processes
rhoz sigmaz mualpha rhoalpha sigmaalpha;

//----------------------------------------------------------------
// 4. Calibration
//----------------------------------------------------------------

// Preferences
beta       = 0.99;
theta      = 0.5;
gamma      = 10;

// Technology
delta      = 0.1;  

// Stochastic processes
rhoz       = 0.95;
sigmaz     = 0.005;
mualpha    = 0.03;
rhoalpha   = 0.9;
sigmaalpha = 0.01;

//----------------------------------------------------------------
// 4. Steady State
//----------------------------------------------------------------


l_ss     = 100;
alpha_ss = 0.3;
k_ss     = (beta*alpha_ss*l_ss^(1-alpha_ss)/(1+beta*(delta-1)))^(1/(1-alpha_ss));
c_ss     = -delta*k_ss + k_ss^alpha_ss*l_ss^(1-alpha_ss);
eta      = (1-alpha_ss)*k_ss^alpha_ss*l_ss^(-1-alpha_ss)/c_ss;
u_ss     = (log(c_ss)-eta*l_ss^2/2);
v_ss     = (log(c_ss)-eta*l_ss^2/2);
ev_ss    = v_ss^(1-gamma);
m_ss     = beta;
r_ss     = (alpha_ss*k_ss^(alpha_ss-1)*l_ss^(1-alpha_ss)+1-delta); 

//----------------------------------------------------------------
// 5. Model
//----------------------------------------------------------------

model;

  // 1. Utility
  u = log(c)-eta*l^2/2;

  // 2. Value function
  v = ((1-beta)*u^theta + beta*ev^(theta/(1-gamma)))^(1/theta);

  // 3. Expected value function
  ev = v(+1)^(1-gamma);

  // 4. Static leisure-consumption
  c = (1-alpha)*exp(z)*k(-1)^alpha/(eta*l^(1+alpha));

  // 5. Pricing Kernel
  m = beta*(c(-1)/c)*(u/u(-1))^(theta-1)*ev(-1)^((theta+gamma-1)/(1-gamma))/v^(theta+gamma-1);
	     
  // 6. Euler equation for  capital
  1 = m(+1)*r(+1);
  
  // 7. Return on capital net of depreciation
  r = (exp(z)*alpha*k(-1)^(alpha-1)*l^(1-alpha)+1-delta);

  // 8. Resource constraint 
  k+c = (1-delta)*k(-1)+exp(z)*k(-1)^alpha*l^(1-alpha);

  // 9. Productivity process 
  z = rhoz*z(-1)+sigmaz*ez;

  // 10. Capital share process
  alpha = mualpha + rhoalpha*alpha(-1)+sigmaalpha*ealpha;
  
end;

//----------------------------------------------------------------
// 6. Computation
//----------------------------------------------------------------

initval;
  c      = c_ss;
  l      = l_ss;
  k      = k_ss;
  r      = r_ss;
  v      = v_ss;
  u      = u_ss;
  ev     = ev_ss;
  m      = m_ss;
  alpha  = alpha_ss;
  ealpha = 0;
  z      = 0;
  ez     = 0;
end;

shocks;
	     var ez = 1;
	     var ealpha = 1;
end;

steady;

check;

stoch_simul(irf = 100, order = 3) k l c z alpha;
