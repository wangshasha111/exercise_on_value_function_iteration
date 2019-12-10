function f = sstate(x,laborSteadyState)

f(1) = 1 - 0.99*(0.3*x(1)^(0.3-1)*laborSteadyState^0.7 + 0.9);
f(2) = 0.7 * x(1)^0.3 * laborSteadyState^(-0.3) - (x(1)^0.3 * laborSteadyState^0.7 - 0.1 * x(1)) * x(2) * laborSteadyState;