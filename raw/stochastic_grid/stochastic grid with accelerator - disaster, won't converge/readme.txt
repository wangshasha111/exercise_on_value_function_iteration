Stochastic grid 
The code is a disaster for two reasons:

I use the accelerator scheme for stochastic grid, which won't make it converge

I randomly pick k' instead of both k' and k grid. But if we keep the k grid fixed, it will not improve efficiency.

So for a nice code, see the upper-level files