# exercise_on_value_function_iteration

This file contains the [project description](https://github.com/wangshasha111/exercise_on_value_function_iteration/blob/master/homework1.pdf), [programming](https://github.com/wangshasha111/exercise_on_value_function_iteration/tree/master/homework1_code_send), and [writeup](https://github.com/wangshasha111/exercise_on_value_function_iteration/blob/master/ShashaWang_homework1.pdf) of the assignment of computational economics on value function iteration.

## Goals
The model is a household life-cycle model with home production and capital accumulation. 
The project aims to test the performance of different methods of value function iteration (VFI hereafter), namely, 

1. Value Function Iteration with a Fixed Grid, i.e., fixing a grid of 250 points of capital, centered around kss with a coverage of 30% of kss and equally spaced.
2. Value Function Iteration with an Endogenous Grid
3. Accelerator, i.e., skipping the max operator in the Bellman equation nine out of each ten times
4. Multigrid, i.e., scheme (Chow and Tsitsiklis, 1991) for a Value function iteration, with the grid centered around k_ss with a coverage of ±30% of k_ss and equally spaced
5. Stochastic Grid (Rust, 1997) for a Value function iteration, with 500 vertex points with a coverage of 25% of k_ss (you can keep the grid of investment fixed). 

## Takeaways
Generally speaking, accelerator, multigrid, and stochastic grid are all effective ways to speed up computation without compromising accuracy. 

