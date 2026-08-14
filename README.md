This package provides the testing and algorithm files to solve the following **Best Subset Selection (BSS)** problem
```math
\min_{ \beta \in R^p} \|y-X \beta \|_{2}^{2}  \text{ subject to } \|\beta\|_{0}\leq k
``` 
where $\beta \in R^p$, $y \in R^n$, $X \in R^{n \times p}$, $`\| \cdot \|_{2}`$ is $`l_{2}`$ norm, $`\|\cdot\|_{0}`$ is a pseudo-norm which counts the no. of non-zero entries of $`\beta`$.\
The following five suboptimal algorithms have been implemented
and compared to solve the BSS problem.
1. Sequential Feature Swapping (SF) 
2. Forward Selection (FS)
3. Sequential Floating Feature Selection (SFFS)
5. Genetic Algorithm (GA)
6. Discrete First Order (DFO)
7. Adaptive Best Subset Selection (ABESS)
8. Mixed Integer Optimization (MIO)
