This package provides the testing and algorithm files to solve the following **Best Subset Selection (BSS)** problem
```math
\min_{ \beta \in R^p} \|y-X \beta \|_{2}^{2}  \text{ subject to } \|\beta\|_{0}\leq k \text{ , } l\leq \beta \leq u
``` 
where $l$, $u$, $\beta \in R^p$, $y \in R^n$, $X \in R^{n \times p}$, $`\| \cdot \|_{2}`$ is $`l_{2}`$ norm, $`\|\cdot\|_{0}`$ is a pseudo-norm which counts the no. of non-zero entries of $`\beta`$.\
The following five suboptimal algorithms have been implemented
and compared to solve the BSS problem.
1. Discrete First Order
2. Forward Selection
3. Genetic Algorithm
4. Sequential Feature Swapping
5. Sequential Floating Feature Selection.
