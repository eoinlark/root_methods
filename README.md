# R_root_methods
Collection of root finding methods implemented in C and callable through R

## Overview
C++ code can be implemented in R to create functions that are optimised in terms of speed relative to default R functions.

Code provides implementations for the following functions in C++ that are callable through R via the Rcpp package.
* Bisection Algorithim
* Netwon-Raphson
* Secant-Method

### Bisection Algorithim
The algorithm takes two input values a and b. It calculates an estimated root m by taking the average of both a and b. The algorithm then calculates the first derivative of the function with respect to a and m and determines if the sign of the function changes over the interval. If the sign changes, this indicates that a local root is present and the function continues to iterate. If the sign doesn’t change, it means that the algorithm has overshot and the value of a is switched with than of b, before continuing to iterate.

### Newton-Raphson
The algorithm calculates the first and second derivative of the function for the input value x1. It then approximates the function and finds the point x2 where the function crosses zero. The algorithm contimues to iterate until the first derivate is zero or the change in the value of the points x1 and x2 is less than zero.

### Secant-Method
This approach is identical to the Newton-Rapson approach, however the second derivative is not explicitly required; rather this is approximated with a finite difference approximation.

## Performance
Benchmarking each of the functions in R, Newton's method proves to be the most efficient, at the expense of requiring an initial accurate "guess" for the root. Otherwise there is a risk that the function will fail to converge.

## Links
[Bisection Algorithim](https://en.wikipedia.org/wiki/Bisection_method)
[Newton-Raphson](https://en.wikipedia.org/wiki/Newton%27s_method)
[Secant Methold](https://en.wikipedia.org/wiki/Secant_method)
