# Student-t Regression-with-Bayesian-inference
A Monte Carlo study comparing the Student-t and the classic regression model with bayesian inference.
This is a part of my thesis to get my degree as a statistician in Peru.

* [Application in House prices Kaggle dataset](https://github.com/diego-renato/Student-t-Regression-with-Bayesian-inference/blob/master/Application.ipynb)

So, what is the main difference of using student-t instead of normal errors in regression problems?
* In real problems, is commom to find outliers observations in the response variable or in the covariates and also is common to use normal errors, maybe because the normal distribution is very simple to interpret and the abundant theory in the last 100 years. But using normal errors the estimated parameters can be biased, so an alternative is to use a distribution with heavy tails than the normal one. In this project the classic regression (normal errors) and the student-t regression are compared by a Monte Carlo study to evaluate the performance of the estimated parameters. The results can be found clicking in [The explanation pdf](https://github.com/diego-renato/Student-t-Regression-with-Bayesian-inference/blob/master/The%20explanation.pdf).

