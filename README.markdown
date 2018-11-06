**Example of using monte carlo methods to produce a a function that returns a p value given a value**

This script works by taking simulating a co-occurrence matrix and performing a significance metric on each cell in this matrix. After you have your matrix of significance factors you can then runthe cdf function to find what level of significance metric has a significant p value.

This script provides an example of how to construct p values from statistics that don't necessarily have obvious mappings to significance
or whose distributions may be heteroscedastic.

There are plenty of statistics methods within the base R package that allow for easy and powerful manipulation of mock null distributions
that allow you to find p values for whatever significance statistic that you want to work with.
