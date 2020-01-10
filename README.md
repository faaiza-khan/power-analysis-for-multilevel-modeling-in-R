# power-analysis-for-multilevel-modeling-in-R
It's difficult to estimate power in multilevel contexts. I planned to run an experiencing sampling study with repeated measures of emotional valence and social context, along with one-time measurements of attachment anxiety and attachment avoidance. I wrote the following R code to figure out (a) sample size and (b) number of experience sampling measurements required to ensure at least 80% power for all main and interaction effects of interest. 

The dependent variable is emotional valence. It is referred to as Valence in the code. It is measured on a scale from 1-7.
The Level 1 independent variable is social context (whether participants were alone or with others at the time of the measurement). It is referred to as Alone in the code and is a binary variable.
The Level 2 independent variables are attachment anxiety and attachment avoidance. They are referred to as Anxiety and Avoidance in the code. They are both measured on a scale from 1-7.
