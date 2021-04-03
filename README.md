# VEwaning
Implements a statistical framework for evaluating the efficacy of vaccines  based on a potential outcomes formulation.


**Usage**
```
veWaning(data, L, ..., 
         lag = 0, 
         modelGam1 = NULL, modelGam2 = NULL, modelEntry = NULL, modelPsiGam1 = NULL, modelPsiGam2 = NULL, 
         gFunc = NULL, v = NULL, 
         minWgt = NULL, maxWgt = NULL, 
         txName = "A", infectionTime = "U", entryTime = "E", Gamma = "Gam", unblindTime = "R", vaccinated = "Psi")
```

**Arguments**
| Input Variable | Description
|----------------|--------------|
|data | A data.frame object containing all relevant data.
|L | A numeric object. The analysis time.
|... |Ignored. Used only to require named inputs.
|lag |A numeric object. The lag time between the initial vaccine dose and full efficacy.
|modelGam1 |A formula object. The coxph model for Gamma = 1. The LHS is set as the appropriate Surv() object internally. If a LHS is provided, it is ignored.
|modelGam2 |A formula object. The coxph model for Gamma = 2. The LHS is set as the appropriate Surv() object internally. If a LHS is provided, it is ignored.
|modelEntry |A formula object. The coxph model for entry times. The LHS is set as the appropriate Surv() object internally. If a LHS is provided, it is ignored.
|modelPsiGam1 |A formula object. The logistic model for vaccination for participants with Gamma = 1. If a LHS is provided, it is ignored.
|modelPsiGam2 |A formula object. The logistic model for vaccination for participants with Gamma = 2. If a LHS is provided, it is ignored.
|gFunc |A character object. The model of infection rates. Must be one of ’lin’, ’piece’, ’splin’, ’spcub’ for the linear, piecewise constant, linear spline, and cubic splinemodels respectively
|v |A numeric vector. The knots or cut-offs to be used by gFunc. If gFunc = ’lin’, this input is ignored. For ’splin’ and ’spcub’, this is the knots of the spline on (0,L). For ’piece’, v is the cut-offs on (0,L). Note that this input should not include the extremes 0 and L.
|minWgt |A numeric object. If not NULL, the minimum non-zero value a weight can have, i.e., weight = max(minWgt, weight). If NULL, no truncation of weights is performed.
|maxWgt |A numeric object. If not NULL, the maximum value a weight can have, i.e., weight = min(maxWgt, weight). If NULL, no truncation of weights is performed.
|txName |A character object. The header of the column of data containing the treatment variable. Default value is ’A’. Treatment must be coded as 0/1, where 1 indicates that participant was vaccinated; 0 otherwise. 
|infectionTime |A character object. The header of the column of data containing the time of infection on the scale of the calendar time. Default value is ’U’.
|entryTime |A character object. The header of the column of data containing the time of entry into the study on the scale of the calendar time. Default value is ’E’.
|Gamma |A character object. The header of the column of data containing the category for the unblinding dynamic. Default value is ’Gam’. Data must be 0/1/2, where 0 indicates infection occurs before requested/ offered unblinding; 1 indicates unblinding was requested by participant prior to the commencement of participant decision clinic visits; and 2 indicates that unblinding occurred after the commencement of participant decision clinic visits unblindTime A character object. The header of the column of data containing the time to requested unblinding, participant decision clinic visit/requested unblinding, or infection, whichever comes first. Default value is ’R’.
|vaccinated |A character object. The header of the column of data containing the indicator of vaccination, where 1 if participant is vaccinated; 0 otherwise. Default value is ’Psi’.

**Details**
Note the infection time, U, can take values NA or a value > L if the participant did not become
infected. All other data must be complete.

**Value**

A list object </br>
**theta** A vector object containing the estimated theta parameters.  </br>
**cov** The covariance estimated using the sandwich estimator. </br>
**SE** The standard error estimated using the sandwich estimator. 

**References**
Tsiatis, A. A. and Davidian, M. (2021) Estimating Vaccine Efficacy Over Time After a Randomized
Study is Unblinded. Submitted.
