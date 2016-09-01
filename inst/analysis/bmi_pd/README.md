# BMI and PD frailty associations

```
exposure ~ snp(s)
mortality ~ age + exposure
outcome ~ age
```

1. simulate population with SNP(s) from Locke et al 2015.
2. simulate influence of BMI on mortality
3. simulate age related incidence of PD that is unrelated to BMI
4. perform MR using surviving population


## Model 1

BMI mortality curves are from [http://www.nejm.org/doi/full/10.1056/NEJMoa1000367](http://www.nejm.org/doi/full/10.1056/NEJMoa1000367)

This has a J-shaped association - both high and low BMI increase mortality


## Model 2

BMI has 5 times increase on mortality if BMI > 27. This is a control model.


## Model 3

Model 1 performed again using the Locke SNPs


## Model 4

Using linear hazard ratio of 1.16 obtained from [Davey Smith et al 2009 BMJ](http://www.bmj.com/content/339/bmj.b5043). This uses the causal effect of BMI on mortality rather than the J-shaped observational effect.