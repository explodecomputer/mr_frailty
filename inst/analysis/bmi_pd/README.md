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


