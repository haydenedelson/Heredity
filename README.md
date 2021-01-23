# Heredity
Given a family tree and some information about population-level gene presence, creates a Bayesian Network for gene inheritance.

Input: family tree csv that indicates whether second-or-older generation family members exhibit a particlar trait (observation), from which genetic information (hidden state) may be inferred.  E.g.:

![Harry Potter family tree](https://github.com/haydenedelson/Heredity/blob/main/assets/Harry_family_tree.png)

Algorithm iterates over all possible combinations of gene values (i.e. 0, 1, or 2) for each family member. It computes the joint probability of each scenario and adds it to the appropriate probability distribution. Finally, it normalizes the distributions.

Output:

```
$ python heredity.py data.family0.csv
Harry:
  Gene:
    2: 0.0092
    1: 0.4557
    0: 0.5351
  Trait:
    True: 0.2665
    False: 0.7335
James:
  Gene:
    2: 0.1976
    1: 0.5106
    0: 0.2918
  Trait:
    True: 1.0000
    False: 0.0000
Lily:
  Gene:
    2: 0.0036
    1: 0.0136
    0: 0.9827
  Trait:
    True: 0.0000
    False: 1.0000
```
