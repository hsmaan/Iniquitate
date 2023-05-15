## Guidelines scripts - skeleton

Present a case where this is imbalance present
based on 1-2 cell-types in one batch but the 
other - probably Tran et al. experiment #7

1. Unsupervised clustering within each batch 
    - Check clusters and proportions 
        - Check if things differ
2. There isn't going to be a reference available for this - skip this step 
3. Start by tuning the tradeoff with some respect (e.g. using a method that better preserves heterogeneity - test 3 of them)
4. Measure degree of batch mixing and compare pre and post integration clusters
5. Assess tradeoff 
6. Compute ARI, Balanced ARI and ARIBatch 
7. Tune something wrt method or method parameters
8. Re-assess tradeoff - see if this is better 
9. Compute ARI, BAri and ARIBatch once again

- Ensure that proper written descriptions are included with everything 
- Need to add some of the OCSA diagnostics libraries in the environment 