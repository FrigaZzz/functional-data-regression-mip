## Proposed Formulation for feature selection without Sparsity

**L0 Norm with Group Constraints**: 

The objective function consists of the sum of squared errors for the regression model, a penalty term proportional to the L0 norm of the coefficients (modeled by the binary variables), and a group sparsity penalty term. The constraints ensure that coefficients are either fixed at zero or allowed to vary freely when their corresponding binary indicators are active, and they maintain the logical connections between individual coefficients and their respective groups.
No coefficient sparsity is considered.

\[
\min_{\gamma, \text{group}} 
\frac{1}{n}\sum_{i=1}^{n}\left(Y_i - \sum_{j=1}^{p}\sum_{k=1}^{r} Z_{ijk} \beta_{jk}\right)^2 + 
\lambda \sum_{j=1}^{p}\sum_{k=1}^{r} | \beta_{jk} |
\]

Subject to the following constraints:
\[
-\text{BIG\_M} \cdot \text{Group}_{j} \leq \beta_{jk} \leq \text{BIG\_M} \cdot \text{Group}_{j} \quad \forall j=1,...,p; \ k=1,...,r
\]
\[
\sum_{j=1}^{p} \text{group}_j \leq \text{group\_limit}
\]

\[
\beta_{jk} \in \mathbb{R},  \quad \forall j=1,...,p; \ k=1,...,r
\]

\[
\text{group}_j \in \{0,1\} \quad \forall j=1,...,p
\]
**Coefficient-Binary Link Constraints:**
\[
-\text{BIG\_M} \cdot \text{Group}_{j} \leq \beta_{jk} \leq \text{BIG\_M} \cdot \text{Group}_{j} \quad \forall j=1,...,p; \ k=1,...,r
\]
These constraints ensure that the regression coefficients $\beta_{jk}$ can only be non-zero when their corresponding group is set to 1.
It is required to correctly link the Group variable to Betas. 

**Group Limit Constraint:**
\[
\sum_{j=1}^{p} \text{group}_j \leq \text{group\_limit}
\]
This constraint enforces that the total number of groups with non-zero coefficients does not exceed a predefined limit $\text{group\_limit}$, thus controlling the group-level sparsity.



And the decision variables are bounded as:

\[
\beta_{jk} \in \mathbb{R},  \quad \forall j=1,...,p; \ k=1,...,r
\]

\[
\text{group}_j \in \{0,1\} \quad \forall j=1,...,p
\]

These constraints define the domains of the decision variables, stating that regression coefficients are real numbers, whereas the binary variables indicating non-zero coefficients and group selections are binary.


In this formulation:

- $\beta_{jk}$ are the continuous decision variables representing the coefficients of the regression model.
- $\text{group}_j$ are binary variables indicating whether any coefficients in the $j$-th group are non-zero.
- $\text{BIG\_M}$ is a sufficiently large positive constant used in the "big M" method to effectively turn constraints on or off based on the value of the binary variables.
- $\text{group\_limit}$ represents the maximum number of groups that can have non-zero coefficients.

