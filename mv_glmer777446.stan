data {
    int n;
    int n_RE;
    int N1;
    int ncx1;
    int id1[N1];
    int RE_ind1[3];
    int colmns_nHC1[4];
    vector[N1] y1;
    matrix[N1, 3] Z1;
    matrix[n, ncx1] Xhc1;
    matrix[N1, ncx1] X1;
    real<lower=0> scale_betas1;
    real<lower=0> scale_sigmas;
    real<lower=0> scale_diag_D;
    real<lower=0> lkj_shape;
}
 
parameters {
    vector[ncx1] betas1;
    real<lower = 0> sigma1;
    matrix[n, n_RE] u;
    vector<lower = 0>[n_RE] L_var_D;
    cholesky_factor_corr[n_RE] L_corr_D;
}
 
transformed parameters {
    vector[N1] eta1;
    matrix[n, n_RE] mu_u;
    for (i in 1:n) {
        mu_u[i, 1] = dot_product(Xhc1[i, 1:2:3], betas1[1:2:3]);
        mu_u[i, 2] = 0.0;
        mu_u[i, 3] = 0.0;
    }
    for (j1 in 1:N1) {
        eta1[j1] = dot_product(Z1[j1, ], u[id1[j1], RE_ind1])
                 + dot_product(X1[j1, colmns_nHC1], betas1[colmns_nHC1]);
    }
}
 
model {
    matrix[n_RE, n_RE] L_D;
    L_D = diag_pre_multiply(L_var_D, L_corr_D);
    L_var_D ~ cauchy(0, scale_diag_D);
    L_corr_D ~ lkj_corr_cholesky(lkj_shape);
    for (i in 1:n) {
        u[i, ] ~ multi_normal_cholesky(mu_u[i, ], L_D);
    }
    for (k1 in 1:ncx1) {
        betas1[k1] ~ normal(0.0, scale_betas1);
    }
    sigma1 ~ cauchy(0, scale_sigmas);
    y1 ~ normal(eta1, sigma1);
}
 
generated quantities {
    matrix[n_RE, n_RE] D;
    matrix[n, n_RE] b;
    D = diag_pre_multiply(L_var_D, L_corr_D) * diag_pre_multiply(L_var_D, L_corr_D)';
    b = u - mu_u;
}
