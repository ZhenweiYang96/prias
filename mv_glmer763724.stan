data {
    int n;
    int n_RE;
    int N1;
    int ncx1;
    int id1[N1];
    int RE_ind1;
    int colmns_nHC1;
    int<lower=0, upper=1> y1[N1];
    matrix[N1, 1] Z1;
    matrix[n, ncx1] Xhc1;
    matrix[N1, ncx1] X1;
    real<lower=0> scale_betas1;
    real<lower=0> scale_diag_D;
}
 
parameters {
    vector[ncx1] betas1;
    matrix[n, n_RE] u;
    real<lower = 0> D;
}
 
transformed parameters {
    vector[N1] eta1;
    matrix[n, n_RE] mu_u;
    for (i in 1:n) {
        mu_u[i, 1] = dot_product(Xhc1[i, 1:2:3], betas1[1:2:3]);
    }
    for (j1 in 1:N1) {
        eta1[j1] = Z1[j1, 1] * u[id1[j1], 1]
                 + X1[j1, 4] * betas1[4];
    }
}
 
model {
    D ~ cauchy(0, scale_diag_D);
    for (i in 1:n) {
        u[i, ] ~ normal(mu_u[i, ], D);
    }
    for (k1 in 1:ncx1) {
        betas1[k1] ~ normal(0.0, scale_betas1);
    }
    y1 ~ bernoulli_logit(eta1);
}
 
generated quantities {
    matrix[n, n_RE] b;
    b = u - mu_u;
}
