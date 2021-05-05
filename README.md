# Implementing a toy ERGM from scratch in R
The `myergm` package tries to estimate the `edges` and `triangle` parameters from scratch as much as possible. 
If you carefully look at the source code, you will realize that `myergm` implements an MPLE using the `ergm` package. 
This is totally attributed to the author's laziness.

This trial is motivated by the following blog posts:
- [Implementing an ERGM from scratch in Python](https://computational-communication.com/ergm-python/)
- [RでERGMを実装したかった（が、失敗した）(Wanted to implement an ERGM in R, but failed)](http://meana0.hatenablog.com/entry/2019/12/21/120043)

## Example
```r
# Install the package.
devtools::install_github("komatsuna4747/myergm@main", subdir = "myergm")

# Load the flomarriage network data.
library(ergm)
data(florentine)

# Estimate the parameters.
theta <- myergm::myergm_MCMLE(model = flomarriage ~ edges + triangle,
                              seed = 334)
                              
# Estimate the parameters by the original package `ergm`.
original_ergm <- ergm(formula = flomarriage ~ edges + triangle, 
                      control = control.ergm(seed = 334))

theta_ergm <- original_ergm$coef

# Compare the estimates.
data.frame("myergm" = theta_myergm,
           "ergm" = theta_ergm)
```

```
             myergm       ergm
edges    -1.6781074 -1.6853490
triangle  0.1548246  0.1600625
```

The estimates by `myergm` seem very close to those by `ergm`, although the computation time of `myergm`is far behind that of `ergm`.
This point needs to be further investigated.

## If you can't install `myergm`...
It would be recommended to use Docker to replicate the example above. Provided that Docker is installed, you can create a Docker container for the replication as follows.

1. Change `.env.example` to `.env` and modify the environment variables as appropriate. You can check USERID (uid) and GROUPID (gid) by typing `id` on Linux.
1. Build a Docker image by `docker-compose build`. This will take a few minutes. Note that `.env`, `Dockerfile`, and `docker-compose.yml` must be all in the current directory.
1. Start a Docker container by `docker-compose up -d`.
1. Enter `http://localhost:8787` on your browser. Then RStudio will be launched. 
