require(covr)
exclusions <- list(
  'src/stanExports_lgp_marginal.h',
  'src/stanExports_lgp_latent.h',
  'src/stanFunctions.cpp'
)
cov <- package_coverage(line_exclusions = exclusions, pre_clean = TRUE)
report(cov)

# upload to codecov.io
#codecov(coverage = cov, token = tok)