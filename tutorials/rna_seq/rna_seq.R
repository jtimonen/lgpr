# Longitudinal RNA-seq data analysis
require(lgpr)
require(ggplot2)

# Setup
N_ITER <- 3000
ADAPT_DELTA <- 0.95
I_GENE <- 425

#Function for reading data
source('create_df.R')

# Control parameters
CNTR  <- list(adapt_delta = ADAPT_DELTA)

# Read data frame from file
dat   <- readRDS(file = 'data_519.rds')
cat(paste0('gene index = ', I_GENE, '\n'))
input <- create_df(I_GENE, dat)
data  <- input$df
gnam  <- dat$gene_names[I_GENE]
cat(paste0('gene_name: ', gnam, '\n'))

# Set normalization
c_hat <- adjusted_c_hat(data$y, input$norm_fac)

# Fit model
my_prior <- list(wrp = igam(shape = 20, scale = 16))

fit <- lgp(formula = y ~ gp(age)*zs(id) + gp(age) + gp_vm(dis_age) + 
             zs(group) + gp(age)*zs(sex),
           data = data,
           prior = my_prior,
           likelihood = "nb",
           c_hat = c_hat,
           iter = N_ITER,
           chains = 4,
           cores = 4,
           control = CNTR,
           save_warmup = FALSE,
           verbose = TRUE)

# Plots
plot_1 <- plot_pred(fit)
draws <- sample.int(1000, 50)
t <- seq(2, 50, by = 2)
x_pred <- new_x(data, t, x_ns = "dis_age")
c_hat_pred <- rep(6.0, nrow(x_pred))
p <- pred(fit, x_pred, draws = draws,
          c_hat_pred = c_hat_pred)
plot_2 <- plot_pred(fit, x = x_pred, pred = p, alpha = 0.2)

cb <- c(NA, NA, "group", "group", "sex", "group")
plot_3 <- plot_components(fit, x = x_pred, pred = p,
                          color_by = cb)
