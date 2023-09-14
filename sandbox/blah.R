quantile_df <- function(x, probs = c(0.025, 0.5, 0.975, 0.99)) {
  tibble(
    value = quantile(x, probs, na.rm = TRUE),
    prob = probs
  )
}


dfc$id <- 1:100

ggplot(dfc, aes(label = id)) + 
  geom_sf() + 
  geom_sf_text()

borders <- c(1:10, 91:100, seq(11,91,10), seq(20,100,10))

dfc_empty <- dfc[borders,] %>% distinct(cellid, .keep_all=TRUE)

# generate connectivity matrix
nb_rook_empty = st_rook(dfc_empty)

nb_rook_empty = as.nb.sgbp(nb_rook_empty)
col.w.empty <- nb2listw(nb_rook_empty)

W.empty <- nb2mat(nb_rook_empty, style="B")
D.empty <- Diagonal(nrow(dfc_empty), rowSums(W))

# Find edges in W
nodes <- which(W==1, arr.ind=TRUE)
nodes <- nodes[nodes[,1] < nodes[,2],]

yc <- rpois(nrow(dfc_empty), 10*dfc_empty$theta)

# Set up data for STAN
data <- list(
  N = nrow(dfc_empty), 
  N_edges = nrow(nodes), 
  #node1 = nodes[,1],
  #node2 = nodes[,2],
  D=as.matrix(D.empty),
  W = as.matrix(W.empty),
  mu = 10,
  counts = yc
)

pois.file <- path_wd("stan/pois",ext="stan")

# Compile and run
pois.mod <- cmdstan_model(pois.file)
fits <- pois.mod$sample(data, chains=2)

draws <- as_draws_matrix(fit.pois)
yrep <- draws[, glue::glue('counts_ppc[{1:36}]')]

ppc_dens_overlay(yc, yrep[1:200, ])

dfc_empty <- dfc[borders,]
# Prepare data for plotting probability of capture draws
plot_df <- map2(dfc_empty$theta, dfc_empty$geometry, \(x,y) c(p = x, x=st_centroid(y)[1], y=st_centroid(y)[2]))
plot_df <- bind_rows(plot_df)
prep <- draws[, glue::glue('prob_ppc[{1:36}]')]

plot_df <- bind_cols(plot_df, t(prep)) %>% mutate(index = row_number())

plot_df <- plot_df %>% 
  rename_with(
    ~ glue("prep_{str_remove(.x, '`')}"),
    matches("[^a-z]")
  )

plot_df <- plot_df %>% 
  pivot_longer(starts_with("p")) %>%
  mutate(fill = ifelse(str_detect(name, 'rep'), "yrep", "y")) %>%
  arrange(name,x) 

cis <- plot_df %>% 
  filter(name != "p") %>%
  reframe(quantile_df(value), .by = index)
  
# 95% CI of draws plotted with true value
p3 <- ggplot(mapping=aes(index, value)) + 
  geom_line(
    data=plot_df %>% filter(fill=='yrep', name %in% paste0("prep_",sample(30))),
    mapping=aes(color=name),
    alpha=0.5,
    #color='#deebf7'
    ) +
  geom_line(
    data=plot_df %>% filter(fill !='yrep'),
    color='#3182bd'
   ) +
  # geom_line(
  #    data = cis %>% filter(prob %in% c(0.025)),
  #    color = 'red'
  # ) +
  # geom_line(
  #   data = cis %>% filter(prob %in% c(0.975)),
  #   color = 'red'
  # ) +
  labs(
    y = "Probability of capture (p)"
  ) + 
  
  guides(color="none") + 
  scale_color_manual(labels=glue("prep_{1:2000}"), values=rep("#deebf7",30)) +
  theme_light() + 
  theme(
    aspect.ratio=1/2,
    panel.grid=element_blank()
  )

ggsave('prob5.png', p3, width=2,height=1.25, scale=2)

hyperpars <- fits$summary(c("car_rho", 'car_scale'))$mean
names(hyperpars) <- c("car_rho", 'car_scale')
hyperpars

library(posterior)
library(bayesplot)

mcmc_hist(fits$draws("car_scale"))
mcmc_intervals(fits$draws(), pars = c("car_scale", "car_rho"))

mcmc_trace(fits$draws(), pars = c('car_scale'))

# Testing Geostan
# library(geostan)
# 
# cp <- prep_car_data(W)
# 
# fit <- stan_car(counts ~ offset(log(rep(10,36))),
#                 car_parts = cp,
#                 data = data,
#                 family = poisson(),
#                 iter = 2000, chains = 1 # for example speed only
# )
