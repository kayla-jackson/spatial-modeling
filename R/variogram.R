# compute max distance

sfe_sub <- sfe[, colData(sfe)$cell_type == ct]
sfe_sub <- sfe_sub[rowSums(counts(sfe_sub)) > 500, ]

ct_sf <- colGeometry(sfe_sub, 'centroids')
dat <- as(t(counts(sfe_sub)), "matrix")

ct_sf  %<>% mutate(data.frame(dat))
min_ct_sf <- ct_sf[,c(1,3)] %>% mutate(row = row_number())

dist2 <- sf::st_distance(min_ct_sf)


D <- dist2 %>% as.vector() %>% max()

# h < D/ 2, distance of reliability 
D/2

# only consider for bins with > 30 measurements, use quantiles to start? 
bins <- dist2 %>% as.vector() %>% unique()%>% quantile(seq(0,1,0.05))
bins <- as.integer(bins)
nb <- length(bins)

intervals <- purrr::map2(bins[-nb], bins[-1], c)

ng <- purrr::map(intervals, ~sum(dist2 > .x[1] & dist2 <= .x[2]))

pairs_ <- purrr::map(intervals, ~ which(dist2 > .x[1] & dist2 <= .x[2], arr.ind = TRUE))

tic()
pairs_rc <- purrr:::map(pairs_, function(x) 
  cbind(x, x=pull(ct_sf, gene)[x[,"row"]],y=pull(ct_sf, gene)[x[,"col"]])
)

pairs_vgm <- purrr::map_dbl(pairs_rc, function(x) mean((x[,3]-x[,4])^2)/2)

toc()
#library(furrr)
# library(tictoc)

tic()
ll <- list()
for (gene in names(ct_sf)[-1][1:10]) {
  plan(sequential)
  
  ll[[gene]] <- furrr:::future_map_dbl(pairs_, function(x) 
    mean((pull(ct_sf, gene)[x[,"row"]]-pull(ct_sf, gene)[x[,"col"]])^2)/2
  )
}
toc()


library(purrr)
ll <- list()
for (gene in names(ct_sf)[-1]) {
  print(gene)
  pairs_rc <- purrr:::map(pairs_, function(x) 
    cbind(x, x=pull(ct_sf, gene)[x[,"row"]],y=pull(ct_sf, gene)[x[,"col"]])
  )

ll %>% walk(function(x){ 
    tiff(glue::glue("plot{rnorm(1)}.tiff"))
    plot(y = c(0,x), x = bins) 
    dev.off()
})
  ll[[gene]] <- purrr::map_dbl(pairs_rc, function(x) mean((x[,3]-x[,4])^2)/2)
  
  #ll[[gene]] <- pairs_vgm
}

