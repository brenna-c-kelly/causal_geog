
library(sf)
library(see)
library(ape)
library(tmap)
library(mgcv)
library(lme4)
library(nlme)
library(ggdag)
library(dagitty)
library(tidyverse)
library(correlation)
library(piecewiseSEM)

keeley <- read.csv("data/keeley.csv")

cor_res <- correlation(keeley)
x <- cor_sort(as.matrix(cor_res))
plot(visualisation_recipe(x))


## models

# saturated model
mod1 <- lm(firesev ~ age, data = keeley)
mod2 <- lm(cover ~ age + firesev, data = keeley)

summary(mod1)
summary(mod2)

keeley_psem1 <- psem(
  mod1, mod2,
  data = keeley
)

summary(keeley_psem1)

# simpler model
mod1 <- lm(firesev ~ age, data = keeley)
mod2 <- lm(cover ~ firesev, data = keeley)

summary(mod1)
summary(mod2)

keeley_psem2 <- psem(
  mod1, mod2,
  data = keeley
)

summary(keeley_psem2)

plot(keeley_psem2, show = "unstd")

# full grace & keeley model
keeley_psem3 <- psem(
  lm(hetero ~ distance, data = keeley),
  lm(abiotic ~ distance, data = keeley),
  lm(age ~ distance, data = keeley),
  lm(firesev ~ age, data = keeley),
  lm(cover ~ firesev, data = keeley),
  lm(rich ~ distance + abiotic + hetero + cover, data = keeley)
)

summary(keeley_psem3, .progressBAr = FALSE)

plot(keeley_psem3, show = "std", ns_dashed = TRUE)

## gam model
mod1 <- gam(firesev ~ s(age), data = keeley)
mod2 <- lm(cover ~ firesev, data = keeley)

keeley_psem4 <- psem(
  mod1, mod2,
  data = keeley
)
summary(keeley_psem4, .progressBar = FALSE)

plot(keeley_psem4[[1]])

AIC(keeley_psem2, keeley_psem4)

# sythetic dataset

shipley_dag <- dagify(DD ~ lat,
                      date ~ DD,
                      growth ~ date,
                      live ~ growth,
                      coords = list(x = c(lat = 1,
                                          DD = 2,
                                          date = 3,
                                          growth = 4,
                                          live = 5), 
                                    y = c(lat = 1,
                                          DD = 2,
                                          date = 3,
                                          growth = 2,
                                          live = 1)
                      ),
                      labels = c(lat = "Latitude",
                                 DD = "Degree Day",
                                 date = "Date",
                                 growth = "Growth",
                                 live = "Live")
)

ggdag(shipley_dag, use_labels = "label",
      text = FALSE) +
  theme_dag()


shipley |>
  mutate(tree = as.factor(tree),
         site = as.factor(site)) |>
  ggplot(aes(x = year, y = Growth, col = tree)) +
  geom_line() +
  facet_wrap(~ site) +
  theme_bw() +
  theme(legend.position = "none")

# build models
mod1 <- lmer(DD ~ lat + (1 | tree) + (1 | site), data = shipley)
mod2 <- lmer(Date ~ DD + (1 | tree) + (1 | site), data = shipley)
mod3 <- lmer(Growth ~ Date + (1 | tree) + (1 | site), data = shipley)

mod4 <- glmer(Live ~ Growth + (1 | tree) + (1 | site),
              family = binomial(link = "logit"), data = shipley)

shipley_psem1 <- psem(
  mod1, mod2, mod3, mod4,
  data = shipley
)

summary(shipley_psem1, .progressBar = FALSE)


# spatial models

boreal <- read.csv("data/boreal.csv")
head(boreal)

boreal_sf <- st_as_sf(boreal, coords = c("x", "y"))

tm_shape(boreal_sf) +
  tm_symbols(col = "NDVI", size = 0.75, alpha = 0.75, 
             palette = "Greens", style = "fisher") +
  tm_layout(legend.outside = TRUE)

mod1 <- lm(richness ~ temp, data = boreal)
mod2 <- lm(NDVI ~ richness + temp + wet, data = boreal)

boreal_psem1 <- psem(
  mod1, mod2,
  data = boreal
)
summary(boreal_psem1)


boreal_sf$richness_res <- residuals(boreal_psem1)[,1]
boreal_sf$ndvi_res <- residuals(boreal_psem1)[,2]

m1 <- tm_shape(boreal_sf) +
  tm_symbols(col = "richness_res", size = 0.75, alpha = 0.75,
             style = "fisher") +
  tm_layout(legend.outside = TRUE)
m2 <- tm_shape(boreal_sf) +
  tm_symbols(col = "ndvi_res", size = 0.75, alpha = 0.75,
             style = "fisher") +
  tm_layout(legend.outside = TRUE)
tmap_arrange(m1, m2)

distMat <- as.matrix(dist(
  cbind(boreal$x, boreal$y))
)

distsInv <- 1/distMat
diag(distsInv) <- 0

Moran.I(boreal_sf$richness_res, distsInv)

Moran.I(boreal_sf$ndvi_res, distsInv)


richness_gls <- gls(richness ~ temp,
                    correlation = spaceCor,
                    data = boreal)
ndvi_gls<- gls(NDVI ~ richness + temp + wet,
               correlation = spaceCor,
               data=boreal)




