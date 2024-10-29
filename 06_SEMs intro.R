


library(GGally)
library(lavaan)
library(tidyverse)
library(lavaanPlot)


keeley <- read.csv("data/keeley.csv")

keeley <- keeley %>%
  select(-elev)

ggpairs(keeley)

cov(keeley)

keeley_formula = 
  'age ~ distance
  firesev ~ distance
'

summary(
  lm(age ~ firesev + distance, data = keeley)
)

keeley_sem1 <- sem(keeley_formula,
                   data = keeley,
                   meanstructure = TRUE)

summary(keeley_sem1,
        standardize = TRUE)

## full model

keeley_formula = 
  '
age ~ distance
hetero ~ distance
abiotic ~ distance
firesev ~ age
cover ~ firesev + hetero
rich ~ cover + hetero + abiotic + distance
'

keeley_sem3 <- sem(keeley_formula, 
                   data = keeley,
                   meanstructure = TRUE)

summary(keeley_sem3)


modificationIndices(keeley_sem3) %>%
  arrange(-mi) %>%
  head(5)


lavaanPlot(model = keeley_sem3, node_options = list(shape = "box", fontname = "Helvetica"), 
           edge_options = list(color = "grey"), coefs = TRUE, covs = TRUE)



### Example 2

mosquito <- read.csv("data/mosquito_nets.csv") |>
  select(-id)

ggpairs(mosquito)



