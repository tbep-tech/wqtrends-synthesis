library(tidyverse)
library(wqtrends)
library(patchwork)

data(datprcorig) # datprc created from branch original-synthesis
data(datprc)

# station 34 --------------------------------------------------------------

datprcold <- datprcorig %>%   
  filter(param == 'gpp') %>% 
  filter(station == 34) 
datprcnew <- datprc %>%   
  filter(param == 'gpp') %>% 
  filter(station == 34) 

# fit model with wqtrends
modold <- anlz_gam(datprcold, trans = 'log10')
modnew <- anlz_gam(datprcnew, trans = 'log10')

lms <- as.Date(c('1990-01-01', '2020-01-01'))
p1 <- show_prdseries(modold, ylab = 'GPP') + labs(title= 'Station 34', subtitle = 'old model') + scale_x_date(limits = lms)
p2 <- show_prdseries(modnew, ylab = 'GPP') + labs(subtitle = 'new model') + scale_x_date(limits = lms) 

p1 + p2 + plot_layout(ncol = 1)

# station 22 --------------------------------------------------------------

datprcold <- datprcorig %>%   
  filter(param == 'gpp') %>% 
  filter(station == 22) 
datprcnew <- datprc %>%   
  filter(param == 'gpp') %>% 
  filter(station == 22) 

# fit model with wqtrends
modold <- anlz_gam(datprcold, trans = 'log10')
modnew <- anlz_gam(datprcnew, trans = 'log10')

lms <- as.Date(c('1983-01-01', '2020-01-01'))
p1 <- show_prdseries(modold, ylab = 'GPP') + labs(title= 'Station 22', subtitle = 'old model') + scale_x_date(limits = lms)
p2 <- show_prdseries(modnew, ylab = 'GPP') + labs(subtitle = 'new model') + scale_x_date(limits = lms)

p1 + p2 + plot_layout(ncol = 1)
