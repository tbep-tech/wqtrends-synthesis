library(tidyverse)
library(lubridate)
library(sf)
library(wqtrends)
library(googlesheets4)

# format raw wq data for use with wqtrends --------------------------------

# all raw files from DS
chlraw <- read.csv('data/raw/sfb_surf_CB_SB_LSB.csv', stringsAsFactors = F)
gppraw <- read.csv('data/raw/sfb_GPP_monthly.csv', stringsAsFactors = F) 
doraw <- read.csv('data/raw/CB_SB_LSB_depthavg_O2.csv', stringsAsFactors = F)

# get chlorophyll 
chldat <- chlraw %>% 
  select(date, station, chl = chl_merge) %>% 
  gather('param', 'value', -date, -station) %>% 
  mutate(
    date = ymd(date),
    doy = yday(date), 
    cont_year = decimal_date(date),
    yr = year(date),
    mo = month(date, label = T),
    param = tolower(param)
  ) %>% 
  filter(yr >= 1990 & yr <= 2019) %>% 
  filter(!is.na(value))

# get gpp
gppdat <- gppraw %>% 
 select(cont_year = dec_date, station, GPP) %>% 
 mutate(
   date = date_decimal(cont_year),
   date = as.Date(date),
   doy = yday(date),
   yr = year(date),
   mo = month(date, label = T)
 ) %>% 
 rename(
   gpp = GPP
 ) %>% 
 gather('param', 'value', gpp) %>% 
 select(date, station, param, value, doy, cont_year, yr, mo) %>% 
 filter(yr >= 1990 & yr <= 2019) %>% 
 filter(!is.na(value))

# get depht-averaged do
dodat <- doraw %>% 
  select(date, station, do, dosat = do_sat) %>% 
  gather('param', 'value', -date, -station) %>% 
  mutate(
    date = ymd(date),
    doy = yday(date), 
    cont_year = decimal_date(date),
    yr = year(date),
    mo = month(date, label = T),
    param = tolower(param)
  ) %>% 
  filter(yr >= 1990 & yr <= 2019) %>% 
  filter(!is.na(value))

# combine new do ests, gpp with datprc
datprc <- bind_rows(chldat, gppdat, dodat) %>% 
 arrange(station, param, date)

save(datprc, file = 'data/datprc.RData', compress = 'xz')

# save separate rdata model files for each station, parameter -------------

data(datprc)

# data to model, same as datprc, params in wide format, nested by station
# crossed with frms
tomod <- datprc %>% 
  filter(param %in% c('chl', 'do', 'dosat', 'gpp')) %>% # add parameters here
  group_by(station, param) %>% 
  nest %>% 
  mutate(
    trans = case_when(
      param %in% c('chl', 'gpp') ~ 'log10', 
      T ~ 'ident'
    )
  )

# create models for every station, gam model eval
modssta <- tomod %>%
  mutate(
    modi = purrr::pmap(list(station, param, trans, data), function(station, param, trans, data){
      
      cat(station, param, '\n')
      out <- anlz_gam(data, trans = trans)
      return(out)
      
    })
  )

# separate models into diff files by parameter and stations (single is too large for git)
tosv <- modssta %>% 
  select(station, param) %>% 
  unique

for(i in 1:nrow(tosv)){
  
  cat(i, 'of', nrow(tosv), '\n')
  
  sta <- tosv[[i, 'station']]
  param <- tosv[[i, 'param']]
  
  fl <- modssta %>% 
    filter(station %in% !!sta) %>% 
    filter(param %in% !!param)
  
  flnm <- paste0('mods_', param, sta)
  
  assign(flnm, fl)
  
  save(list = flnm, file = paste0('data/', flnm, '.RData'), compress = 'xz')
  
}

# get model estimates for synthesis indicators ----------------------------

# # file with lookup combos
# # https://docs.google.com/spreadsheets/d/13qYK33QKpfm6_jVKXB8J_T3Q8K8ut8lZUM1CTzXqQ5o/edit#gid=0
# sht <- read_sheet('13qYK33QKpfm6_jVKXB8J_T3Q8K8ut8lZUM1CTzXqQ5o') %>% 
#   select(var, metfun, doystr, doyend)

# metric lookup table
mets <- tibble(
    param = c("chl", "gpp", "dosat", "do", "gpp", "chl", "chl", "chl", "gpp", "gpp"), 
    metfun = c("mean", "mean", "mean", "mean", "max", "max", "var", "var", "var", "var"), 
    doystr = c(32, 32, 182, 182, 32, 32, 32, 152, 32, 152), 
    doyend = c(121, 152, 305, 305, 152, 152, 152, 274, 152, 274)
  )

# estimates, takes about four hours to run
metests <- list.files('data', pattern = '^mods', full.names = T) %>% 
  tibble(
    fl = .,
    station = gsub('\\D+', '', .),
    param = gsub('^data/mods_|[0-9]*|\\.RData$', '', .)
  ) %>% 
  left_join(mets, by = 'param') %>% 
  group_by(fl) %>% 
  nest() %>% 
  mutate(
    mod = purrr::map(fl, function(x){
      
      load(file = x)
      
      nm <- basename(x)
      nm <- gsub('\\.RData', '', nm)
      
      out <- get(nm) %>% 
        pull(modi) %>% 
        deframe()
      
      return(out)
      
    })
  ) %>% 
  unnest(c('data')) %>% 
  group_by(metfun, doystr, doyend, mod, fl) %>% 
  nest() %>% 
  mutate(
    metseas = purrr::pmap(list(fl, metfun = metfun, mod = mod, doystr = doystr, doyend = doyend), function(fl, metfun, mod, doystr, doyend){
      
      cat(fl, metfun, '\n')
      
      toprs <- paste0('anlz_metseason(mod, ', metfun, ', ', doystr, ', ', doyend, ', nsim = 1e5, na.rm = TRUE)')
      out <- eval(parse(text = toprs))
    
      return(out)
      
    })
  ) %>% 
  unnest('data') %>% 
  ungroup() %>% 
  select(-mod) %>% 
  unnest('metseas') %>% 
  select(station, param, metfun, doystr, doyend, yr, upr = bt_upr, lwr = bt_lwr, met = bt_met)

write.csv(metests, 'results/metests.csv', row.names = F)
