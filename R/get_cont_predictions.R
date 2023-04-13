library(wqtrends)
library(tidyverse)

# download and load data file with all parameters, loaded as datprc object
fl <- paste0(tempdir(), '/datproc.RData')
download.file('https://github.com/tbep-tech/wqtrends-synthesis/raw/main/data/datprc.RData', destfile = fl)
load(file = fl)

# parameter to model, one of chl, do, dosat, gpp, or kd
prm <- 'chl'

# location to save file
outdir <- '~/Desktop'

# data to model, same as datprc, params in wide format, nested by station
# transformation assigned to parameter
tomod <- datprc %>% 
  filter(param %in% prm) %>% 
  mutate(
    trans = case_when(
      param %in% c('chl', 'gpp', 'kd') ~ 'log10', 
      T ~ 'ident'
    )
  ) %>% 
  nest(.by = c('station', 'trans'))

# create models for every station for selected parameter
# ignore the warnings
modssta <- tomod %>%
  mutate(
    modi = purrr::pmap(list(station, trans, data), function(station, trans, data){
      
      cat(station, '\n')
      out <- anlz_gam(data, trans = trans)
      return(out)
      
    })
  )

# get daily model predictions and CI
prds <- modssta %>% 
  mutate(
    preds = purrr::map(modi, function(x){

      prddat <- x$model
      trans <- x$trans
    
      # prediction matrix for daily values
      rng <- prddat$cont_year %>%
        range(na.rm  = T) %>% 
        date_decimal() %>% 
        as.Date
      newdat <- seq.Date(floor_date(rng[1], 'year'), ceiling_date(rng[2], 'year'), by = 'days') %>%
        tibble(date = .) %>%
        mutate(
          doy = yday(date),
          cont_year = decimal_date(date),
          yr = year(date)
        )
      
      # get predictions and se, calc CI and backtransform
      prd <- predict(x, newdata = newdat, se.fit = T)
      prd <- tibble(fit = prd$fit, se = prd$se.fit)
      
      if(trans == 'log10')
        out <- prd %>% 
          mutate(
            hiv = 10 ^ (fit + 1.96 * se), 
            lov = 10 ^ (fit - 1.96 * se), 
            val = 10 ^ fit
          ) %>% 
          select(val, hiv, lov)
        
      if(trans == 'ident')
        out <- prd %>% 
          mutate(
            hiv = fit + 1.96 * se, 
            lov = fit - 1.96 * se, 
            val = fit
          ) %>% 
          select(val, hiv, lov)
        
      out <- bind_cols(newdat, out)
      
      return(out)
      
    }) 
    
  ) %>% 
  select(station, preds) %>% 
  unnest(preds)

# save file to location at outdir
flnm <- paste0(outdir, '/preds', prm, '.csv')
write.csv(prds, file = flnm, row.names = F)

##
# make a plot for run

# observed data
pts <- tomod %>% 
  unnest('data')

# axis transformation
trans <- ifelse(unique(tomod$trans) == 'ident', 'identity', 'log10')

ggplot(prds, aes(x = date)) + 
  geom_ribbon(aes(ymin = lov, ymax = hiv), fill = 'grey') +
  geom_line(aes(y = val), col = 'blue') + 
  geom_point(data = pts, aes(y = value), size = 0.1) + 
  scale_y_continuous(trans = trans) + 
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
  ) + 
  facet_wrap(~station, ncol = 2, scales = 'free_y') + 
  labs(
    y = prm,
    x = NULL
  )
  
