# load packges
library(tidyverse)
library(lubridate)

# ── Attaching packages ───────────────────────────────────────────────────────── tidyverse 1.2.1 ──
# ✔ ggplot2 3.0.0     ✔ purrr   0.2.5
# ✔ tibble  1.4.2     ✔ dplyr   0.7.6
# ✔ tidyr   0.8.1     ✔ stringr 1.3.1
# ✔ readr   1.1.1     ✔ forcats 0.3.0
# ── Conflicts ──────────────────────────────────────────────────────────── tidyverse_conflicts() ──
# ✖ dplyr::filter() masks stats::filter()
# ✖ dplyr::lag()    masks stats::lag()

# load data
pit = read_csv("data-raw/myvatn_pitfalls.csv")
inf = read_csv("data-raw/myvatn_infalls.csv")

# clean piftall data
pit_clean = pit %>%
    # only include data from lakeid = "myv" (the other lake has gaps) and dist = 200 (which is an error)
    # remove trans "kal2" which is an error
    filter(lakeid == "myv", dist != 200, trans != "kal2") %>%
    # create variable for year
    # create variables for aggregated abundances of certain groups
    mutate(year = year(coldate),
           sheet = aran_other + liny) %>%
    # rename 'acar_total" as "acar" and "lyco_total" and "lyco
    rename(acar = acar_total,
           lyco = lyco_total) %>%
    # select columns to keep
    select(trans, dist, year, lyco, sheet, gnap, opil, acar, stap, cara) %>%
    # convert to long form
    # keep specific taxa that make most sense for analysis
    gather(taxon, count, lyco, sheet, gnap, opil, acar, stap, cara) %>%
    # calculate total year abundance for each taxon and site
    # a few "counts" are averages, so round
    group_by(trans, dist, year, taxon) %>%
    summarize(count = sum(round(count), na.rm=T)) %>%
    # sort by year within site and taxon
    arrange(trans, dist, taxon, year)

# clean infall data
inf_clean = inf %>%
    # only include data from lakeid = "myv" (the other lake has gaps)
    filter(lakeid=="myv") %>%
    # create variable for year
    # correct large (bgch) and small (smch) midges for subsampling; calculate total midge abundance
    mutate(year = year(coldate),
           bgch = bgch/fract_bgch,
           smch = smch/fract_smch,
           midges = bgch + smch) %>%
    # calculate total year abundance for each taxon and site
    group_by(trans, dist, year) %>%
    summarize(midges = sum(midges, na.rm=T))

# merge data frames and save
myv_arth = left_join(pit_clean, inf_clean) %>% as.data.frame()
usethis::use_data(myv_arth)
