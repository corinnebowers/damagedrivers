
###################################################################################################
## CRS FOIA request results
## Created by Corinne Bowers 11/23/20
## Last updated 3/13/23
###################################################################################################

#### functions & packages #########################################################################

## packages
require(readxl)
require(dplyr)

## working directory
setwd('D:/Research/_data/losses/CRS_FOIA')

## functions
toNumber <- function(x) as.numeric(paste(x))


#### read in data #################################################################################

crs1998 <- read_xls('(1) CRS_Historical_Rating_Data_1998 to 1999.zip5.xls', sheet = 'Oct98')
crs1999 <- read_xls('(1) CRS_Historical_Rating_Data_1998 to 1999.zip5.xls', sheet = 'Oct99')

crs2000 <- read_xls('(2) CRS_Historical_Rating_Data_2000 to 2004.zip4.xls', sheet = 'Oct00')
crs2001 <- read_xls('(2) CRS_Historical_Rating_Data_2000 to 2004.zip4.xls', sheet = 'Oct01')
crs2002 <- read_xls('(2) CRS_Historical_Rating_Data_2000 to 2004.zip4.xls', sheet = 'Oct02')
crs2003 <- read_xls('(2) CRS_Historical_Rating_Data_2000 to 2004.zip4.xls', sheet = 'Oct03')
crs2004 <- read_xls('(2) CRS_Historical_Rating_Data_2000 to 2004.zip4.xls', sheet = 'Oct04')

crs2005 <- read_xls('(3) CRS_Historical_Rating_Data_2005 to 2009_.zip3.xls', sheet = 'Oct05')
crs2006 <- read_xls('(3) CRS_Historical_Rating_Data_2005 to 2009_.zip3.xls', sheet = 'Oct06')
crs2007 <- read_xls('(3) CRS_Historical_Rating_Data_2005 to 2009_.zip3.xls', sheet = 'Oct07')
crs2008 <- read_xls('(3) CRS_Historical_Rating_Data_2005 to 2009_.zip3.xls', sheet = 'May08')
crs2009 <- read_xls('(3) CRS_Historical_Rating_Data_2005 to 2009_.zip3.xls', sheet = 'May09')

crs2010 <- read_xlsx('(4) CRS_Historical_Rating_Data_2010 to 20152.xlsx', sheet = 'Oct10')
crs2011 <- read_xlsx('(4) CRS_Historical_Rating_Data_2010 to 20152.xlsx', sheet = 'Oct11')
crs2012 <- read_xlsx('(4) CRS_Historical_Rating_Data_2010 to 20152.xlsx', sheet = 'Oct12')
crs2013 <- read_xlsx('(4) CRS_Historical_Rating_Data_2010 to 20152.xlsx', sheet = 'Oct13')
# crs2014 <- read_xlsx('(4) CRS_Historical_Rating_Data_2010 to 20152.xlsx', sheet = 'May14')
crs2015a <- read_xlsx('(4) CRS_Historical_Rating_Data_2010 to 20152.xlsx', sheet = 'Oct15 (2013Manual)')
crs2015b <- read_xlsx('(4) CRS_Historical_Rating_Data_2010 to 20152.xlsx', sheet = 'Oct15 (2007Manual)')

crs2014a <- read_xlsx('(5) CRS_Historical_Rating_Data_2014 to 20188.xlsx', sheet = 'Oct14(2007CM)')
crs2014b <- read_xlsx('(5) CRS_Historical_Rating_Data_2014 to 20188.xlsx', sheet = 'Oct14(2013CM)')
# crs2015 <- read_xlsx('(5) CRS_Historical_Rating_Data_2014 to 20188.xlsx', sheet = 'Oct15')
crs2016a <- read_xlsx('(5) CRS_Historical_Rating_Data_2014 to 20188.xlsx', sheet = 'Oct16(2007CM)')
crs2016b <- read_xlsx('(5) CRS_Historical_Rating_Data_2014 to 20188.xlsx', sheet = 'Oct16(2013CM)')
crs2017a <- read_xlsx('(5) CRS_Historical_Rating_Data_2014 to 20188.xlsx', sheet = 'Oct17(2007CM)')
crs2017b <- read_xlsx('(5) CRS_Historical_Rating_Data_2014 to 20188.xlsx', sheet = 'Oct17(13&17CM)')
crs2018a <- read_xlsx('(5) CRS_Historical_Rating_Data_2014 to 20188.xlsx', sheet = 'Oct18 (2007Manual)')
crs2018b <- read_xlsx('(5) CRS_Historical_Rating_Data_2014 to 20188.xlsx', sheet = 'Oct18 (13&17CM)')

crs2015c <- read_xlsx('(6) 2015 thru 2019 CRS Classifications7.xlsx', skip = 2, sheet = '2015')
crs2016c <- read_xlsx('(6) 2015 thru 2019 CRS Classifications7.xlsx', skip = 2, sheet = '2016')
crs2017c <- read_xlsx('(6) 2015 thru 2019 CRS Classifications7.xlsx', skip = 2, sheet = '2017')
crs2018c <- read_xlsx('(6) 2015 thru 2019 CRS Classifications7.xlsx', skip = 2, sheet = '2018')
crs2019a <- read_xlsx('(6) 2015 thru 2019 CRS Classifications7.xlsx', skip = 2, sheet = '2019')

crs2019b <- read_xlsx('(7) CRS_Historical_Rating_Data_Oct_20196.xlsx', skip = 3, sheet = 'Oct19 (2013&2017 Manual)')


#### combine files ################################################################################

crs <- rbind(cbind(Year = 1998, crs1998), cbind(Year = 1999, crs1999))
crs <- crs %>%
  rbind(cbind(Year = 2002, crs2002)) %>%
  rbind(cbind(Year = 2003, crs2003)) %>%
  rbind(cbind(Year = 2004, crs2004)) %>%
  rbind(cbind(Year = 2005, crs2005)) %>%
  rbind(cbind(Year = 2006, crs2006)) %>%
  rbind(cbind(Year = 2007, crs2007)) %>%
  rbind(cbind(Year = 2008, crs2008)) %>%
  rbind(cbind(Year = 2009, crs2009)) %>%
  rbind(cbind(Year = 2010, crs2010) %>% rename('Community Name' = 'COMMUNITY NAME')) %>%
  rbind(cbind(Year = 2011, crs2011 %>% select(-CGA)) %>% rename('Community Name' = 'COMMUNITY')) %>%
  rbind(cbind(Year = 2012, crs2012) %>% rename('Community Name' = 'Name')) %>%
  rbind(cbind(Year = 2013, crs2013) %>% rename('Community Name' = 'Name')) %>%
  rbind(cbind(Year = 2014, rbind(crs2014a, crs2014b)) %>% rename('Community Name' = 'Name'))
crs <- rbind(crs2015a, crs2015b, by = c('CID', 'COMMUNITY NAME', 'STATE')) %>%
  setNames(c('CID', 'Community Name', 'State', 'Class', 'cTot')) %>%
  cbind(Year = 2015, .) %>%
  subset(!is.na(toNumber(CID))) %>%
  rbind(crs, .)
crs <- rbind(crs2016a %>% rename(cTOT = cTOTAL), crs2016b, by = c('CID', 'COMMUNITY NAME', 'STATE')) %>%
  setNames(c('CID', 'Community Name', 'State', 'Class', 'cTot')) %>%
  cbind(Year = 2016, .) %>%
  subset(!is.na(toNumber(CID))) %>%
  rbind(crs, .)
crs <- rbind(crs2017a %>% rename(cTOT = cTOTAL), crs2017b, by = c('CID', 'COMMUNITY NAME', 'STATE')) %>%
  setNames(c('CID', 'Community Name', 'State', 'Class', 'cTot')) %>%
  cbind(Year = 2017, .) %>%
  subset(!is.na(toNumber(CID))) %>%
  rbind(crs, .)
crs <- rbind(crs2018a, crs2018b, by = c('CID', 'COMMUNITY NAME', 'STATE')) %>%
  setNames(c('CID', 'Community Name', 'State', 'Class', 'cTot')) %>%
  cbind(Year = 2018, .) %>%
  subset(!is.na(toNumber(CID))) %>%
  rbind(crs, .)
crs <- full_join(crs2019a, crs2019b %>% select(-REGION, -CLASS, -'COMMUNITY NAME'), by = c('CID', 'STATE')) %>%
  setNames(c('CID', 'Community Name', 'State', 'Class', 'cTot')) %>%
  cbind(Year = 2019, .) %>%
  subset(!is.na(toNumber(CID))) %>%
  subset(!(CID == 135161 & State == 'PA')) %>%
  rbind(crs, .)
crs <- crs2000[,c('Community Name', 'State')] %>%
  full_join(crs[,c('Community Name', 'State', 'CID')], by = c('Community Name', 'State')) %>%
  full_join(crs2000, by = c('Community Name', 'State')) %>%
  mutate(Year = 2000) %>%
  select(match(names(crs), names(.))) %>%
  rbind(crs, .)
crs <- crs2001[,c('Community Name', 'State')] %>%
  full_join(crs[,c('Community Name', 'State', 'CID')], by = c('Community Name', 'State')) %>%
  full_join(crs2001, by = c('Community Name', 'State')) %>%
  mutate(Year = 2001) %>%
  select(match(names(crs), names(.))) %>%
  rbind(crs, .)
crs <- crs[order(crs$Year),]


###################################################################################################
