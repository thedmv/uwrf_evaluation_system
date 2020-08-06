---
Title: "Product01-Forecast-Hour-Based Evaluation"
Author: "David Melecio-Vazquez"
Date: "2/18/2020"
Output: html_document
---

<style>

table, td, th {
  border: none;
  padding-left: 1em;
  padding-right: 1em;
  min-width: 50%;
  margin-left: auto;
  margin-right: auto;
  margin-top: 1em;
  margin-bottom: 1em;
}

</style>

```{r Initialization, echo = FALSE}
# Introduction  
#In this document I will explore how to create the first part of the evaluation system I proposed. The working title of this is the "Forecast-Hour Evaluation." The idea here is that we are looking at the performance of the model by looking at how it performed with different start times (using the most recent 00-hr forecast as input).  

rm(list = ls())
# Libraries to be used are loaded.
source("./scripts/call_libraries.R")
# Functions to use 
source("./scripts/ccnymwr_functions.R")
source("./scripts/ccny_windlidar_functions.R")
source("./scripts/wrf_functions.R")

# Functions for Evaluation System
source("./scripts/evalsys_functions.R")
```


```{r Read WRF Data, echo = FALSE}
# Read Model and Observation Data  
## Read WRF Data  
#For this evaluation system we need to look at three different output folders. Here we use the folders named, `forecast_day_minus_0`, `forecast_day_minus_1`, `forecast_day_minus_2`. The contents of each of these folders will be similar: wrfout files for 86 forecast hours and time-series data for different locations of interest. Here we will first read the forecast data.  

# Set the forecast start times for each time-series output
#fcst_start0 = "2020-01-02 00:00:00"
fcst_start0 = "CHANGE_DATE_HERE"
fcst_start1 = as.character(ymd_hms(fcst_start0) - days(1))
fcst_start2 = as.character(ymd_hms(fcst_start0) - days(2))

# Read the TS data and then calculate the wind speed and wind direction from the u,v model data
tsjfkd0 = read.wrf_tslistTS("./forecast_day_minus_0_cuerg/kjfk.d03.TS", fcst_start0) %>% ts_wswd()
tsjfkd1 = read.wrf_tslistTS("./forecast_day_minus_1_cuerg/kjfk.d03.TS", fcst_start1) %>% ts_wswd()
tsjfkd2 = read.wrf_tslistTS("./forecast_day_minus_2_cuerg/kjfk.d03.TS", fcst_start2) %>% ts_wswd()

tslgad0 = read.wrf_tslistTS("./forecast_day_minus_0_cuerg/klga.d03.TS", fcst_start0) %>% ts_wswd()
tslgad1 = read.wrf_tslistTS("./forecast_day_minus_1_cuerg/klga.d03.TS", fcst_start1) %>% ts_wswd()
tslgad2 = read.wrf_tslistTS("./forecast_day_minus_2_cuerg/klga.d03.TS", fcst_start2) %>% ts_wswd()

tsnycd0 = read.wrf_tslistTS("./forecast_day_minus_0_cuerg/knyc.d03.TS", fcst_start0) %>% ts_wswd()
tsnycd1 = read.wrf_tslistTS("./forecast_day_minus_1_cuerg/knyc.d03.TS", fcst_start1) %>% ts_wswd()
tsnycd2 = read.wrf_tslistTS("./forecast_day_minus_2_cuerg/knyc.d03.TS", fcst_start2) %>% ts_wswd()

# Change column names and add Station column
#tscolnames_old = c("Date.Time", "year", "mon", "day", "hour", "min", "sec", "t", "q", "u", "v", "wspd", "wdir")
#tscolnames_new = c("Date.Time", "year", "mon", "day", "hour", "min", "sec", "Temperature", "Mixing.Ratio", "U_WIND", "V_WIND", "Wind.Speed", "Wind.Direction")
tsjfkd0 = tsjfkd0 %>% rename(Temperature = "t", Mixing.Ratio = "q", U_WIND = "u", V_WIND = "v", Wind.Speed = "wspd", Wind.Direction = "wdir") %>% mutate(Station = "JFK")
tsjfkd1 = tsjfkd1 %>% rename(Temperature = "t", Mixing.Ratio = "q", U_WIND = "u", V_WIND = "v", Wind.Speed = "wspd", Wind.Direction = "wdir") %>% mutate(Station = "JFK")
tsjfkd2 = tsjfkd2 %>% rename(Temperature = "t", Mixing.Ratio = "q", U_WIND = "u", V_WIND = "v", Wind.Speed = "wspd", Wind.Direction = "wdir") %>% mutate(Station = "JFK")
tslgad0 = tslgad0 %>% rename(Temperature = "t", Mixing.Ratio = "q", U_WIND = "u", V_WIND = "v", Wind.Speed = "wspd", Wind.Direction = "wdir") %>% mutate(Station = "LGA")
tslgad1 = tslgad1 %>% rename(Temperature = "t", Mixing.Ratio = "q", U_WIND = "u", V_WIND = "v", Wind.Speed = "wspd", Wind.Direction = "wdir") %>% mutate(Station = "LGA")
tslgad2 = tslgad2 %>% rename(Temperature = "t", Mixing.Ratio = "q", U_WIND = "u", V_WIND = "v", Wind.Speed = "wspd", Wind.Direction = "wdir") %>% mutate(Station = "LGA")
tsnycd0 = tsnycd0 %>% rename(Temperature = "t", Mixing.Ratio = "q", U_WIND = "u", V_WIND = "v", Wind.Speed = "wspd", Wind.Direction = "wdir") %>% mutate(Station = "NYC")
tsnycd1 = tsnycd1 %>% rename(Temperature = "t", Mixing.Ratio = "q", U_WIND = "u", V_WIND = "v", Wind.Speed = "wspd", Wind.Direction = "wdir") %>% mutate(Station = "NYC")
tsnycd2 = tsnycd2 %>% rename(Temperature = "t", Mixing.Ratio = "q", U_WIND = "u", V_WIND = "v", Wind.Speed = "wspd", Wind.Direction = "wdir") %>% mutate(Station = "NYC")

# Combine all WRF data to a single data frame per init forecast
wrfd0.df = rbind(tsjfkd0, tslgad0, tsnycd0)
wrfd1.df = rbind(tsjfkd1, tslgad1, tsnycd1)
wrfd2.df = rbind(tsjfkd2, tslgad2, tsnycd2)
```

```{r Read OBS data, echo = FALSE}  
## Read the OBS Data  
#Now we will read the observation data from the ASOS stations. The script that downloads the data is in `./obs_station_day_minus_0/dl_ny_asos.py`. The lines for the dates to download need to be changed before running it. Once the files are download, the lines below reads the data and adds column names.  

# Read the ASOS data for each station
asosdates = data.frame(begin = paste0(as.character(ymd_hms(fcst_start0) - days(1), format = "%Y%m%d"), "0000"),
                       end   = paste0(as.character(ymd_hms(fcst_start0) + days(1), format = "%Y%m%d"), "0000") )
obkjfk = read.asos(paste0("./obs_station_day_minus_0/JFK_", asosdates$begin, "_", asosdates$end, ".txt"))
obklga = read.asos(paste0("./obs_station_day_minus_0/LGA_", asosdates$begin, "_", asosdates$end, ".txt"))
obknyc = read.asos(paste0("./obs_station_day_minus_0/NYC_", asosdates$begin, "_", asosdates$end, ".txt"))

# Combine into one data frame. 
asos.df = rbind(obkjfk, obklga, obknyc)

# Rename the columns
names(asos.df) = c("Station", "Date.Time", "Temperature", "Relative.Humidity", "Wind.Direction", "Wind.Speed")

# Add columns for date-times to select the date of interest using filters
asos.df = asos.df %>% dateTimeCol()
```


```{r Unit Conversion, echo = FALSE}
## Unit Conversion  

#Model and observation data do not share the same units for the same variable. For temperature, WRF is in Kelvin and ASOS is in degreesF. For winds, WRF is in m/s and ASOS is in knots. The formulas used to convert the numbers to a common system is shown here. For temperature I will use Kelvin, and m/s for wind speeds.

# WRF Temperature conversions
# wrfd0.df$Temperature = wrfd0.df$Temperature - 273.15
# wrfd1.df$Temperature = wrfd1.df$Temperature - 273.15
# wrfd2.df$Temperature = wrfd2.df$Temperature - 273.15

# ASOS Temperature conversions
asos.df$Temperature = (asos.df$Temperature - 32) * 5/9 + 273.15

# ASOS Wind speed conversions
asos.df$Wind.Speed = asos.df$Wind.Speed/1.944
```

```{r Combined DFs, echo = FALSE, eval = FALSE}
## Combined Data Frames

Now we have one data frame for all the observations, and three (3) data frames of the WRF data (one data frame per forecast init time). The lines below provide a visual of the data frames.

# WRF Data
head(wrfd0.df)
head(wrfd1.df)
head(wrfd2.df)

# Observations from ASOS
head(asos.df)
```

```{r Plot Image Locations, echo = FALSE}
## Locations for Plots

# Save the output plots
knitr::opts_chunk$set(fig.path = paste0("./", fcst_start0 %>% str_replace(" 00:00:00", ""), "/") , dev = "png")

```


```{r Filter for Day of Interest, echo = FALSE}
# Match Times for Model and Observations  

## Select Day of Interest  

# Time-matching is performed using a routine that can be found in `Analysis01-Time_Matching_Problem.Rmd`. The time matching will be done per variable. For the Forecast-Hour Evaluation product, we will focus on the temperature, wind speed and wind direction variables. Also, now that we have read all the TS data and ASOS data, we need to extract the day of interest, or `doi` for the time-series.  

# Note that for this product the "day of interest" will always be the UTC date of the day before.

# Select a day of interest
doi = data.frame(Date.Time = as.POSIXct(fcst_start0, tz = "UTC")) %>% dateTimeCol()

# Day of interest for WRF data
wrfd0.df_doi = wrfd0.df %>% filter(year == doi$year) %>% filter(mon == doi$mon) %>% filter(day == doi$day)
wrfd1.df_doi = wrfd1.df %>% filter(year == doi$year) %>% filter(mon == doi$mon) %>% filter(day == doi$day)
wrfd2.df_doi = wrfd2.df %>% filter(year == doi$year) %>% filter(mon == doi$mon) %>% filter(day == doi$day)

# Day of interest for ASOS data
asos.df_doi = asos.df %>% filter(year == doi$year) %>% filter(mon == doi$mon) %>% filter(day == doi$day)
```
  
```{r Comments, echo = FALSE}
# We now have filtered data frames for the observations and model data for the day of interest.  

# Next, we will select only the temperature data for comparing the model and observations. This needs to be done on a per station basis. Note that we use the function `drop_na()` to drop rows which contain NaN or NA data. Since each variable is measured at different intervals, not all variables will have data available at every time step in the ASOS data. The functions may be too sensitive to missing data and thus we take care to remvove it here from the observations, after we have isolated a particular variable.

## Temperature Time-Matching  
```

```{r Temperature Time Matching - JFK, echo = FALSE}
mystation = "JFK"
# Observational data used to probe WRF date-times.
asos.df_doi_jfk = asos.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% drop_na() %>% mutate(Source = "ASOS")

# WRF data for JFK for each forecast init time.
wrfd0.df_doi_jfk = wrfd0.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% mutate(Source = "WRF D-0")
wrfd1.df_doi_jfk = wrfd1.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% mutate(Source = "WRF D-1")
wrfd2.df_doi_jfk = wrfd2.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% mutate(Source = "WRF D-2")

# Perform Time-matching operation
wrfd0.df_doi_jfk_match = nearest.dataframe(asos.df_doi_jfk, wrfd0.df_doi_jfk, fcst_start0)
wrfd1.df_doi_jfk_match = nearest.dataframe(asos.df_doi_jfk, wrfd1.df_doi_jfk, fcst_start0)
wrfd2.df_doi_jfk_match = nearest.dataframe(asos.df_doi_jfk, wrfd2.df_doi_jfk, fcst_start0)

# Combine into a single data frame
jfk_temp = rbind(asos.df_doi_jfk, wrfd0.df_doi_jfk_match, wrfd1.df_doi_jfk_match, wrfd2.df_doi_jfk_match)
```


```{r Temperature Time Matching - LGA, echo = FALSE}
mystation = "LGA"
# Observational data used to probe WRF date-times.
asos.df_doi_lga = asos.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% drop_na() %>% mutate(Source = "ASOS")

# WRF data for LGA for each forecast init time.
wrfd0.df_doi_lga = wrfd0.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% mutate(Source = "WRF D-0")
wrfd1.df_doi_lga = wrfd1.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% mutate(Source = "WRF D-1")
wrfd2.df_doi_lga = wrfd2.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% mutate(Source = "WRF D-2")

# Perform Time-matching operation
wrfd0.df_doi_lga_match = nearest.dataframe(asos.df_doi_lga, wrfd0.df_doi_lga, fcst_start0)
wrfd1.df_doi_lga_match = nearest.dataframe(asos.df_doi_lga, wrfd1.df_doi_lga, fcst_start0)
wrfd2.df_doi_lga_match = nearest.dataframe(asos.df_doi_lga, wrfd2.df_doi_lga, fcst_start0)

# Combine into a single data frame
lga_temp = rbind(asos.df_doi_lga, wrfd0.df_doi_lga_match, wrfd1.df_doi_lga_match, wrfd2.df_doi_lga_match)
```

### Location: NYC

```{r Temperature Time Matching - NYC, echo = FALSE}
mystation = "NYC"
# Observational data used to probe WRF date-times.
asos.df_doi_nyc = asos.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% drop_na() %>% mutate(Source = "ASOS")

# WRF data for NYC for each forecast init time.
wrfd0.df_doi_nyc = wrfd0.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% mutate(Source = "WRF D-0")
wrfd1.df_doi_nyc = wrfd1.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% mutate(Source = "WRF D-1")
wrfd2.df_doi_nyc = wrfd2.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Temperature) %>% mutate(Source = "WRF D-2")

# Perform Time-matching operation
wrfd0.df_doi_nyc_match = nearest.dataframe(asos.df_doi_nyc, wrfd0.df_doi_nyc, fcst_start0)
wrfd1.df_doi_nyc_match = nearest.dataframe(asos.df_doi_nyc, wrfd1.df_doi_nyc, fcst_start0)
wrfd2.df_doi_nyc_match = nearest.dataframe(asos.df_doi_nyc, wrfd2.df_doi_nyc, fcst_start0)

# Combine into a single data frame
nyc_temp = rbind(asos.df_doi_nyc, wrfd0.df_doi_nyc_match, wrfd1.df_doi_nyc_match, wrfd2.df_doi_nyc_match)
```

```{r Wind Speed Time Matching - JFK, echo = FALSE}
## Wind Speed Time-Matching  

mystation = "JFK"
# Observational data used to probe WRF date-times.
asos.df_doi_jfk = asos.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% drop_na() %>% mutate(Source = "ASOS")

# WRF data for JFK for each forecast init time.
wrfd0.df_doi_jfk = wrfd0.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% mutate(Source = "WRF D-0")
wrfd1.df_doi_jfk = wrfd1.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% mutate(Source = "WRF D-1")
wrfd2.df_doi_jfk = wrfd2.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% mutate(Source = "WRF D-2")

# Perform Time-matching operation
wrfd0.df_doi_jfk_match = nearest.dataframe(asos.df_doi_jfk, wrfd0.df_doi_jfk, fcst_start0)
wrfd1.df_doi_jfk_match = nearest.dataframe(asos.df_doi_jfk, wrfd1.df_doi_jfk, fcst_start0)
wrfd2.df_doi_jfk_match = nearest.dataframe(asos.df_doi_jfk, wrfd2.df_doi_jfk, fcst_start0)

# Combine into a single data frame
jfk_wspd = rbind(asos.df_doi_jfk, wrfd0.df_doi_jfk_match, wrfd1.df_doi_jfk_match, wrfd2.df_doi_jfk_match)
```


```{r Wind Speed Time Matching - LGA, echo = FALSE}
### Location: LGA

mystation = "LGA"
# Observational data used to probe WRF date-times.
asos.df_doi_lga = asos.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% drop_na() %>% mutate(Source = "ASOS")

# WRF data for LGA for each forecast init time.
wrfd0.df_doi_lga = wrfd0.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% mutate(Source = "WRF D-0")
wrfd1.df_doi_lga = wrfd1.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% mutate(Source = "WRF D-1")
wrfd2.df_doi_lga = wrfd2.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% mutate(Source = "WRF D-2")

# Perform Time-matching operation
wrfd0.df_doi_lga_match = nearest.dataframe(asos.df_doi_lga, wrfd0.df_doi_lga, fcst_start0)
wrfd1.df_doi_lga_match = nearest.dataframe(asos.df_doi_lga, wrfd1.df_doi_lga, fcst_start0)
wrfd2.df_doi_lga_match = nearest.dataframe(asos.df_doi_lga, wrfd2.df_doi_lga, fcst_start0)

# Combine into a single data frame
lga_wspd = rbind(asos.df_doi_lga, wrfd0.df_doi_lga_match, wrfd1.df_doi_lga_match, wrfd2.df_doi_lga_match)
```

```{r Wind Speed Time Matching - NYC, echo = FALSE}
### Location: NYC

mystation = "NYC"
# Observational data used to probe WRF date-times.
asos.df_doi_nyc = asos.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% drop_na() %>% mutate(Source = "ASOS")

# WRF data for NYC for each forecast init time.
wrfd0.df_doi_nyc = wrfd0.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% mutate(Source = "WRF D-0")
wrfd1.df_doi_nyc = wrfd1.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% mutate(Source = "WRF D-1")
wrfd2.df_doi_nyc = wrfd2.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Speed) %>% mutate(Source = "WRF D-2")

# Perform Time-matching operation
wrfd0.df_doi_nyc_match = nearest.dataframe(asos.df_doi_nyc, wrfd0.df_doi_nyc, fcst_start0)
wrfd1.df_doi_nyc_match = nearest.dataframe(asos.df_doi_nyc, wrfd1.df_doi_nyc, fcst_start0)
wrfd2.df_doi_nyc_match = nearest.dataframe(asos.df_doi_nyc, wrfd2.df_doi_nyc, fcst_start0)

# Combine into a single data frame
nyc_wspd = rbind(asos.df_doi_nyc, wrfd0.df_doi_nyc_match, wrfd1.df_doi_nyc_match, wrfd2.df_doi_nyc_match)
```

```{r Wind Direction Time Matching - JFK, echo = FALSE}
## Wind Direction Time-Matching  

### Location: JFK

mystation = "JFK"
# Observational data used to probe WRF date-times.
asos.df_doi_jfk = asos.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% drop_na() %>% mutate(Source = "ASOS")

# WRF data for NYC for each forecast init time.
wrfd0.df_doi_jfk = wrfd0.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% mutate(Source = "WRF D-0")
wrfd1.df_doi_jfk = wrfd1.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% mutate(Source = "WRF D-1")
wrfd2.df_doi_jfk = wrfd2.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% mutate(Source = "WRF D-2")

# Perform Time-matching operation
wrfd0.df_doi_jfk_match = nearest.dataframe(asos.df_doi_jfk, wrfd0.df_doi_jfk, fcst_start0)
wrfd1.df_doi_jfk_match = nearest.dataframe(asos.df_doi_jfk, wrfd1.df_doi_jfk, fcst_start0)
wrfd2.df_doi_jfk_match = nearest.dataframe(asos.df_doi_jfk, wrfd2.df_doi_jfk, fcst_start0)

# Combine into a single data frame
jfk_wdir = rbind(asos.df_doi_jfk, wrfd0.df_doi_jfk_match, wrfd1.df_doi_jfk_match, wrfd2.df_doi_jfk_match)
```

```{r Wind Direction Time Matching - LGA, echo = FALSE}
### Location: LGA  

mystation = "LGA"
# Observational data used to probe WRF date-times.
asos.df_doi_lga = asos.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% drop_na() %>% mutate(Source = "ASOS")

# WRF data for NYC for each forecast init time.
wrfd0.df_doi_lga = wrfd0.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% mutate(Source = "WRF D-0")
wrfd1.df_doi_lga = wrfd1.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% mutate(Source = "WRF D-1")
wrfd2.df_doi_lga = wrfd2.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% mutate(Source = "WRF D-2")

# Perform Time-matching operation
wrfd0.df_doi_lga_match = nearest.dataframe(asos.df_doi_lga, wrfd0.df_doi_lga, fcst_start0)
wrfd1.df_doi_lga_match = nearest.dataframe(asos.df_doi_lga, wrfd1.df_doi_lga, fcst_start0)
wrfd2.df_doi_lga_match = nearest.dataframe(asos.df_doi_lga, wrfd2.df_doi_lga, fcst_start0)

# Combine into a single data frame
lga_wdir = rbind(asos.df_doi_lga, wrfd0.df_doi_lga_match, wrfd1.df_doi_lga_match, wrfd2.df_doi_lga_match)
```

```{r Wind Direction Time Matching - NYC, echo = FALSE}
### Location: NYC

mystation = "NYC"
# Observational data used to probe WRF date-times.
asos.df_doi_nyc = asos.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% drop_na() %>% mutate(Source = "ASOS")

# WRF data for NYC for each forecast init time.
wrfd0.df_doi_nyc = wrfd0.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% mutate(Source = "WRF D-0")
wrfd1.df_doi_nyc = wrfd1.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% mutate(Source = "WRF D-1")
wrfd2.df_doi_nyc = wrfd2.df_doi %>% filter(Station == mystation) %>% select(Date.Time, Wind.Direction) %>% mutate(Source = "WRF D-2")

# Perform Time-matching operation
wrfd0.df_doi_nyc_match = nearest.dataframe(asos.df_doi_nyc, wrfd0.df_doi_nyc, fcst_start0)
wrfd1.df_doi_nyc_match = nearest.dataframe(asos.df_doi_nyc, wrfd1.df_doi_nyc, fcst_start0)
wrfd2.df_doi_nyc_match = nearest.dataframe(asos.df_doi_nyc, wrfd2.df_doi_nyc, fcst_start0)

# Combine into a single data frame
nyc_wdir = rbind(asos.df_doi_nyc, wrfd0.df_doi_nyc_match, wrfd1.df_doi_nyc_match, wrfd2.df_doi_nyc_match)
```

# Forecast Hour Evaluation for JFK  

```{r JFK Eval Stats, echo = FALSE, echo = FALSE}
# For the temperature data I will use Bias, RMSE and MAE for the comparison statistics

eeval_temp = evalsys_temp(jfk_temp)
eeval_wspd = evalsys_wspd(jfk_wspd)
eeval_wdir = evalsys_wdir(jfk_wdir)
write.table(eeval_temp$eeval_df, "./jfk_temp_evaltimeseries.csv", row.names = FALSE, sep = ",")
write.table(eeval_wspd$eeval_df, "./jfk_wspd_evaltimeseries.csv", row.names = FALSE, sep = ",")
write.table(eeval_wdir$eeval_df, "./jfk_wdir_evaltimeseries.csv", row.names = FALSE, sep = ",")
write.table(eeval_temp$eeval_table, "./jfk_temp_evaltable.csv", row.names = FALSE, sep = ",")
write.table(eeval_wspd$eeval_table, "./jfk_wspd_evaltable.csv", row.names = FALSE, sep = ",")
write.table(eeval_wdir$eeval_table, "./jfk_wdir_evaltable.csv", row.names = FALSE, sep = ",")
```

```{r Temp Eval Table jfk, results='asis', echo = FALSE}
# Table of Daily Statistics
stargazer::stargazer(eeval_temp$eeval_table, type = 'html', summary = FALSE, column.sep.width = "100pt", 
                     title = "JFK - WRF 2-m Temperature (K) Performance")
```

```{r WSPD Eval Table jfk, results='asis', echo = FALSE}
# Table of Daily Statistics
stargazer::stargazer(eeval_wspd$eeval_table, type = 'html', summary = FALSE, column.sep.width = "100pt", title = "JFK - WRF 10-m Wind Speed (m/s) Performance")
```

```{r WDIR Eval Table jfk, results='asis', echo = FALSE}
# Table of Daily Statistics
stargazer::stargazer(eeval_wdir$eeval_table, type = 'html', summary = FALSE, column.sep.width = "100pt", 
                     title = "JFK - WRF 10-m Wind Direction (degN) Performance")
```

# Forecast Hour Evaluation for LGA

```{r LGA Eval Stats, echo = FALSE}
eeval_temp = evalsys_temp(lga_temp)
eeval_wspd = evalsys_wspd(lga_wspd)
eeval_wdir = evalsys_wdir(lga_wdir)
write.table(eeval_temp$eeval_df, "./lga_temp_evaltimeseries.csv", row.names = FALSE, sep = ",")
write.table(eeval_wspd$eeval_df, "./lga_wspd_evaltimeseries.csv", row.names = FALSE, sep = ",")
write.table(eeval_wdir$eeval_df, "./lga_wdir_evaltimeseries.csv", row.names = FALSE, sep = ",")
write.table(eeval_temp$eeval_table, "./lga_temp_evaltable.csv", row.names = FALSE, sep = ",")
write.table(eeval_wspd$eeval_table, "./lga_wspd_evaltable.csv", row.names = FALSE, sep = ",")
write.table(eeval_wdir$eeval_table, "./lga_wdir_evaltable.csv", row.names = FALSE, sep = ",")

```

```{r Temp Eval Table lga, results='asis', echo = FALSE}
# Table of Daily Statistics
stargazer::stargazer(eeval_temp$eeval_table, type = 'html', summary = FALSE, column.sep.width = "100pt", title = "LGA - WRF 2-m Temperature (K) Performance")
```

```{r WSPD Eval Table lga, results='asis', echo = FALSE}
# Table of Daily Statistics
stargazer::stargazer(eeval_wspd$eeval_table, type = 'html', summary = FALSE, column.sep.width = "100pt", title = "LGA - WRF 10-m Wind Speed (m/s) Performance")
```

```{r WDIR Eval Table lga, results='asis', echo = FALSE}
# Table of Daily Statistics
stargazer::stargazer(eeval_wdir$eeval_table, type = 'html', summary = FALSE, column.sep.width = "100pt", 
                     title = "LGA - WRF 10-m Wind Direction (degN) Performance")
```

# Forecast Hour Evaluation for NYC

```{r NYC Eval Stats, echo = FALSE}
eeval_temp = evalsys_temp(nyc_temp)
eeval_wspd = evalsys_wspd(nyc_wspd)
eeval_wdir = evalsys_wdir(nyc_wdir)
write.table(eeval_temp$eeval_df, "./nyc_temp_evaltimeseries.csv", row.names = FALSE, sep = ",")
write.table(eeval_wspd$eeval_df, "./nyc_wspd_evaltimeseries.csv", row.names = FALSE, sep = ",")
write.table(eeval_wdir$eeval_df, "./nyc_wdir_evaltimeseries.csv", row.names = FALSE, sep = ",")
write.table(eeval_temp$eeval_table, "./nyc_temp_evaltable.csv", row.names = FALSE, sep = ",")
write.table(eeval_wspd$eeval_table, "./nyc_wspd_evaltable.csv", row.names = FALSE, sep = ",")
write.table(eeval_wdir$eeval_table, "./nyc_wdir_evaltable.csv", row.names = FALSE, sep = ",")
```

```{r Temp Eval Table nyc, results='asis', echo = FALSE}
# Table of Daily Statistics
stargazer::stargazer(eeval_temp$eeval_table, type = 'html', summary = FALSE, column.sep.width = "100pt", title = "NYC - WRF 2-m Temperature (K) Performance")
```

```{r WSPD Eval Table nyc, results='asis', echo = FALSE}
# Table of Daily Statistics
stargazer::stargazer(eeval_wspd$eeval_table, type = 'html', summary = FALSE, column.sep.width = "100pt", title = "NYC - WRF 10-m Wind Speed (m/s) Performance")
```

```{r WDIR Eval Table nyc, results='asis', echo = FALSE}
# Table of Daily Statistics
stargazer::stargazer(eeval_wdir$eeval_table, type = 'html', summary = FALSE, column.sep.width = "100pt", title = "NYC - WRF 10-m Wind Direction (degN) Performance")
```

# uWRF and ASOS Time-Series Visualization
## Temperature  

```{r All_Locations_Temperature, echo = FALSE, fig.align = "center"}
# theme_set(theme_pubr())
jfk_temp$Location = "JFK"
lga_temp$Location = "LGA"
nyc_temp$Location = "NYC"
ggtemp = rbind(jfk_temp, lga_temp, nyc_temp)

ggplot(ggtemp, aes(x = Date.Time, y = Temperature, color = Source)) + 
       geom_point() + geom_line() + ylim(290, 320) +
       xlab("Time of Day") + ylab("Temperature (K)") + ggtitle("Forecast-Hour Evaluation: Temperature") + 
       facet_grid(Location ~ .) + theme(legend.position = "top") + scale_color_discrete(name = "")
```
  
## Wind Speed  

```{r All_Locations_Wind_Speed, echo = FALSE, fig.align = "center"}
jfk_wspd$Location = "JFK"
lga_wspd$Location = "LGA"
nyc_wspd$Location = "NYC"
ggwspd = rbind(jfk_wspd, lga_wspd, nyc_wspd)

ggplot(ggwspd, aes(x = Date.Time, y = Wind.Speed, color = Source)) + 
       geom_point() + geom_line() + ylim(0, 15) +
       xlab("Time of Day") + ylab("Wind Speed (m/s)") + ggtitle("Forecast-Hour Evaluation: Wind Speed (m/s)") + 
       facet_grid(Location ~ .) + theme(legend.position = "top") + scale_color_discrete(name = "")
```
  
## Wind Direction  

```{r All_Locations_Wind_Direction, echo = FALSE, fig.align = "center"}
jfk_wdir$Location = "JFK"
lga_wdir$Location = "LGA"
nyc_wdir$Location = "NYC"
ggwdir = rbind(jfk_wdir, lga_wdir, nyc_wdir)

ggplot(ggwdir, aes(x = Date.Time, y = Wind.Direction, color = Source)) + 
       geom_point() + geom_line() + ylim(0, 360) +
       xlab("Time of Day") + ylab("Wind Direction (degN)") + ggtitle("Forecast-Hour Evaluation: Wind Direction (degN)") + 
       facet_grid(Location ~ .) + theme(legend.position = "top") + scale_color_discrete(name = "")
```
