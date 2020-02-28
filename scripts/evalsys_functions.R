# This script contains the functions used for the Evaluation System forecast.

#### Functions for TS Lists ####
read.wrf_tslistTS = function(tsfilename, fcst_hr0) {
    # This function reads TS files created from the tslist runtime configurations. 
    # Only works for the files in which the extension *.TS. First seen in `Load01-WRF_tslist.Rmd.`
    # INPUT
    # tsfilename  --  string of the file name; 
    # fcst_hr0    --  string of the first forecast hour; for example: "2018-06-30 18:00:00"
    # OUTPUT
    # tsdata  -- data frame of the tslist; only columns of date and time, temperature, and winds (u, v, wspd, wdir)
    
    # Read the WRF tslist file 
    tsdata = read.table(tsfilename, skip = 1, stringsAsFactors = FALSE)
    
    # Add the column names for the  Time series (TS) output of surface variables
    names(tsdata) = c("id", "ts_hour", "id_tsloc", "ix", "iy", "t", "q", "u", "v", "psfc", "glw", "gsw", "hfx", "lh", "tsk", "tslb_1", "rainc", "rainnc", "clw")
    
    # Set the start of the forecast hour
    tsdata$Date.Time = as.POSIXct(fcst_hr0, tz = "UTC") + tsdata$ts_hour * 3600
    
    # Here we add columns for the date and times so that the wspd/wdir finctuon can work
    tsdata = tsdata %>% dateTimeCol() %>% select(Date.Time, year:sec, t, q, u, v)
    
    return(tsdata)
}

ts_wswd = function(tsdata) {
    # This function returns the wind speed and wind directions given u and v components.
    # INPUT
    # tsdata -- data.frame; read into memory by read.wrf_tslistTS
    # OUTPUT
    # tsdata -- data.frame with wind speed (wspd) and wind direction (wdir) columns added.
    
    # Set the u and v wind vector components
    uwind = tsdata$u
    vwind = tsdata$v
    
    # Calculate the wind speed
    tsdata$wspd = sqrt(uwind^2 + vwind^2)
    
    # Calculate the wind direction (https://stackoverflow.com/questions/1311049/how-to-map-atan2-to-degrees-0-360/25398191)
    tsdata$wdir = ((atan2(-uwind, -vwind) * 180/pi) + 360) %% 360
    
    # Return the new data frame with truncated columns from the tslist AND calculated wind speed and wind direction
    return(tsdata)
}

timeavg.ts = function(tsdata, agg.time = "1 hour", statistic = "mean") {
    # This function reads the "tsdata" object read into memory by `read.wrf_tslistTS` and 
    # calculates a statistic based on the desired time aggregate.
    # INPUT
    # tsdata  --  data.frame; read *.TS files using read.wrf_tslistTS
    # agg.time -- time used to form the aggregates that the statistics will be calculated for
    # statistic -- statistic operation to be used on the aggregates picked by "agg.time"
    #
    # NOTE: The `agg.time` and `statistic` inputs are used in the timeAverage function from the 
    #       openair library.
    #
    # OUTPUT
    # tsdata -- data.frame; averaged (or other statistic) for the agg.time chosen
    
    # 1. For timeAverage to work we need to change the Date.Time column name
    names(tsdata)[names(tsdata) == "Date.Time"] = "date"
    
    # 2. Save the heights in the variables
    timenames = c("date", "year", "mon", "day", "hour", "min", "sec")
    varnames = names(tsdata)[names(tsdata) %w/o% timenames]
    
    # 3. Calculate the time average
    tsdata = timeAverage(tsdata %>% select(date, varnames), avg.time = agg.time, statistic =  statistic)
    names(tsdata)[names(tsdata) == "date"] = "Date.Time"
    tsdata = tsdata %>% dateTimeCol() %>% select(Date.Time, year:sec, varnames)
    
    return(tsdata)
}


#### Functions for ASOS Stations ####
read.asos = function(fASOS) {
    # This function reads the ASOS data downloaded using the dl_ny_asos.py script.
    # In this function I am selecting only the columns that correspond to the temperature,
    # relative humidity, wind direction, and wind speed, in addition to the station name and
    # the time stamp.
    #
    # INPUT
    # fASOS -- string; file name of the ASOS data
    # OUTPUT
    # asos.df -- data.frame; Data frame of the ASOS data.
    #
    
    # Read ASOS data
    asos.df = read.table(fASOS, header = T, stringsAsFactors = FALSE, skip = 5, sep = ",")
    
    # Pick columns of interest
    select_columns = c("station", "valid", "tmpf", "relh", "drct", "sknt")
    asos.df = asos.df %>% select(select_columns)
    
    # Replace M with NaN
    asos.df[asos.df == "M"] = NaN
    
    # Rename columns to match WRF tslist columns
    rename_columns = c("station", "Date.Time", "t", "rh", "wdir", "wspd")
    names(asos.df) = rename_columns
    
    # Change variable types for every column
    asos.df$Date.Time = as.POSIXct(asos.df$Date.Time, tz = "UTC") # Change the string date.time to date-time objects
    asos.df$t = as.numeric(asos.df$t)                             # Temperature variables (Fahrenheit)
    asos.df$rh = as.numeric(asos.df$rh)                           # Relative Humidity (%)
    asos.df$wdir = as.numeric(asos.df$wdir)                       # Wind Direction (degrees N)
    asos.df$wspd = as.numeric(asos.df$wspd)                       # Wind Speed (in knots)
    
    # Return the new data frame
    return(asos.df)
}

#### Functions for Time-Matching ####

next.day = function(select.date) {
    # This function calculates the next date and outputs it as a string
    #
    # INPUT
    # select.date - string of the desired date in "YYY-mm-dd" format
    # OUTPUT
    # next.date - string representing the next day from select.date
    
    next.date = as.character(str_split(select.date, "-")[[1]] %>% as.numeric() + c(0, 0, 1)) %>% 
        str_pad(width = 2, side = "left", pad = 0) %>% str_flatten("-")
    return(next.date)
}

prev.day = function(select.date) {
    # This function calculates the previous date and outputs it as a string
    #
    # INPUT
    # select.date - string of the desired date in "YYY-mm-dd" format
    # OUTPUT
    # next.date - string representing the next day from select.date
    prev.date = as.character(str_split(select.date, "-")[[1]] %>% as.numeric() - c(0, 0, 1)) %>% 
        str_pad(width = 2, side = "left", pad = 0) %>% str_flatten("-")
    return(next.date)
}

nearest <- function(probe, target, ends=c(-Inf,Inf)) {
    # Return an array `i` of indexes into `target`, parallel to array `probe`.
    # For each index `j` in `target`, probe[i[j]] is nearest to target[j].
    # Both `probe` and `target` must be vectors of numbers in ascending order.
    #
    glb <- function(u, v) {
        n <- length(v)
        z <- c(v, u)
        j <- i <- order(z)
        j[j > n] <- -1
        k <- cummax(j)
        return (k[i > n])
    }
    y <- c(ends[1], target, ends[2])
    
    i.lower <- glb(probe, y)
    i.upper <- length(y) + 1 - rev(glb(rev(-probe), rev(-y)))
    y.lower <- y[i.lower]
    y.upper <- y[i.upper]
    lower.nearest <- probe - y.lower < y.upper - probe
    i <- ifelse(lower.nearest, i.lower, i.upper) - 1
    i[i < 1 | i > length(target)] <- NA
    return (i)
}

nearest.dataframe = function(probe, target, select.date) {
    # This function uses the function 'nearest' to return one data frame of the target whose time-series is a closest match
    # to the probe time-series. The date-times columns must be named "Date.Time" or else this function will not work.
    # In the first implementation of this code we use ASOS date times to probe the WRF TS list target.
    
    ############
    # INTPUT
    # probe - data frame of the observatiopn data
    # target - data frame of the target data
    # select.date - string in YYY-mm-dd format of the date of interest
    # OUTPUT
    # near.df - subset of target data frame that has matching time-series to the probe
    #
    
    # Select a single day to filter the data
    daybegin = select.date %>% as.POSIXct(., tz = "UTC")
    dayend   = daybegin + days(1)
    
    # Create the date-time arrays that will be used to match the times. 
    myprobe  = probe %>% filter(Date.Time >= daybegin) %>% filter(Date.Time <= dayend)
    mytarget = target %>% filter(Date.Time >= daybegin) %>% filter(Date.Time <= dayend)
    
    # Convert to numeric, from Date.Time object, so function works.
    myprobets = as.numeric(myprobe$Date.Time)
    mytargetts = as.numeric(mytarget$Date.Time)
    # Find the WRF (target) values that match with the observations (probe) of ASOS stations
    ii = nearest(myprobets, mytargetts)
    new_target = mytarget[ii, ]
    
    # Return the new target data frame
    return(new_target)
}

evalsys_temp = function(data.for.eval) {
    # This function calculates the evaluation statistics for temperature.
    # The output is a data frame of the data before the evaluation statistics,
    # As of 2020-02-26 the performance statistics used are Bias, RMSE and MAE.
    # INPUT
    # data.for.eval - data frame with three columns: Date.Time (UTC), Temperature, and
    #                 Source (ASoS, WRF D-0, WRF D-1, WRF D-2)
    # OUTPUT
    # eeval_table  -  list containing a table of the performance statistics for temperature,
    #                 and a data frame with the difference columns.z
    eeval = data.frame(Date.Time = (data.for.eval %>% filter(Source == "ASOS"))$Date.Time, 
                       OBS = (data.for.eval %>% filter(Source == "ASOS"))$Temperature,
                       WRF_D0 =  (data.for.eval %>% filter(Source == "WRF D-0"))$Temperature,
                       WRF_D1 =  (data.for.eval %>% filter(Source == "WRF D-1"))$Temperature,
                       WRF_D2 =  (data.for.eval %>% filter(Source == "WRF D-2"))$Temperature)
    eeval$DIFF_D0 = eeval$WRF_D0 - eeval$OBS
    eeval$DIFF_D1 = eeval$WRF_D1 - eeval$OBS
    eeval$DIFF_D2 = eeval$WRF_D2 - eeval$OBS
    d = data.frame(D0 = eeval$DIFF_D0, D1 = eeval$DIFF_D1, D2 = eeval$DIFF_D2)
    eeval_table = data.frame(Forecast.Init = c("WRF D-0", "WRF D-1", "WRF D-2"),
                             BIAS = c(sum(d$D0)/length(d$D0), 
                                      sum(d$D1)/length(d$D1), 
                                      sum(d$D2)/length(d$D2)),
                             RMSE = c(sqrt(sum(d$D0^2)/length(d$D0)), 
                                      sqrt(sum(d$D1^2)/length(d$D1)), 
                                      sqrt(sum(d$D2^2)/length(d$D2))),
                             MAE  = c(sum(abs(d$D0))/length(d$D0),
                                      sum(abs(d$D1))/length(d$D1),
                                      sum(abs(d$D2))/length(d$D2)) )
    eeval = list(eeval_df = eeval, eeval_table = eeval_table)
    return(eeval)
}

evalsys_wspd = function(data.for.eval) {
    # This function calculates the evaluation statistics for wind speed.
    # The output is a data frame of the data before the evaluation statistics,
    # As of 2020-02-26 the performance statistics used are Bias, RMSE and MAE.
    # INPUT
    # data.for.eval - data frame with three columns: Date.Time (UTC), Temperature, and
    #                 Source (ASoS, WRF D-0, WRF D-1, WRF D-2)
    # OUTPUT
    # eeval_table  -  list containing a table of the performance statistics for wind speed,
    #                 and a data frame with the difference columns.
    eeval = data.frame(Date.Time = (data.for.eval %>% filter(Source == "ASOS"))$Date.Time, 
                       OBS = (data.for.eval %>% filter(Source == "ASOS"))$Wind.Speed,
                       WRF_D0 =  (data.for.eval %>% filter(Source == "WRF D-0"))$Wind.Speed,
                       WRF_D1 =  (data.for.eval %>% filter(Source == "WRF D-1"))$Wind.Speed,
                       WRF_D2 =  (data.for.eval %>% filter(Source == "WRF D-2"))$Wind.Speed)
    eeval$DIFF_D0 = eeval$WRF_D0 - eeval$OBS
    eeval$DIFF_D1 = eeval$WRF_D1 - eeval$OBS
    eeval$DIFF_D2 = eeval$WRF_D2 - eeval$OBS
    d = data.frame(D0 = eeval$DIFF_D0, D1 = eeval$DIFF_D1, D2 = eeval$DIFF_D2)
    eeval_table = data.frame(Forecast.Init = c("WRF D-0", "WRF D-1", "WRF D-2"),
                             BIAS = c(sum(d$D0)/length(d$D0), 
                                      sum(d$D1)/length(d$D1), 
                                      sum(d$D2)/length(d$D2)),
                             RMSE = c(sqrt(sum(d$D0^2)/length(d$D0)), 
                                      sqrt(sum(d$D1^2)/length(d$D1)), 
                                      sqrt(sum(d$D2^2)/length(d$D2))),
                             MAE  = c(sum(abs(d$D0))/length(d$D0),
                                      sum(abs(d$D1))/length(d$D1),
                                      sum(abs(d$D2))/length(d$D2)) )
    eeval = list(eeval_df = eeval, eeval_table = eeval_table)
    return(eeval)
}

evalsys_wdir = function(data.for.eval) {
    # This function calculates the evaluation statistics for wind direction.
    # The output is a data frame of the data before the evaluation statistics,
    # As of 2020-02-26 the performance statistics used are Bias, RMSE and MAE.
    # INPUT
    # data.for.eval - data frame with three columns: Date.Time (UTC), Temperature, and
    #                 Source (ASoS, WRF D-0, WRF D-1, WRF D-2)
    # OUTPUT
    # eeval_table  -  list containing a table of the performance statistics for wind speed,
    #                 and a data frame with the difference columns.
    eeval = data.frame(Date.Time = (data.for.eval %>% filter(Source == "ASOS"))$Date.Time, 
                       OBS = (data.for.eval %>% filter(Source == "ASOS"))$Wind.Direction,
                       WRF_D0 =  (data.for.eval %>% filter(Source == "WRF D-0"))$Wind.Direction,
                       WRF_D1 =  (data.for.eval %>% filter(Source == "WRF D-1"))$Wind.Direction,
                       WRF_D2 =  (data.for.eval %>% filter(Source == "WRF D-2"))$Wind.Direction)
    
    # Calculate the difference between simulated and observed.
    eeval$DIFF_D0= eeval$WRF_D0 - eeval$OBS
    eeval$DIFF_D0[eeval$DIFF_D0 > 180]  = eeval$DIFF_D0[eeval$DIFF_D0 > 180]  - 360
    eeval$DIFF_D0[eeval$DIFF_D0 < -180] = eeval$DIFF_D0[eeval$DIFF_D0 < -180] + 360
    
    eeval$DIFF_D1= eeval$WRF_D1 - eeval$OBS
    eeval$DIFF_D1[eeval$DIFF_D1 > 180]  = eeval$DIFF_D1[eeval$DIFF_D1 > 180]  - 360
    eeval$DIFF_D1[eeval$DIFF_D1 < -180] = eeval$DIFF_D1[eeval$DIFF_D1 < -180] + 360
    
    eeval$DIFF_D2= eeval$WRF_D2 - eeval$OBS
    eeval$DIFF_D2[eeval$DIFF_D2 > 180]  = eeval$DIFF_D2[eeval$DIFF_D2 > 180]  - 360
    eeval$DIFF_D2[eeval$DIFF_D2 < -180] = eeval$DIFF_D2[eeval$DIFF_D2 < -180] + 360
    
    d = data.frame(D0 = eeval$DIFF_D0, D1 = eeval$DIFF_D1, D2 = eeval$DIFF_D2)
    eeval_table = data.frame(Forecast.Init = c("WRF D-0", "WRF D-1", "WRF D-2"),
                             RMSE = c(sqrt(sum(d$D0^2)/length(d$D0)), 
                                      sqrt(sum(d$D1^2)/length(d$D1)), 
                                      sqrt(sum(d$D2^2)/length(d$D2))),
                             MAE  = c(sum(abs(d$D0))/length(d$D0),
                                      sum(abs(d$D1))/length(d$D1),
                                      sum(abs(d$D2))/length(d$D2)) )
    eeval = list(eeval_df = eeval, eeval_table = eeval_table)
    return(eeval)
}
