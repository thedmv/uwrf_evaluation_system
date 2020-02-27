## CCNY WindLidar -- ANNIE
# Functions

# Organize the data file from the csvs into a single data frame representing the UTC date of interest
# heights 


read.wlccny = function(mm, dd, yyyy) {
    # This function reads the csv files for the corresponding date
    # and grouping the data from each file into a single data frame
    # for that day.
    
    # Input: mm, dd, yyyy, are strings of the month (01-12), day (01-31),
    # and year (with century) of the date of interest.
    
    # NOTE: There are files that have no data and therefore yield
    # a file.size() of zero. The code uses that information to 
    # select only the files with data. Otherwise the read.table
    # does not work.
    # NOTE: If the day is filled with csv files with no data, this function
    # may not generate an error
    
    # Parse the selected date into the path name
    date.path = paste0("C:/Users/melev/Desktop/FORCE/WindLidar_CCNY/leosphere_annie/", mm, dd, "_reconstruction_wind_data", yyyy, "/", mm, "/", dd, "/")
    filelist = list.files(date.path, full.names = TRUE)
    
    # Valid files do not have a 0 KB file size, ie they contain data.
    validfilelist = filelist[file.size(filelist) > 0]
    groupedday = list()
    cnt = 1 # file counter
    for (ff in validfilelist) {
        groupedday[[cnt]] = read.table(ff, header = TRUE, stringsAsFactors = FALSE, sep = ";")
        cnt = cnt + 1
    }
    
    # Group all the data in the list into one data frame.
    dfwl = do.call(rbind, groupedday)
    
    return(dfwl)
}

# Fix the elevation and azimuth angles and consolidate columns

wlccny_clean_elevazim = function(dfwl) {
    # This function does two primary tasks:
    # 1. Creates a Elevation-Azimuth column that consolidates the information
    #    from the individual columns of Elevation and Azimuth. Keep only the
    #    resulting ElevAzim column.
    # 2. Eliminates the columns that deal with the scan and/or configuration.
    
    # Task 1
    # Round the elevation and azimuth data.
    dfwl$Azimuth.... = round(dfwl$Azimuth...., digits = 0)
    dfwl$Elevation.... = round(dfwl$Elevation...., digits = 0)
    dfwl$ElevAzim = paste0(dfwl$Elevation...., "-", dfwl$Azimuth....)
    
    # Change the 75-0 to 75-360
    dfwl$ElevAzim[dfwl$ElevAzim %in% c("75-0")] = "75-360" # Change the angles of 75-0 to 75-360
    
    # Change the 90-360 to 90-0; change 90-180 to 90-0
    dfwl$ElevAzim[dfwl$ElevAzim %in% c("90-360", "90-180")] = "90-0" # Change the angles of 75-0 to 75-360
    
    # Task 2
    # Select columns of interest.
    dfwl = dfwl %>% select(c(Timestamp, ElevAzim, Range..m., Horizontal.Wind.Speed..m.s.,
                         Vertical.Wind.Speed..m.s., Horizontal.Wind.Direction....,
                         CNR..dB., Confidence.Index.Status))
    
    return(dfwl)
}

xtimelimits <- function(xva, tz) {
    # This function computes the two element vector used to define the limits
    # on the date-time for data that has data over a 24-hour period.
    # xva - This a data frame with columns of Date.Time (a POSIXct class object),
    #       year, mon, day, hour, min, sec columns. The columns can be generated just
    #       Date.Time by using the function 'dateTimeCol'.
    # fc_xlim - This is the x axis limits.
    
    startdate = paste0(unique(xva$year), 
                       sprintf("%02d", unique(xva$mon)), #pad single digit with 0's
                       # sprintf("%02d", day),
                       sprintf("%02d", unique(xva$day)), # Need this for plotter_quiver_obs.R
                       "-", "00:00")
    enddate = paste0(unique(xva$year), 
                     sprintf("%02d", unique(xva$mon)),
                     # sprintf("%02d", day),
                     sprintf("%02d", unique(xva$day)), # Need this for plotter_quiver_obs.R
                     "-", "23:59")
    xtlim = c(as.POSIXct(startdate, format = "%Y%m%d-%H:%M", tz = tz),
              as.POSIXct(enddate,   format = "%Y%m%d-%H:%M", tz = tz))
    xtlim = c(head(xtlim, 1), tail(xtlim, 1)) # If time contains minutes from previous day, take the two extreme dates as the xlims
    attr(xtlim, "tzone") <- tz
    return(xtlim)
}

# Perform time average of the CCNY wind lidar data.

timeavg.ccnywindlidar <- function(wlvar, agg.time, statistic) {
    # This function perfomrs the time 'statistics' on the data as well as setting
    # up the rest of the columns. The difference from the other function is the addition
    # of the 'angpairs' column.
    
    # 1. For the timeAverage to work we need to change the Date.Time column name
    names(wlvar)[names(wlvar) == "Date.Time"] = "date"
    
    # 2. Save the heights in the variables
    timenames = c("date", "year", "mon", "day", "hour", "min", "sec")
    htsnames = names(wlvar)[names(wlvar) %w/o% timenames]
    
    # 4. Calculate the time average
    wlvar.aggregate = timeAverage(wlvar %>% select(date, htsnames), avg.time = agg.time, statistic =  statistic)
    names(wlvar.aggregate)[names(wlvar.aggregate) == "date"] = "Date.Time"
    wlvar.aggregate = wlvar.aggregate %>% dateTimeCol() %>% select(Date.Time, year:sec, htsnames)
    
    return(wlvar.aggregate)
} 

# Clean the wind speed, wind direction and vertical wind speed reconstruction data.
wlccny_fullclean_wvar <- function(mm, dd, yyyy, wvar) {
	# This function performs the cleaning operation on the CCNY Wind Lidar CSV'scan
	# of the reconstruction data.
	# Input:
	# mm, dd, yyyy - UTC date of interest. mm- month, dd - day, yyyy - year
	# wvar - String indicating reconstruction variable of interest: wspd, wdir, vesp
	
	# Select correct column from reconstruction data
	wvarname = switch(wvar,
	                  wspd = "Horizontal.Wind.Speed..m.s.",
	                  wdir = "Horizontal.Wind.Direction....",
	                  vesp = "Vertical.Wind.Speed..m.s.")
					  
	# Read, clean, and spread the data
	wvardf = read.wlccny(mm, dd, yyyy) %>%
             wlccny_clean_elevazim() %>%
             select(c(Timestamp:Range..m., wvarname)) %>%
             spread(Range..m., wvarname) 
			 
	# Filter depending if horizontal or vertical componenets selected.
	if ( (wvarname == "Horizontal.Wind.Speed..m.s.") | 
	     (wvarname == "Horizontal.Wind.Direction....") ) {
			wvardf = wvardf %>% filter(ElevAzim != "90-0")
		 } else	if(wvarname == "Vertical.Wind.Speed..m.s.") {
			wvardf = wvardf %>% filter(ElevAzim == "90-0")
		 } else {
			print("Something went horribly wrong in wlccny_clean_wvar!!")
		 }
		 
	# Change height names for the individual variable
	heightcols = names(wvardf) %w/o% c("Timestamp", "ElevAzim")
	names(wvardf)[heightcols] = paste0("X", names(wvardf)[heightcols])

	# Apply Time Zone to the individual variable
	wvardf$Timestamp = as.POSIXct(wvardf$Timestamp, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
	names(wvardf)[names(wvardf) == "Timestamp"] = "Date.Time"
	wvardf = wvardf %>% dateTimeCol() %>% select(Date.Time, year:sec, ElevAzim, X50:X2975)
	
	return(wvardf)
}

# TODO: Perform a clean on the CNR

# Perform a clean on the Quality Flag data
wlccny_fullclean_qflag <- function(mm, dd, yyyy, wvar) {
    # This function assembles a clean data frame of the Confidence Index Status.
    # Input:
    # mm, dd, yyyy - UTC date of interest. mm- month, dd - day, yyyy - year
    # wvar - String indicating which qflag of interest: wswd, vesp
    
    # Select correct column from reconstruction data
    wvarname = "Confidence.Index.Status"
    
    # Spread the data
    qflagdf = read.wlccny(mm, dd, yyyy) %>%
             wlccny_clean_elevazim() %>%
             select(c(Timestamp:Range..m., wvarname)) %>%
             spread(Range..m., wvarname)
    
    # Filter depending if horizontal or vertical
    if (wvar == "wswd") {
        qflagdf = qflagdf %>% filter(ElevAzim != "90-0")
    } else	if(wvar == "vesp") {
        qflagdf = qflagdf %>% filter(ElevAzim == "90-0")
    } else {
        print("Something went horribly wrong in wlccny_fullclean_qflag!!")
    }
    
    # Convert zeros to NaN
    qflagdf[qflagdf == 0] = NaN
    
    # Change height names
    heightcols = names(qflagdf) %w/o% c("Timestamp", "ElevAzim")
    names(qflagdf)[heightcols] = paste0("X", names(qflagdf)[heightcols])
    
    # Apply Time Zone
    qflagdf$Timestamp = as.POSIXct(qflagdf$Timestamp, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
    names(qflagdf)[names(qflagdf) == "Timestamp"] = "Date.Time"
    qflagdf = qflagdf %>% dateTimeCol() %>% select(Date.Time, year:sec, ElevAzim, X50:X2975)
    
    return(qflagdf)
}

wlccny_adjustheight <- function(clean_wvardf, baseheight = 65) {
    # This function takes a "clean" dataset adjusts the heights to reflect
    # heights from mean sea-level. The default `baseheight` is set to 65-m
    # to reflect the height of the CCNY Engineering building.
    # This can be used after the functions `wlccny_fullclean_wvar` and `wlccny_fullclean_qflag`.
    heightcols = names(clean_wvardf) %w/o% c("Date.Time", "year", "mon", "day", "hour", "min", "sec", "ElevAzim")
    names(clean_wvardf)[heightcols] = paste0("X", (char2height2(names(clean_wvardf)[heightcols]) + baseheight))
    return(clean_wvardf)
}

wlccny_adjustqc <- function(wvardf, qflagdf) {
    # This function applies the quality flag data.frame `qflagdf`
    # to the data.frame wvardf. Both data frames need to have gone through
    # the "fullclean" and "adjustheight" routines in order for the CCNY Wind
    # Lidar to be ready for analysis after this function.
    
    # Select the height columns
    heightcols = names(wvardf) %w/o% c("Date.Time", "year", "mon", "day", "hour", "min", "sec", "ElevAzim")
    
    
    # Quality control the data
    wvardf.qf = cbind(wvardf %>% select(Date.Time:ElevAzim), 
                      wvardf[heightcols] * qflagdf[heightcols], stringsAsFactors = FALSE)
    
    # Remove rows full of nans
    fullnanrowcnt = length(wvardf.qf[heightcols] %>% names())
    wvardf.qf$nansums = wvardf.qf[heightcols] %>% as.matrix() %>% is.nan() %>% rowSums()
    wvardf.qf = wvardf.qf %>% filter(nansums != fullnanrowcnt) %>% select(-nansums) # remove column of nansums
    
    return(wvardf.qf)
}


# Perform all the cleaning and adjusting procedures for the horizontal wind speed

wlccny_fullclean_fulladjust_wspd = function(mm, dd, yyyy) {
    # This function performs all the necessary procedures to load the horizontal
    # wind speed from the reconstruction data of the CCNY Wind lidar.
    
    qflag = wlccny_fullclean_qflag(mm, dd, yyyy, "wswd") %>% wlccny_adjustheight()
    wspd  = wlccny_fullclean_wvar(mm, dd, yyyy, "wspd") %>% wlccny_adjustheight() %>%
           wlccny_adjustqc(qflag)
    
    return(wspd)
}

wlccny_fullclean_fulladjust_wdir = function(mm, dd, yyyy) {
    # This function performs all the necessary procedures to load the horizontal
    # wind direction from the reconstruction data of the CCNY Wind lidar.
    
    qflag = wlccny_fullclean_qflag(mm, dd, yyyy, "wswd") %>% wlccny_adjustheight()
    wdir  = wlccny_fullclean_wvar(mm, dd, yyyy, "wdir") %>% wlccny_adjustheight() %>%
            wlccny_adjustqc(qflag)
    
    return(wdir)
}

wlccny_fullclean_fulladjust_vesp = function(mm, dd, yyyy) {
    # This function performs all the necessary procedures to load the vertical
    # wind speed from the reconstruction data of the CCNY Wind lidar.
    
    qflag = wlccny_fullclean_qflag(mm, dd, yyyy, "vesp") %>% wlccny_adjustheight()
    vesp  = wlccny_fullclean_wvar(mm, dd, yyyy, "vesp") %>% wlccny_adjustheight() %>%
            wlccny_adjustqc(qflag)
    
    return(vesp)
}
