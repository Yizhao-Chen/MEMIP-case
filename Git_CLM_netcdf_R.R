######################################################################################
################### LOOKING AT NETCDF FILES IN R ########################################
# R. Quinn Thomas (2016) and Kyla M. Dahlin (2014) CLM Tutorial
# this is written to be run line by line (or nearly so) and to think about the output
# NOT to be run as a script

# on yellowstone you can either run R by typing
# >module load R/3.2.2
# >R
# (but only do this for very simple tasks)
# or by opening an interactive job in geyser as we did with NCL
# >bsub -Ip -q geyser -W 2:00 -n 1 -P UCGD0002 xterm
# then in that interactive window type 
# >module load R/3.2.2
# >R

###################### FIRST INSTALL and LOAD PACKAGES #######################
install.packages(c("ncdf4", "raster", "rasterVis", "rgdal"))
# you'll be asked to select a CRAN mirror (I like CA1, personally) then you'll see lots of text go by
# these packages should be placed in a temporary directory. Also, you should only have to do this once
# on a given machine (yellowstone, your desktop, etc)

# then load the packages
library(ncdf4)
library(raster)
library(rasterVis)
library(rgdal)

############################# SET UP YOUR SPACE ##############################
# set a working directory (CHANGE THIS TO WHERE YOUR OUTPUT IS!!!!!)
setwd('/glade/scratch/dll/CLMTutorial2016_DataForAnalysis/I1850CLM50_001/')

############### READ IN AND LOOK AT A SINGLE .nc file ################

# open a single .nc file - you probably need to change this file name.
clm.tutorial.0001.01 <- nc_open('I1850CLM50_001.clm2.h0.0001-01.nc')
print(clm.tutorial.0001.01)  # look at all the data in one file

# extract a single variable
clm.0001.01.2m_air_temp <- ncvar_get(clm.tutorial.0001.01, varid = "TSA")

# turn this variable into a raster object
clm.0001.01.2m_air_temp.raster <- raster(clm.0001.01.2m_air_temp)

# create a "blank" version of the same raster (you'll use this later)
clm.in.raster.blank <- clm.0001.01.2m_air_temp.raster
clm.in.raster.blank[!is.na(clm.in.raster.blank)] <- NA

# open a plotting window
x11()

# plot the air temp raster
plot(clm.0001.01.2m_air_temp.raster)
# oh no! it's flipped on it's side!

# assign an extent (corners of the map)
extent(clm.0001.01.2m_air_temp.raster) <- c(-90, 90, 0, 360)

# flip the map counter-clockwise and rotate so it's split at 180 latitude instead of 0
clm.0001.01.2m_air_temp.raster <- rotate(flip(t(clm.0001.01.2m_air_temp.raster), direction = "y"))

# convert from Kelvin to degrees C
clm.0001.01.2m_air_temp.raster <- clm.0001.01.2m_air_temp.raster - 273.15

# plot it again to make sure it looks right
plot(clm.0001.01.2m_air_temp.raster)

######### NOW DO THE SAME THING IN A LOOP TO READ IN THE TEMP DATA FROM EACH MONTH / YEAR #########

# first make a vector of month nums
month.num <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

# create an empty output raster object
out.2m_air_temp <- clm.in.raster.blank

# loop through each of the files, reading them in and extracting the temp layer (each loop should take about 4 sec)
# note copy and paste the whole thing from "for" to the last "}"
for (y in 1:5) {
  for (m in 1:12) {
    # read in a file (CHANGE THIS IF YOUR OUTPUT FILES HAVE DIFFERENT NAMES!!!!)
    in.file <- nc_open(paste("I1850CLM50_001.clm2.h0.000", y, "-", month.num[m], ".nc", sep = ""))
    
    # get the temperature variable
    in.2m_air_temp <- ncvar_get(in.file, varid = "TSA")
    
    # turn that temperature variable into a raster
    in.2m_air_temp <- raster(in.2m_air_temp)
    
    # use the raster "brick" command to stack up the air temp files from each year/month
    out.2m_air_temp <- stack(out.2m_air_temp, in.2m_air_temp)
    
    # remove the giant .nc file just so it's not taking up memory
    rm(in.file)
    
    # make sure it's working... (each month should take a few seconds to read in)
    print(paste("done with year", y, "month", m))        
  }    
}

# as with the single month, we need to change the extent and rotate the data

# first remove the "blank" layer
out.2m_air_temp <- out.2m_air_temp[[-1]]

# assign an extent (corners of the map)
extent(out.2m_air_temp) <- c(-90, 90, 0, 360)

# flip the map counter-clockwise and rotate so it's split at 180 latitude instead of 0
out.2m_air_temp <- rotate(flip(t(out.2m_air_temp), direction = "y"))

# assign a projection (aka set the "coordinate reference system" or CRS)
crs(out.2m_air_temp) <- "+proj=longlat +ellps=WGS84" 

# convert from Kelvin to degrees C
out.2m_air_temp <- out.2m_air_temp - 273.15

# plot the first four layers to make sure they looks right
plot(out.2m_air_temp[[1:4]])

######################## CREATE ANNUAL TIME SERIES ###########################################

# calculate the mean for the first year
annual.mean.2m_air_temp <- mean(out.2m_air_temp[[1:12]], na.rm = T)

# calculate and stack the means for the 4 subsequent years
for (i in 1:4) {
  start <- i * 12
  end <- start + 11
  annual.mean.1 <- mean(out.2m_air_temp[[start:end]], na.rm = T)
  annual.mean.2m_air_temp <- stack(annual.mean.2m_air_temp, annual.mean.1)     
}

# look at all the years' means
plot(annual.mean.2m_air_temp)

# look at a difference map between two years, just to confirm that they're not the same
plot(annual.mean.2m_air_temp[[1]] - annual.mean.2m_air_temp[[2]])


############## CALCULATE ANNUAL MEAN AND INTERANNUAL STANDARD DEVIATION OF MEANS ##################

# get the number of years in your analysis (should be 5 here, but that may change for future stuff)
year.ct <- dim(annual.mean.2m_air_temp)[3] 

# just calculate the mean of the annual means (you could just calculate the mean of the whole
# stack, out.2m_air_temp, too)
mean.2m_air_temp <- mean(annual.mean.2m_air_temp, na.rm = T)

## two ways of calculating the standard deviation ##
# sadly this won't work: 
# stdev.2m_air_temp <- sd(annual.mean.2m_air_temp, na.rm = T)

# method 1: do the math yourself
stdev1.2m_air_temp <- sqrt(sum(((annual.mean.2m_air_temp - mean.2m_air_temp)^2), na.rm = T)/year.ct)

# method 2: use the raster package's function "calc" (more on this belo)
stdev2.2m_air_temp <- calc(annual.mean.2m_air_temp, fun = sd, na.rm = T)

# these two options don't actually give you identical numbers, but they're very similar!

#################################### MAKE DECENT LOOKING MAPS ######################################

# plot the means
x11() # open a new plotting window if you need to 
levelplot(mean.2m_air_temp, main = "Mean Annual Temperature", col.regions = topo.colors(20),  margin = FALSE)

# plot the standard deviations
x11() # open a second plotting window
levelplot(stdev2.2m_air_temp, main = "Interannual Standard Deviation of Means", col.regions = colorRampPalette(c("orange", "white", "blue"))(50), 
          margin = FALSE)

##################### CALCULATE SOMETHING MORE INTERESTING THAN MEAN AND SD #########################
# you can create any function and run it through "calc" to apply it to each gridcell
# here we build a function to identify which year had the hottest annual temp for each gridcell
# (if 2 years have the exact same annual mean, this will return the first one)

hottest.calc <- function(x, na.rm = T) {
  max.x <- max(x, na.rm = T)
  as.numeric(which(x == max.x)[1])
}

hottest.mean.year <- calc(annual.mean.2m_air_temp, fun = hottest.calc)

x11()
levelplot(hottest.mean.year, col.region = heat.colors(5), margin = F, main = "Hottest Avg Year of Series")

############################## NOW EXTRACT A REGION OF THE DATA #####################################
# we're going to look at Australia!

# define the extent of the area you're interested in
australia <- extent(c(110,155,-45,10))

# crop out the raster or raster stack
mean.2m_air_temp.AUS <- crop(mean.2m_air_temp, australia)

levelplot(mean.2m_air_temp.AUS, margin = F, main = "Mean Annual Temperature in Australia")


###################################### FINISHED? ######################################################

# to close R just type q() - R will then ask if you want to "save your workspace image" - I usually don't
# so type "n" (not in quotes)







