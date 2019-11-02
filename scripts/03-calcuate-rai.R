library(here)
library(tidyverse)
library(dplyr)
library(plyr)
library(progress)
library(sjmisc)

# read in cleaned record table (which already has dropped dates)
record.table <- read_csv(here::here('data', 'recordtable_year1and2_15min.csv'))
record.table$Date <- as.Date(record.table$Date)

# bring in master matrix for camera records - this is more info that we need, but it's used to extract the operation dates
# Side note: I recently learned that I should be using read_csv rather than read.csv
# read_csv actually keeps dates at top of columns rather than changing to weird format ("X.6.23.16") but all of the code below depends on that weird format, so we stick with it for now
matrix.all <- read.csv(here::here('data', 'All_species_by_date_year1and2_082919.csv'))

# put dates of interest into lists (one as character, one as list)
start.date.list <- list()
start.date.asdate.list <- list()
end.date.list <- list()
end.date.asdate.list <- list()
start.date.list[[1]] <- "X6.23.16" # no cameras up before this date
start.date.asdate.list[[1]] <- as.Date("6/23/16", format = "%m/%d/%y")
end.date.list[[1]] <- "X8.31.16"
end.date.asdate.list[[1]] <- as.Date("8/31/16", format = "%m/%d/%y")
start.date.list[[2]] <- "X6.1.17"
start.date.asdate.list[[2]] <- as.Date("6/1/17", format = "%m/%d/%y")
end.date.list[[2]] <- "X8.31.17"
end.date.asdate.list[[2]] <- as.Date("8/31/17", format = "%m/%d/%y")
start.date.list[[3]] <- "X6.1.18"
start.date.asdate.list[[3]] <- as.Date("6/1/18", format = "%m/%d/%y")
end.date.list[[3]] <- "X8.31.18"
end.date.asdate.list[[3]] <- as.Date("8/31/18", format = "%m/%d/%y")

# and make list of years for naming
year.list <- list()
year.list[[1]] <- 2016
year.list[[2]] <- 2017
year.list[[3]] <- 2018

# make a loop and calculate RAI for 2016, 2017, 2018

RAI.table.list <- list() # make list for storing output

for (i in 1:length(start.date.list)) {
  
  start.date <- start.date.list[[i]]
  end.date <- end.date.list[[i]]
  
  # create camera operation table
  camop <- dplyr::select(matrix.all, Camera, Species, start.date:end.date) # selects columns within specified dates
  
  zeroes <- row_count(camop, start.date:end.date, count = 0, append = FALSE) # count the columns with zeroes
  ones <- row_count(camop, start.date:end.date, count = 1, append = FALSE) # count the columns with ones
  Operation <- zeroes + ones # add zeroes and ones together to figure out how many days it operated for
  names(Operation) <- "Operation"
  
  # combine count dataframes with the matrix subset
  camop <- cbind(camop, Operation) 
  
  # change "Camera" name to "StudySite"
  names(camop) <- c("StudySite", names(camop)[2:ncol(camop)])
  
  # get just operation
  camop <- dplyr::select(camop, StudySite, Operation)
  
  # remove any cameras operating for <10 days (set Operation = 0)
  for (j in 1:nrow(camop)) {
    if (camop$Operation[j] < 10) {
      camop$Operation[j] <- 0
    }
  }
  
  # extract only the top 60 rows, since they are duplicated 42 times (one for each species)
  camop <- camop[1:60,] 

  ############# NOW RECORD TABLE
  
  start.date.asdate <- start.date.asdate.list[[i]]
  end.date.asdate <- end.date.asdate.list[[i]]
  
  # subset record table to date of interest
  record.table.subset <- record.table[record.table$Date >= start.date.asdate & record.table$Date <= end.date.asdate,]
  
  # just take the columns we need (camera and species)
  record.table.subset <- record.table.subset[,2:3]
  
  # change column names to match other datasets
  colnames(record.table.subset) = c("StudySite", "CommName")
  
  # remove species that we don't care about (have to remove Eland because it was only present in 1 of the years and is therefore messing things up far down the line)
  for (k in c("Bat", "Bird_other", "Duiker_unknown", "Eland", "Fire", "Ghost", "Ghosts Part 1", 
              "Ghosts Part 2", "Ground_hornbill", "Guineafowl_crested", 
              "Guineafowl_helmeted", "Hornbill_ground", "Hornbill_ground 2", 
              "Human", "Insect", "Mongoose_other", "Mongoose_unknown", 
              "Monitor_lizard", "Rain", "Reptile", "Rodent", "Setup",
              "Snake", "Unknown", "Unknown_antelope")) {
    record.table.subset <- record.table.subset[record.table.subset$CommName != k, ]
  }
  
  # calculates number of observations of each species at each camera
  records <- record.table.subset %>%
    dplyr::group_by(CommName, StudySite) %>%
    dplyr::summarise(Detections = n()) %>%     # counts number of observations of each species
    spread(key = CommName, value = Detections)  # gets from long to wide format
  
  # replace NA values with 0 (no detections)
  records[is.na(records)] <- 0
  
  # add together major species
  records$MajorSpecies <- NA
  for (l in 1:nrow(records)) {
    records$MajorSpecies[l] <- sum(c(records$Bushbuck[l], records$Impala[l], records$Nyala[l], 
                                     records$Oribi[l], records$Warthog[l], records$Waterbuck[l]), na.rm=TRUE)
  }

  # join camera operation dates and species observations
  RAI.table <- join(camop, records)
  
  # gather data so each species-camera is its own row again
  RAI.table <- RAI.table %>% gather(3:ncol(RAI.table), key = "CommName", value = "Count")
  
  # replace 0 with NA (not sure why they got un-replaced...)
  RAI.table[is.na(RAI.table)] <- 0
  
  # calculate RAI
  RAI.table$RAI <- RAI.table$Count / RAI.table$Operation  
  

  
  # add column for year
  RAI.table$Year <- year.list[[i]]
  
  # replace NaN with NA
  RAI.table[is.na(RAI.table)] <- NA
  
  # store table in list
  RAI.table.list[[i]] <- RAI.table
  
}

# now join them all together
RAI.table.master <- rbind(RAI.table.list[[1]], RAI.table.list[[2]], rbind(RAI.table.list[[3]]))

# now calculate overall RAI across all three years
RAI.table.grouped <- RAI.table.master %>%
  dplyr::group_by(StudySite, CommName) %>%
  dplyr::summarise(Operation = sum(Operation),
                   Count = sum(Count))
for (o in 1:nrow(RAI.table.grouped)){
  RAI.table.grouped$RAI[o] <- RAI.table.grouped$Count[o] / RAI.table.grouped$Operation[o]
}
RAI.table.grouped$RAI <- NA
RAI.table.grouped$RAI <- RAI.table.grouped$Count / RAI.table.grouped$Operation
RAI.table.grouped$Year <- "all_years"

# combine back in with master
RAI.table.master <- rbind.fill(RAI.table.grouped, RAI.table.master)

# export csv
write.csv(RAI.table.master, here::here('data', 'RAI_LOD_by_year_long.csv'), row.names=F)

# make it wide
RAI.table.master.wide <- RAI.table.master
RAI.table.master.wide$CommName_Year <- paste(RAI.table.master.wide$CommName, RAI.table.master.wide$Year, sep = "_")
RAI.table.master.wide <- RAI.table.master.wide %>% dplyr::select(StudySite, CommName_Year, RAI) %>% spread(key = CommName_Year, value = RAI)


# export csv
write.csv(RAI.table.master.wide, "Matt_LOD/RAI_LOD_by_year_wide_091919.csv", row.names=F)
write.csv(RAI.table.master.wide, here::here('data', 'RAI_LOD_by_year_wide.csv'), row.names=F)





# old code - back when I was doing it a stupider way

# change names (Operation, Count, RAI) to reflect year
#for (m in c(2,4:5)) {
#  names(RAI.table)[m] <- paste(names(RAI.table[m]), year.list[[i]], sep = "_")
#}

# now join them all together
#RAI.table.master <- full_join(RAI.table.list[[1]], RAI.table.list[[2]], by = c("StudySite", "CommName")) %>% full_join(RAI.table.list[[3]], by = c("StudySite", "CommName"))

# make new columns for all three years combined
#for (n in 1:nrow(RAI.table.master)) {
#  RAI.table.master$Operation_all_years[n] <- sum(RAI.table.master$Operation_2016[n], RAI.table.master$Operation_2017[n], RAI.table.master$Operation_2018[n], na.rm = TRUE)
#  RAI.table.master$Count_all_years[n] <- sum(RAI.table.master$Count_2016[n], RAI.table.master$Count_2017[n], RAI.table.master$Count_2018[n], na.rm = TRUE)
#  RAI.table.master$RAI_allyears[n] <- RAI.table.master$Count_all_years[n] / RAI.table.master$Operation_all_years[n]
#}


