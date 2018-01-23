

# necessary packages
library(sp) 	# provides classes for spatial data in R
library(raster)	# provides classes and methods for raster datasets

# necessary - but attention for Mac users!
library(rgdal)	# interface to the Geospatial Data Abstraction Library to read and write different spatial file formats

# usefull packages
library(maptools)	# some specific spatial methods and conversion methods for spatial data 
library(fields)		# contains methods for spatial interpolation and statistics but also nice color palettes fro mapping
library(RColorBrewer) # color palettes for mapping
library(marxan)
library(rasterVis)
# library(viridis)
#library(gstat)
library(marxan)
library(stringr)

library(ggplot2)
library(gridExtra)

crs.geo <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84")    # geographical, datum WGS84
sps<-list.files("data")
sp.list<-list()
e<- extent(-77.5, -71.1, 0, 11.5)
# for(i in 1:length(sps)){
#   sp_x<-raster(paste("data/", sps[i], "/", sps[i], ".tif", sep = ""))
#   proj4string(sp_x) <- crs.geo     # define projection system of our data
#   #extent (sp_x) <- e #put max extent
#   sp.list[[i]] <- sp_x
#   print(res(sp.list[[i]]))
#   print(extent(sp_x))
# }

raster_data <- list.files("data")
# w <- raster("data/Allobates_juanii/Allobates_juanii.tif", crs="+proj=longlat +datum=WGS84 +ellps=WGS84")
w <- raster("shp/boundbox1.tif", crs="+proj=longlat +datum=WGS84 +ellps=WGS84")


sp.list<-list()


# all values >= -1 and <= 0.1 become 0
m <- c(-1, 0.1, 0, 0.2, 2, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
# rc <- reclassify(w, rclmat)

fun <- function(x) { x[is.na(x)] <- 0; return(x)} # funtion to convert NA to O


for (i in 1:length(raster_data)){
  r <- raster(paste("data/", raster_data[i], "/", raster_data[i], ".tif", sep = ""), crs=crs.geo)
  rp <- projectRaster(from = r, to = w,
                      filename = file.path ("crop/", raster_data[i]),
                      method = "bilinear",
                      format = "raster",
                      overwrite = TRUE)
  rc <- reclassify(rp, rclmat) # corrige ceros
  rm <- merge(rc, w) # extiende todos los chiquitos a w
  rc3 <- calc(rm, fun) # remueve NAs
  sp.list[[i]] <- rc3 # adiciona a la lista
}


sp_staked<- stack(sp.list)

levelplot(sp_staked) # slow

# load marxan R package
library(marxan)

# load example data
data(taspu, tasinvis)

# real data
cons_unit<-readOGR("C:/Users/Diego/Documents/CodigoR/Aldemar_anuros/shp",
                layer="hexagons_dia01_clip", stringsAsFactors=FALSE)
cons_unit$id<-cons_unit$id+1 #start in 1
cons_unit$cost<-as.numeric(cons_unit$cost) #make cost numeric
cons_unit$id<-as.integer(cons_unit$id)
cons_unit$status<-as.integer(cons_unit$status)

# make reserve systems 
# lets asume first equal cost and no protected areas
# copy taspu
cons_unit_scn1<-cons_unit

# set costs
cons_unit_scn1@data$cost<-1

# set status values
# note the 'L' after the zero is used to indicate
# that we mean the integer zero and not the decimal
# place number zero
cons_unit_scn1@data$status<-0L

# Here, we will generate a portfolio of reserve systems that represent 25% of each species. 
# results<-marxan(taspu, tasinvis, targets="20%", NUMREPS=100L, BLM=0)
result_scn1<-marxan(cons_unit_scn1, sp_staked, 
                     targets="25%", # level of protection
                     spf=1, # species penalty factor for species
                     NUMREPS=100L, # controls the number of solutions in our portfolio
                     NCORES=2L, # number of threads for parallel processing
                     BLM=0)

# histogram of proportion of species adequately
# represented in each solution
# if many of the solutions adequately represented the classes
# most of the bins would be close to 1, whereas if 
# the solutions failed to represent the classes most of the
# bins would be close to zero
hist(rowMeans(targetsmet(result_scn1)), freq=TRUE, xlim=c(0,1), las=1,
     main='Histogram of representation in portfolio',
     ylab='Frequency of solutions',
     xlab='Proportion of frog species adequately represented'
    )

# geoplot for best solution
plot(result_scn1, 0)

# geoplot for selection frequencies
# plot(results)

# plot distribution of sp 5
spplot(result_scn1, 5, var='occ')

# plot richness in planning units
spplot(result_scn1, var='occ')

# geoplot richness in planning units
# with a satellite base map
spplot(result_scn1, var='occ', basemap='satellite')



# if all of the solutions in our portfolio have failed to meet the targets for the all species. 
# To fix this, we need to increase the species penalty factors (SPFs), and rerun MARXAN to generate a new portfolio of solutions.

# copy the MARXAN parameters and pre-processed data in results,
# update the SPF parameter for all species,
# run MARXAN,
# load the solutions back into R,
# store the solutions in results2
results2<-update(result_scn1, ~spp(1:42, spf=rep(100,42)))

# Lets compare two portfolios
# get levels of representation in each portfolio
results.repr<-rowMeans(targetsmet(result_scn1))
results2.repr<-rowMeans(targetsmet(results2))

# create 2 plotting areas in the one window
par(mfrow=c(1,2))

  # histogram of first portfolio
  hist(results.repr, freq=TRUE, xlim=c(0,1), las=1,
       ylab='Frequency of solutions',
       xlab='Proportion of species adequately represented',
       main="Level of representation with SPF=1"
  )
  
  # print best level of representation
  print(max(results.repr))
  
  # histogram of second portfolio
  # if you see a giant single rectangle this means
  # all the solutions have the same level of representation
  hist(results2.repr, freq=TRUE, xlim=c(0,1), las=1,
       ylab='Frequency of solutions',
       xlab='Proportion of species adequately represented',
       main="Level of representation with SPF=100"
  )
  
  # print best level of representation
  print(max(results2.repr))

  # make a geoplot of the best solution
  plot(result_scn1, 0)
  plot(results2, 0)
  
  # planing uniit selection frecuency
  plot(result_scn1, colramp='YlGnBu')
  plot(results2, colramp='YlGnBu')  
  
# The solutions in this portfolio are fairly fragmented. 
# All the solutions so far were made under the assumption that all planning units have equal acquisition costs 
# and that Col does not have any protected areas. Not true!
  
# get planning unit ids
pu.ids<-cons_unit@data$id
  
# get planning unit costs from ArcGis Table
costtable<-read.csv("C:/Users/Diego/Documents/CodigoR/Aldemar_anuros/shp/tablepop4.txt")
pu.costs<-costtable$MEAN # this is the column
  
# get planning unit statuses from original file
pu.status<-cons_unit@data$status
  
# copy input parameters and data in results2, 
# change planning unit costs and statuses
# rerun MARXAN,
# load outputs into R and store them in results3
results3<-update(results2, ~pu(pu.ids, cost=pu.costs, status=pu.status))

# Now we have a third portfolio. Let's compare the previous portfolio based on unrealistic planning unit data 
# with the new portfolio based on realistic planning unit data.
# geoplot showing differences between the best solution in each portfolio
plot(results2, results3, i=0, j=0)

# geoplot showing differences between the third solution
# in results2 and the fifth solution in results3
plot(results2, results3, i=3, j=5)

# geoplot showing difference in selection frequencies between the two objects
# white colors indicate that units are already in a protected area
# blue colours indicate that units were more often selected in results2
# red  colours indicate that units were more often selected in results3
plot(results2, results3)

# All the solutions in our current portfolio, results3, seem to be fairly fragmented. 
# If implemented as protected areas, these solutions might be associated with poor connectivity 
# and high maintenance costs. To reduce fragmentation, we can increase the boundary length
# modifier (BLM). However, in order to maintain adequate levels of representation for the 
# species, MARXAN will select more expensive planning units. 
# How can we pick an appropriate BLM while still making sure the acquisition costs are adequate cost? 
# Let's generate six more portfolios, each using a different BLM, and plot the trade-off 
# between acquisition cost and fragmentation using the best solutions in each portfolio.

## generate list of portfolios with different BLMS
# make vector BLM parameters to use
blm.pars=c(0, 100, 250, 500, 750, 1000)

# create list with different portfolio for each BLM
results4<-list()
for (i in seq_along(blm.pars)) {
  results4[[i]]<-update(results3, ~opt(BLM=blm.pars[i], NUMREPS=10L))
}

## extract data from portfolios
# create empty vectors to store values
cost<-c()
con<-c()
blm<-c()

# extract values for best solutions
for (i in seq_along(blm.pars)) {
  cost<-append(cost, summary(results4[[i]])[["Cost"]])
  con<-append(con, summary(results4[[i]])[["Connectivity"]])
  blm<-append(blm, rep(blm.pars[i], nrow(summary(results4[[i]]))))
}

## plot trade-off between shortfall and connectivity
# get colours for legend
legend.cols<-c("#FFFFB2", "#FED976", "#FEB24C", "#FD8D3C", "#F03B20", "#BD0026")
pt.cols<-legend.cols[match(blm, blm.pars)]

# reset plotting window
par(mfrow=c(1,1))

# plot trade-off data
# higher shortfall values means worse representation
# higher connectivity values mean more fragmentation
plot(cost~con, bg=pt.cols, col='black', ylab='Cost', xlab='Connectivity', pch=21,
     main='Trade-off between cost and connectivity')
abline(lm(cost~con))

# add legend
legend("topright", legend=blm.pars, col='black', pt.bg=legend.cols, pch=21, title='BLM')

# Looking at this curve and depending on the total budget, you might decide that second portfolio 
# in results4 achieves an acceptable level of representation and fragmentation. 
# Let's generate another portfolio of solutions with BLM=0.0001, and then 
# make some geoplots to compare it to results3.


# make new solutions with BLM=0.0001
results5<-update(results3, ~opt(BLM=0.0001))

# geoplot showing differences between the best solution in each portfolio
plot(results5, results3, i=0, j=0)

# geoplot showing difference in selection frequencies between the two objects
# black colours indicate that units are already in a protected area
# blue colours indicate that units were more often selected in results4[[2]],
# and red colours indicate that they were often selected in results3
plot(results5, results3)

# So now, in our final portfolio, we have one hundred solutions. 
# How can we compare them all and decide on a final prioritisation to implement? 
# We don't time time to make 100 maps; while the maps this R package makes are pretty, 
# they are not that pretty. Instead, we could make some dotcharts that let us compare various 
# properties of the solutions.

# make dotchart showing the score of each solution
# the score describes the overall value of the prioritisations based on our criteria
# the lower the value, the better the solution
# the best solution is coloured in red
dotchart(results5, var='score')

# make dotchart showing the connectivity of the solutions
# solutions with lower values are more clustered
# solutions with higher values are more fragmented
# argument to n specifies the number of solutions to plot
# argument to nbest specifies number of solutions to colour in red
dotchart(results5, var='con', nbest=5, n=50)


# How can we summarise the main themes of variation in our portfolio?
# Fortunately, statisticians solved this problem a long time ago. 
# We can use ordination techniques to create a few variables that describe commonalities 
# among the solutions, and visualise the main sources of variation in a small number of 
# dimensions.

## dendrogram showing differences between solutions based on which planning units 
## were selected (using Bray-Curtis distances by default)
# the solutions are shown at the (bottom) tips of the tree.
# solutions that occupy nearby places in tree
# have similar sets of planning units selected.
# the best prioritisation is coloured in red.
dendrogram(results5, type='dist', var='selections')

## same dendrogram as above but with the best 10 prioritisations coloured in red
# if all the red lines connect together at the bottom of the dendrogram
# this means that all the best prioritisations are really similar to each other,
# but if they connect near the top of the dendrogram then this means that
# some of the best prioritisations have totally different sets of planning units
# selected for protection.
dendrogram(results5, type='dist', var='selections', nbest=10)

## ordination plot showing differences between solutions based on the number of units
## occupied by each vegetation class (using MDS with Bray-Curtis distances)
# we can also use multivariate techniques to see how the solutions vary
# based on how well they represent different vegetation classes.
# the numbers indicate solution indices.
# solutions closer to each other in this plot have more
# similar levels of representation for the same species.
# the size of the numbers indicate solution quality,
# the bigger the number, the higher the solution score.
ordiplot(results5, type='mds', var='occheld', method='bray')

# ordination plot showing differences between solutions based on the amount held 
# by each vegetation class (using a principle components analysis)
# labels are similar to the previous plot.
# the arrows indicate the variable loadings.
ordiplot(results5, type='pca', var='amountheld')




