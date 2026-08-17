######
### High res data
### Author: DHare
### contact dhare@umass.edu

Start_time <- Sys.time()

setwd("Proposal_2024/R")

#install.packages("nhdplusTools")

library(igraph)
library(tidyverse)
library(sf)

library(sp)
library(dataRetrieval)
library(nhdplusTools)
#Plotting Functions
library(cowplot)
library(scales)
#Using Larger Model for Moss for Automatic 
sf_use_s2(FALSE)
#Steele, Distribution of moss percentages
moss_p <- data.frame(percent_moss = c(0, 41, 36, 36, 9, 10, 41, 24, 20))
# prob_m <- prop.table(moss_p) #find probability 
# sample(moss_p$percent_moss,1000,replace=TRUE, prob = prob_m$percent_moss)
# ggplot(moss_p)+
#   geom_density(aes(x =percent_moss))

#create continuous density plot
pdf_of_data <- density(moss_p$percent_moss, from= 0, to=60, bw=10)
#plot for review
plot(pdf_of_data)
#create reach distirbution example 
random.points <- approx(
  cumsum(pdf_of_data$y)/sum(pdf_of_data$y),
  pdf_of_data$x,
  runif(10000)
)$y
#plot random example for review. 
hist(random.points, 100)

ggsave("histogram_streamMossDistirbution.png")

n_loading_rate = 0.001
#Read in nearest USGS station for ease
site <- dataRetrieval::readNWISsite("01075000")
#193 sq miles = 500 km2 = 50000 ha
start_point <- st_sfc(st_point(c(site$dec_long_va, site$dec_lat_va)), crs = 4326)
start_comid <- discover_nhdplus_id(start_point)


hu <- get_huc(start_point, type = "huc08") 
#(hu <- substr(hu$huc8, 1, 2)) download_nhdplushr(tempdir(), c(hu, "0203"), download_files = FALSE)

flowline <- navigate_nldi(list(featureSource = "comid", 
                               featureID = start_comid), 
                          mode = "upstreamTributaries", 
                          distance_km = 1000)

ggplot()+
  geom_sf(data = sf::st_geometry(flowline$UT_flowlines))


#if doesnt work redownload NHDPlusTools
subset_file <- tempfile(fileext = ".gpkg")
subset <- subset_nhdplus(comids = as.integer(flowline$UT_flowlines$nhdplus_comid),
                         output_file = subset_file,
                         nhdplus_data = "download", 
                         flowline_only = FALSE,
                         return_data = TRUE, overwrite = TRUE)

#saveRDS(subset, "Upper_Pecw_HB.RDS")
subset <- readRDS("Upper_Pecw_HB.RDS")


Total_WatershedArea <- sum(subset$NHDFlowline_Network$areasqkm) # matches the online resources 500 km2 == 506


#calculate_total_drainage_area(sf::st_set_geometry(subset$NHDFlowline_Network, NULL))

#how many stream reaches are 1-2 order 
n_Smallstreams = nrow(dplyr::filter(subset$NHDFlowline_Network, streamorde < 3))
set.seed(123)
random.mossD <- approx(
  cumsum(pdf_of_data$y)/sum(pdf_of_data$y),
  pdf_of_data$x,
  runif(n_Smallstreams)
)$y

perc_moss_reach = replace_na(c(random.mossD, rep(0, nrow(subset$NHDFlowline_Network) - n_Smallstreams)), 0)

# Create the node metadata for igraph
fline_nodes <- subset$NHDFlowline_Network %>%
  as.data.frame() %>%
  dplyr::select(id = hydroseq,
         to = dnhydroseq,
         comid,
         vpuid,
         length_km = lengthkm,
         stream_order = streamorde,
         total_drain_area_km = totdasqkm,
         area_km = areasqkm,
         qa_ma,
         nhd_id = comid,
         reachcode = reachcode)%>%
  arrange(stream_order,length_km)%>% #added length so sort is consistent with nodes and lines
  #how do I just apply the random distribution to a subset of values while keeping the distribution 
  #my approach right now is to create order by stream order, then knowing the number of 1 & 2 already (n_smallstreams)
  #cbind that list with zeros at the end. I imagine there is a more sufficated way. 
  mutate(
    perc_moss_reach = perc_moss_reach
    #perc_moss_reach = 0
    # perc_moss_reach = replace(perc_moss_reach, sample(row_number(),  
    #                         size = ceiling(0.36 * n()), replace = FALSE), 1) #0.36 # ifelse(stream_order==1, 0.36, 0), #percent moss coverage each reach 
    
  )
  
fline_link <- subset$NHDFlowline_Network %>%
  as.data.frame() %>%
  dplyr::select(from = hydroseq,
         id = hydroseq,
         to = dnhydroseq,
         comid,
         vpuid,
         length_km = lengthkm,
         stream_order = streamorde,
         total_drain_area_km = totdasqkm,
         area_km = areasqkm,
         qa_ma,
         nhd_id = comid,
         reachcode = reachcode,
         -geometry)%>%
  arrange(stream_order,length_km)%>% #stream order has to be first
  mutate(
    perc_moss_reach = perc_moss_reach)
    #perc_moss_reach = 0)
  # mutate(
  #   perc_moss_reach = 0,
  #   perc_moss_reach = replace(perc_moss_reach, sample(moss_p$percent_moss,row_number(),replace=TRUE), 1)#sample(row_number(),  
  #                                                     #size = ceiling(0.36 * n()), replace = FALSE), 1) #0.36 # ifelse(stream_order==1, 0.36, 0), #percent moss coverage each reach 
  #                                                     
  #   )
#find outlet data specific
outlet_code <- subset$NHDFlowline_Network$dnhydroseq[!subset$NHDFlowline_Network$dnhydroseq %in% subset$NHDFlowline_Network$hydroseq]
#network without outlet
network_links <- fline_link[!fline_link$to %in% outlet_code,]

# Convert NHD to igraph
net <- graph_from_data_frame(d=network_links, vertices=fline_nodes, directed=T)

# Save file
write_rds(net, paste0('../data/nhd_files/igraph_network_PewcHR.rds'))


#############################################################################################################################
### set up network N model HR
#############################################################################################################################

#keep original
network <- net
# Convert igraph to dataframe
net_data_frame <- igraph::as_data_frame(network)%>%
  mutate(id = from) 

# Get the outlet of the basin
last_flowline <- net_data_frame %>%
  filter(total_drain_area_km == max(total_drain_area_km)) %>%
  pull(from)

all_ids <- net_data_frame$id

## ---Function identifies all nodes upstream of a selected node--- ##
get_number_upstream_nodes <- function(ids){
  
  network_connections <- tibble(id=all_ids,
                                upstream_connections = NA)
  get_length <- function(ids) {
    length(head(unlist(ego(network,order=length(V(network)),nodes=as.character(ids),mode=c("in"),mindist=0)),-1))
  }
  
  connections <- unlist(map(ids, get_length))
  
  network_connections$upstream_connections <- connections
  
  return(network_connections)
}
network_order <- get_number_upstream_nodes(all_ids) %>%
  arrange(upstream_connections)

# network_order <- rbind(network_order, tibble(id = as.double(last_flowline), upstream_connections = NA))


## ------------- Function for reach Q ------------------- ##
q_catchment <- function(area, yeild = (7.69e-09 * 86400)){
  #  7.69 × 10−9 m3 m−2 s−1 is the yield from Helton et al. 2017
  # Area must be in km2
  
  # Convert from km to m2 and multiply times yield to get m3/d of water
  q_catch <- area*1e+6*yeild
  
  return(q_catch)
}

## ------------- Function for reach N ------------------- ##
n_catchemnt <- function(area, yeild = n_loading_rate){
  # 0.0001 to 100 kg N km−2 day−1 from helton et al 2017
  # Area must be in km
  
  n_catch <- area*yeild
  
  return(n_catch)
}

## ------------- N removal Model ------------------------------- ##
n_uptake_velocity <- function(discharge, n_load, c, d, SurfArea, perc_moss_reach){

  #adding "normal denitrification" + moss 
    Vfp <- c * (n_load/discharge)^d + 
    perc_moss_reach * SurfArea * 0.1 * 1E-5 * 24 #kgN/m2/d (1ug/cm2 = 1E-5 kg/m2)
    #STEELE REMOVE 0.` - 0.1 ugN-NO3/cm2/hr
  
}




n_uptake <- function(n_load, uptake_velocity, surface_area, discharge){
  
  # Ratio of discharge to streambed surface area (width * length)
  Hlp <- discharge/surface_area
  
  Nrp = n_load * (1 - exp(1)^(-uptake_velocity/Hlp))
}


## ---------- Functions used in the mass balance model --------- ##

# Function to implement the N mass-balance model:
solveMB_N <- function(network, network_order){
  ##debug
  #network <- net#net_BF
  # Define input parameters:
  
  # Inflow discharge from local catchment (m3/d):
  V(network)$Qlocal <- 0
  # Inflow discharge from upstream reaches (m3/d):
  V(network)$Qnet <- 0
  # Outflow discharge from each reach (m3/d):
  V(network)$Qout <- NA
  # N load from upstream catchment (kg/d):
  V(network)$Nload <- 0
  # N load from upstream reaches (kg/d):
  V(network)$Nnet <- 0
  # Total N in reach (kg/d)
  V(network)$Nloc <- 0
  # Surface Area (m2)
  V(network)$SurfArea <- 0
  # N Uptake Velocity (m/d)
  V(network)$Nv <- 0
  # N removed (kg/d)
  V(network)$Nlost <- NA
  # Exported N load from each reach (kg/d):
  V(network)$Nout <- 0

  
  #hedrology/n removal Parameters Helton 2017
  a=7.17
  b=0.348
  c=3.196e-4
  d=-0.49311
  
  # Calculate mass-balance for each reach moving down the network from headwaters to mouth:
  # i = 0
  for(i in 1:length(network_order$id)){
    #for(i in 1:52){
    
    # Get location in igraph network
    p <- grep(network_order$id[i], names(V(network)))
    n_ID <- as.character(network_order$id[i]) #node ID being run
    e_ID <- incident(network, n_ID, mode = c("in")) #get edges so the data is shared - this should replace the mutate to nodes.
    #e_ID$SEGMENT_ID
    #e_ID <-
    # length_reach <- sum(n_ID$length_km)
    # length_reach <- ifelse(length_reach == 0, 10, length_reach) # locations within any edges in are springs and thus have 10 meter contributing
    # V(network)$length_reach[i] <- length_reach
    
    # Find neighboring reaches upstream that flow in to the reach:
    up <- igraph::neighbors(network, p, mode=c("in")) #only single up
    up.all.nodes <- head(unlist(ego(network,order=length(V(network)),nodes=n_ID,mode=c("in"),mindist=0)),-1)
    
    # Define hydrologic inflows/outflows for each reach (m3 d-1)
    # Discharge inflow from local catchment (m3 d-1):
    #V(network)$Qlocal[i] <- V(network)$runoff_mday[i] * (V(network)$areasqkm[i]*10^6)
    
    if(V(network)$area_km[p] == 0){
      V(network)$area_km[p] <- 0.1
    }
    
    # m3/d
    V(network)$Qlocal[p] <-  q_catchment(V(network)$area_km[p])
    
    # Discharge inflow from upstream network (m3 d-1):
    if(length(up)>0){
      V(network)$Qnet[p] <- sum(V(network)$Qout[up]) #only one or junction will have two
      V(network)$Qout[p] <- sum(V(network)$Qlocal[p], V(network)$Qnet[p], na.rm = T)
    }else{ #if not inputs then the Q is just 0 to the local input
      V(network)$Qout[p] <- V(network)$Qlocal[p]
    }
    
    #w = aQb
    # Discharge outflow to downstream reach (m3 d-1):

    
    # Q must be in m3/s
    width_m <- a * ((V(network)$Qout[p])/86400)^b
    network <- set_vertex_attr(network, "width_m", p, width_m)
    #V(network)$width_m[p] =  width_m #wasnt working 10/17/2024, could it be because object name is the same as the colum name? Anyway this adjustment works
    
    # Calc surface area and convert length in km to m (m2)
    V(network)$SurfArea[p] <-  V(network)$width_m[p] * (V(network)$length_km[p]) * 1000
    
    # Define nitrogen inflows/outflows for each reach (kg d-1)
    # nitrogen inflow from local catchment (kg d-1):
    V(network)$Nload[p] <- n_catchemnt(V(network)$area_km[p])
    
    # Nitrogen inflow from upstream network (kg d-1):
    if(length(up)>0){
      V(network)$Nnet[p] <- sum(V(network)$Nout[up], na.rm = T)
    }
    
    # Total N in reach
    V(network)$Nloc[p] <- V(network)$Nload[p] +  V(network)$Nnet[p]
    
    ## Nitrogen Removed
    # N uptake velocity
    V(network)$Nv[p] <- n_uptake_velocity(V(network)$Qout[p],  V(network)$Nloc[p], c = c, d = d,
                                          V(network)$SurfArea[p], V(network)$perc_moss_reach[p])
    
    # Removed from channel
    n_removed <- n_uptake(n_load = V(network)$Nloc[p],
                          uptake_velocity =  V(network)$Nv[p],
                          surface_area = V(network)$SurfArea[p],
                          discharge =  V(network)$Qout[p])
    
    V(network)$Nlost[p] <- n_removed
    
    # Total N exported downstream
    V(network)$Nout[p] <- V(network)$Nloc[p] - V(network)$Nlost[p]
    
  }
  
  # Get list with attributes
  out <- get.vertex.attribute(network)
  
  # Export network:
  return(out)
  
}


########################################################################################################################
### model run
########################################################################################################################

## J.Steele 2019 Values ##
# NO3 uptake rate: 2.3 g N y-1
# percent reach range: 9.1%–46.7%
# percent network range: 4%–40%
# median coverage of that stream is around 36% (Vought et al. 2019).
#By scaling nitrate assimilation rates from our jar incubations using our estimates oftotal stream bryophyte mass,
#we estimated whole-stream N uptake rates due to bryophyte
# #associated nitrate assimilation to be 2.3 g N m-2 y-1 for Bear Brook and 1.4 g N m-2 y-1 for 327 Paradise Brook.
# By scaling the average, 1-day nitrate uptake from 442 our Bear Brook bryophyte incubations to the stream-level bryophyte 
# cover, we estimate that 443 bryophyte-associated NO3- uptake (2.3 g N m-2 yr-1) is similar in
# magnitude to previous estimates 444 of whole-stream NO3- uptake rates in this stream (12 g N m-2 yr-1) (Bernhardt et al. 2003).

#Areal Specfic 
#0.1 - 0.2 ug N/cm2/hr

#surveyed 50m reaches

n_loading_rate <- 0.001
#perc_moss <-  0.36
  
# Run network model
n_network <- solveMB_N(network = network,
                       network_order = network_order)

n_net_data <- as.data.frame(do.call(cbind, n_network))

n_net_data <- as.data.frame(lapply(n_net_data, unlist, recursive = TRUE))

# Add in final flowline # For some reason, the last nodel has to be treated separately, look at iinital code if desired

write_rds(n_net_data, paste0("../data/nhd_files/model_output_NLoad", n_loading_rate,'Moss_SteeleDistribution_1_2.rds'))
#write_rds(n_net_data, paste0("../data/nhd_files/model_output_NLoad", n_loading_rate,'Moss0.rds'))



#######################################################################################################################
### plotting
#######################################################################################################################

color_set = c("335C67", "yellow", "orange", "9E2A2B", "540B0E" )
#Hydrologic Model NHD Geometry/Framework
fline_og <- subset$NHDFlowline_Network %>%
  as.data.frame() ## has geometry

#Moss Model 
model_output <- read_rds(paste0("../data/nhd_files/model_output_NLoad0.001Moss_SteeleDistribution_1_2.rds")) %>% 
  mutate_at(5:21, as.numeric)%>%
  #mutate_if(is.character,as.numeric, -vp)%>% #this does all, which is okay, but first few columns can be char for now
  mutate(frac_removed = (Nlost/Nloc)*100)%>%
  dplyr::filter(!is.na(width_m)) #removed outlet as has NAs

# flowlines <- read_rds(paste0('data/nhd_files/', 6109553, '/flowlines.rds')) %>%
#   select(comid, geometry)



n_net_data_sf <- merge(fline_og, model_output, by.x = "comid", by.y = "nhd_id")

n_net_data_sf$Nout_g <- cut(n_net_data_sf$Nout *365,
               breaks = c(-Inf, 1, 5, 10, 25, Inf),
               labels = c("<1","1-5", "5-10", "10-25", ">25"))

theme_set(theme_map())
##%percent 
p_NperRemoved <- ggplot()+
                    geom_sf(data = n_net_data_sf$geometry, 
                            aes(color = n_net_data_sf$frac_removed), linewidth = 1.25)+
                    geom_sf(data = start_point, shape = 21)+
  geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.015, color = "black")+
  scale_color_viridis_c(option = "G", limits = c(0, 100), breaks = seq(25,100, by =25),
                        labels = c("25", 50, 75, 100), begin=0.2, end = 0.8, direction = -1)+
                    labs(x = "", y = "", color = "% N Load Removed")+
  theme(
    axis.text.x=element_blank()#,
    #panel.background = element_rect(fill = 'black')
  )
             
##amount exported
p_Nexport <- ggplot()+
                geom_sf(data = n_net_data_sf$geometry, 
                        aes(color = n_net_data_sf$Nout_g), linewidth = 1.25)+
                geom_sf(data = start_point, shape = 21)+
                geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.015, color = "black")+
  scale_color_viridis_d(option = "F", begin = 0.3, end = 0.9)+
                labs(subtitle = "Aquatic Bryophytes", x = "", y = "", color = "")+ #subtitle = paste("N exported (kg/d), N load =", n_loading_rate)
  theme(   

    axis.text.x=element_blank()#,
    #panel.background = element_rect(fill = 'black')
  )+
   guides(color=guide_legend(title="N exported (kg/yr)"))
p_Nexport

n_net_data_0moss <- read_rds(paste0("../data/nhd_files/model_output_NLoad0.001Moss0.rds"))%>%
  mutate_at(5:21, as.numeric)%>%
  #mutate_if(is.character,as.numeric, -vp)%>% #this does all, which is okay, but first few columns can be char for now
  mutate(frac_removed = (Nlost/Nloc)*100)%>%
  dplyr::filter(!is.na(width_m)) 
n_net_data_0moss <- merge(fline_og, n_net_data_0moss, by.x = "comid", by.y = "nhd_id")

n_net_data_0moss$Nout_g <- cut(n_net_data_0moss$Nout *365,
                            breaks = c(-Inf, 1, 5, 10, 25,Inf),
                            labels = c("<1","1-5", "5-10", "10-25", ">25"))

p_NperRemoved_0moss <- ggplot()+
  geom_sf(data = n_net_data_0moss$geometry, 
          aes(color = n_net_data_0moss$frac_removed), linewidth = 1.25)+
  geom_sf(data = start_point, shape = 21)+
  geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.015, color = "black")+
  scale_color_viridis_c(option = "G", limits = c(0, 100), breaks = seq(25,100, by =25),
                        labels = c("25", 50, 75, 100), begin=0.2, end = 0.8, direction = -1)+
  labs(subtitle = "", x = "", y = "", color = "% N Load Removed")+
  theme(
   legend.position = "none",
    axis.text.x=element_blank()
  #  panel.background = element_rect(fill = 'black')
  )

#n_net_data_0moss <- dplyr::filter(n_net_data_0moss, lengthkm > 0.1)
p_Nexport_nomoss <- ggplot()+
  geom_sf(data =n_net_data_0moss$geometry, 
          aes(color = n_net_data_0moss$Nout_g), linewidth = 1.25)+
  geom_sf(data = start_point, shape = 21)+
  geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.015, color = "black")+
  scale_color_viridis_d(option = "F", begin = 0.3, end = 0.9)+
  labs(subtitle = "No Aquatic Bryophytes", x = "", y = "", color = "")+ #subtitle = paste("N exported (kg/d), N load =", n_loading_rate)
  theme(
    legend.position = "none",
    axis.text.x=element_blank()
   # panel.background = element_rect(fill = 'white')
  )+
  guides(color=guide_legend(title="N exported (kg/yr)"))
p_Nexport_nomoss

# These are the output values used in the proposal 

max(n_net_data_0moss$Nout *365)#90.6
max(n_net_data_sf$Nout * 365)#28.17

max(n_net_data_0moss$Nout *365)/(Total_WatershedArea * 100) #(ha) 
max(n_net_data_sf$Nout * 365)/ (Total_WatershedArea * 100) #(ha)


max(n_net_data_0moss$Nout *365)/(Total_WatershedArea) #(km2) 
max(n_net_data_sf$Nout * 365)/ (Total_WatershedArea) #(km2)


max(n_net_data_0moss$Nout *365)/(Total_WatershedArea)- max(n_net_data_sf$Nout * 365)/ (Total_WatershedArea)#(km2) 
 #(km2)

n_loading_rate





#library(cowplot)

plot_grid(p_Nexport_nomoss, p_Nexport,
          p_NperRemoved_0moss, p_NperRemoved,
          align = "hv", nrow = 2)


ggsave(paste0("C:/Users/hared/Dropbox/ResearchProjects/100_HubbardBrookProjects/1010_Mosses/Proposal_2024/figs/Network_Plots_CompareMoss.png"), width = 6, height = 6, units = "in")
ggsave(paste0("C:/Users/hared/Dropbox/ResearchProjects/100_HubbardBrookProjects/1010_Mosses/Proposal_2024/figs/Network_Plots_CompareMoss.svg"), width = 6, height = 6, units = "in")
plot_grid(p_NperRemoved_0moss , p_Nexport_nomoss, nrow = 2, align = "hv")
ggsave(paste0("C:/Users/hared/Dropbox/ResearchProjects/100_HubbardBrookProjects/1010_Mosses/Proposal_2024/figs/Network_Plots_Nload_Moss36pRandom1.png"), width = 6, height = 8, units = "in")

### Compare Plots

model_compare_12_Steele <- read_rds(paste0("../data/nhd_files/model_output_NLoad0.001Moss0.rds"))%>%
  left_join(read_rds(paste0("../data/nhd_files/model_output_NLoad0.001Moss_SteeleDistribution_1_2.rds")), by = "nhd_id")%>%
  mutate_if(is.character,as.numeric)%>%
  mutate(
    Nout_pdiff = (Nout.y - Nout.x)/abs(Nout.x) * 100
  )%>%
  dplyr::filter(!is.na(width_m.y)) #removed outlet as has NAs


model_compare_36full <- read_rds(paste0("../data/nhd_files/model_output_NLoad0.001Moss0.rds"))%>%
  left_join(read_rds(paste0("../data/nhd_files/model_output_NLoad0.001Moss36pFull.rds")), by = "nhd_id")%>%
  mutate_if(is.character,as.numeric)%>%
  mutate(
    Nout_pdiff = (Nout.y - Nout.x)/abs(Nout.x) * 100
  )

model_compare_36_1stonly <- read_rds(paste0("../data/nhd_files/model_output_NLoad0.001Moss0.rds"))%>%
  left_join(read_rds(paste0("../data/nhd_files/model_output_NLoad0.001Moss0361st.rds")), by = "nhd_id")%>%
  mutate_if(is.character,as.numeric)%>%
  mutate(
    Nout_pdiff = (Nout.y - Nout.x)/abs(Nout.x) * 100
  )


model_compare_36_random <- read_rds(paste0("../data/nhd_files/model_output_NLoad0.001Moss0.rds"))%>%
  left_join(read_rds(paste0("../data/nhd_files/model_output_NLoad0.001Moss36pRandom1.rds")), by = "nhd_id")%>%
  mutate_if(is.character,as.numeric)%>%
  mutate(
    Nout_pdiff = (Nout.y - Nout.x)/abs(Nout.x) * 100,
    perc_moss_reach.y = ifelse(length_km.y <0.013,0.36, perc_moss_reach.y) #for legend
  )



p_percDiff_12_Steele <- ggplot()+
                            geom_sf(data = n_net_data_sf$geometry, 
                                    aes(color = model_compare_12_Steele$Nout_pdiff), linewidth = 1.25)+
                            geom_sf(data = start_point, size = 0.5, shape = 21)+
                            geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.015)+
                            scale_color_viridis_c(direction =-1, option = "C", limits = c(NA, 0), oob = scales::squish)+
                            labs(subtitle =  paste("D"), x = "", y = "", color = "")+
  theme(
    axis.text.x=element_blank(),
    legend.position = "right",
  )+
  guides(color=guide_legend(title="% Change"))

p_percDiff_12_Steele_mossD <- ggplot()+
  geom_sf(data = n_net_data_sf$geometry, 
          aes(color = model_compare_12_Steele$perc_moss_reach.y), linewidth = 1.25)+
  geom_sf(data = start_point, shape = 21)+
  scale_color_gradient(low = "grey", high = "darkgreen", na.value = NA)+
  geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
  labs(#title ="Observed Moss Distribution Applied to 1st and 2nd Streams",
       subtitle = "A", #
                x = "", y = "", color = "")+
  theme(
    axis.text.x=element_blank(),
    legend.position = "right"
  )+
  guides(color=guide_legend(title="% Moss Coverage"))

 n_loading_rate
 

 
 p_percDiff_12_Steele_mossD <- ggplot()+
   geom_sf(data = n_net_data_sf$geometry, 
           aes(color = model_compare_12_Steele$perc_moss_reach.y), linewidth = 1.25)+
   geom_sf(data = start_point, shape = 21)+
   scale_color_gradient(low = "grey", high = "darkgreen", na.value = NA)+
   geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
   labs(#title ="Observed Moss Distribution Applied to 1st and 2nd Streams",
     subtitle = "A", #
     x = "", y = "", color = "")+
   theme(
     axis.text.x=element_blank(),
     legend.position = "right"
   )+
   guides(color=guide_legend(title="% Moss Coverage"))
 
 
ggplot()+
   geom_sf(data = n_net_data_sf$geometry, 
           aes(color = model_compare_12_Steele$Nout.x), linewidth = 1.25)+
   geom_sf(data = start_point, shape = 21)+
   scale_color_gradient(low = "grey", high = "darkgreen", na.value = NA)+
   geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
   labs(#title ="Observed Moss Distribution Applied to 1st and 2nd Streams",
     subtitle = "A", #
     x = "", y = "", color = "")+
   theme(
     axis.text.x=element_blank(),
     legend.position = "right"
   )+
   guides(color=guide_legend(title="N exported (kg/yr)"))
 
 
  ggplot()+
    geom_sf(data = n_net_data_sf$geometry, 
            aes(color = model_compare_12_Steele$Nout.y), linewidth = 1.25)+
    geom_sf(data = start_point, shape = 21)+
    scale_color_gradient(low = "grey", high = "darkgreen", na.value = NA)+
    geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
    labs(#title ="Observed Moss Distribution Applied to 1st and 2nd Streams",
      subtitle = "A", #
      x = "", y = "", color = "")+
    theme(
      axis.text.x=element_blank(),
      legend.position = "right"
    )+
    guides(color=guide_legend(title="N exported (kg/yr)"))
  
 # ggplot()+
 #   geom_sf(data = n_net_data_sf$geometry, 
 #           aes(color = model_compare_12_Steele$stream_order.y), linewidth = 1.25)+
 #   geom_sf(data = start_point, size = 0.5)+
 #   geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
 #   labs(subtitle =  paste("Observed Moss Distribution \n 1st and 2nd Only"), x = "", y = "", color = "")+
 #   theme(
 #     axis.text.x=element_blank(),
 #     legend.position = "right"
 #   )
 
 
 plot_grid(p_percDiff_12_Steele_mossD, p_Nexport,p_NperRemoved, 

           #p_percDiff_12_Steele,
           nrow = 3, align = "hv")
 
 ggsave("Network_Plots_Nload_Moss_Steele1-2.png", width = 4, height = 6, units = "in")
 
 

p_percDiff_36p1st <- ggplot()+
                              geom_sf(data = n_net_data_sf$geometry, 
                                      aes(color = model_compare_36_1stonly$Nout_pdiff), linewidth = 1.25)+
                              geom_sf(data = start_point, size = 0.5)+
                              geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
                              scale_color_viridis_c(direction =-1, option = "C", limits = c(-50, 0), oob = scales::squish)+
                              labs(subtitle = "36% Bryophyte Coverage in 1st Order Streams", x = "", y = "", color = "")+
                              theme(
                                axis.text.x=element_blank(),
                                legend.position = "none"
                              )

p_percDiff_36random <- ggplot()+
                            geom_sf(data = n_net_data_sf$geometry, 
                                    aes(color = model_compare_36_random$Nout_pdiff), linewidth = 1.25)+
                            geom_sf(data = start_point, size = 0.5)+
                            geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
                            scale_color_viridis_c(direction =-1, option = "C", limits = c(-50, 0), oob = scales::squish)+
                            labs(subtitle = "Random Coverage Placement Across Network", x = "", y = "", color = "")+
                            theme(
                              axis.text.x=element_blank(),
                              legend.position = "right"
                            )


plot_grid(p_percDiff_12_Steele, #p_percDiff_36full,  p_percDiff_36p1st, p_percDiff_36random, align = "hv",
          title = "Percent Change from Background N Export")
ggsave(paste0("Network_Plots_perNexportChange.png"), width = 6, height = 4, units = "in")


p_percMcover_36full <- ggplot()+
  geom_sf(data = n_net_data_sf$geometry, 
          aes(color = as.factor(model_compare_36full$perc_moss_reach.y)), linewidth = 1.25)+
  geom_sf(data = start_point, size = 0.5)+
  geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
  scale_color_manual(values = c("0" = "#dad7cd", "0.36" = "#606c38", "1" = "#3a5a40"))+
  labs(subtitle =  paste("Coverage Across All Reaches"), x = "", y = "", color = "")+
  theme(
    axis.text.x=element_blank(),
    legend.position = "none"
  )

p_percMcover_36p1st <- ggplot()+
  geom_sf(data = n_net_data_sf$geometry, 
          aes(color = as.factor(model_compare_36_1stonly$perc_moss_reach.y)), linewidth = 1.25)+
  geom_sf(data = start_point, size = 0.5)+
  geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
  scale_color_manual(values = c("0" = "#dad7cd", "0.36" = "#606c38", "1" = "#3a5a40"))+
  labs(subtitle = " Bryophyte Coverage in 1st Order Streams Only", x = "", y = "", color = "")+
  theme(
    axis.text.x=element_blank(),
    legend.position = "none"
  )

p_percMcover_36random <- ggplot()+
  geom_sf(data = n_net_data_sf$geometry, 
          aes(color = as.factor(model_compare_36_random$perc_moss_reach.y)), linewidth = 1.25)+
  geom_sf(data = start_point, size = 0.5)+
  geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
  scale_color_manual(values = c("0" = "#dad7cd", "0.36" = "#606c38", "1" = "#3a5a40"))+
  #scale_color_viridis_d(direction =-1, option = "A")+
  labs(subtitle = "Coverage Random Across Network", x = "", y = "", color = "")+
  theme(
    axis.text.x=element_blank()
  )

plot_grid(p_percMcover_36full,  p_percMcover_36p1st,p_percMcover_36random, align = "hv",
          axis = "bt",
          title = "Percent Change from Background N Removal")

ggsave(paste0("Network_Plots_perMoss.png"), width = 6, height = 4, units = "in")

model_compare_36_random[which.max(model_compare_36_random$qa_ma.x),]$



ggplot(model_compare_36_1stonly)+
  geom_density(aes(Nout_pdiff))+
  theme_bw()

ggplot(model_compare_36_random)+
  geom_density(aes(Nout_pdiff))+
  theme_bw()

ggplot(model_compare_36full)+
  geom_density(aes(Nout_pdiff))+
  theme_bw()


ggplot()+
  geom_sf(data = n_net_data_sf$geometry, 
          aes(color = as.factor(model_compare_36_random$perc_moss_reach.y)), linewidth = 1.25)+
  geom_sf(data = start_point)+
  geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
  scale_color_viridis_d(direction =-1, option = "G")+
  labs(title = "Random Moss Locations", x = "", y = "", color = "")+
  theme(
    axis.text.x=element_blank()
  )



p_StreamOrder <- ggplot()+
                        geom_sf(data = n_net_data_sf$geometry, 
                                aes(color = model_compare_36_random$stream_order.x), linewidth = 1.25)+
                        geom_sf(data = start_point)+
                        geom_sf_text(data = start_point, aes(label = "Outlet"), nudge_y = -0.01)+
                        scale_color_viridis_c(direction =-1)+
                        labs( x = "", y = "", color = "")+
                        theme(
                          axis.text.x=element_blank(),
                          legend.position = "none"
                        )


#percent of streams count 1st order
ggplot(model_compare_36_random)+
  geom_density(aes(stream_order.x))+
  theme_bw()

#percent of stream length 1st order
total_net_SA <- sum(model_compare_36_random$SurfArea.x)

p_surfaceArea <- ggplot(model_compare_36_random)+
  geom_col(aes(stream_order.x, SurfArea.x/total_net_SA, fill = stream_order.x))+
  theme_bw()+
  scale_fill_viridis_c(direction =-1)+
  scale_y_continuous(labels = scales::percent)+
  labs(x = "Stream Order", y = "% of Network Surface Area")+
  theme(
    legend.position = "none"
  )


plot_grid(p_StreamOrder, p_surfaceArea, nrow = 1, align = "hv", rel_widths = c(2, 1))

ggsave(paste0("Network_Plots_StreamOrder.png"), width = 6, height = 4, units = "in")





#######################################################################################################
### Try to get NLCD for each stream segment

library(FedData)
library(raster)
library(Matrix)
library(data.table)
library(plyr)


df <- st_as_sf(n_net_data_sf, wkt = "geom")
df$midpoint <- st_centroid(df)

bbox <- st_bbox(df)
# bbox <- list(bbox)
bbox_df <- data.frame(lat = c(bbox[2] - .2 , bbox[4] + .2),
                      lon = c(bbox[1] + .5, bbox[3] - .5))
coordinates(bbox_df) <- ~lon + lat
proj4string(bbox_df) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

buffer_dist <- 300

nlcd_raster <- get_nlcd(bbox_df, label = "NLCD_2019", year = 2019, extraction.dir = '.')
ssurgo <- get_ssurgo(bbox_df, label = "ssurgo",
                     raw.dir = '.',
                     extraction.dir = '.')



# projectRaster(nlcd_raster, crs = projection(df))


# landcover <- data.frame()

coords <- do.call(rbind, st_geometry(df$midpoint)) %>% as_tibble() %>% setNames(c("lon", "lat"))

coords <- cbind(as.data.frame(df$COMID), coords)
coordinates(coords) <- ~lon + lat
proj4string(coords) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
coords <- spTransform(coords, projection(nlcd_raster))

# landcover <- data.frame()

# for (q in 1:nrow(coords)) {
#   landcover[q] <- extract(nlcd_raster, coords[q], buffer = buffer_dist)
# }

landcover <- raster::extract(nlcd_raster, coords, buffer = buffer_dist, na.rm = F)

landcover_proportions <- lapply(landcover, function(x) {
  counts_x <- table(x)
  proportions_x <- prop.table(counts_x)
  # sort(proportions_x)
})
# l <- unlist(landcover_proportions)

# lc <- rbindlist(landcover_proportions, fill = TRUE)
# 

##### START HERE!!! ######

# dfs <- lapply(landcover_proportions, data.frame, stringsAsFactors = F)
# 
# full_df <- do.call("rbind.fill", dfs)
# 
# 
# test <- bind_rows(lapply(dfs, bind_rows))
# 
# lcs <- bind_rows(dfs)
# 
# lcs <- unnest_wider(as.data.frame(landcover_proportions))

max.length <- max(sapply(landcover_proportions, length))
l <- lapply(landcover_proportions, function(v) {c(v, rep(NA, max.length-length(v)))})
# lc_dfs <- lapply(l, data.frame, stringsAsFactors = F)

lcs <- do.call("rbind", l)





#####################################################################################################

blank_map <- tm_shape(n_net_data_sf) +
  tm_lines() +
  tm_scale_bar()

tmap_save(blank_map, 'figs/FarmingtonHR_lines.png', width = 6, height = 6)


frac_removed_map <- tm_shape(n_net_data_sf) +
  tm_lines(col = 'frac_removed', lwd = 2, title.col = 'N Percent Removed, N load = 0.001') +
  tm_scale_bar()

tmap_save(frac_removed_map, 'figs/FarmingtonHR_frac_removed_load_0.001.png', width = 6, height = 6)

n_exported_map <- tm_shape(n_net_data_sf) +
  tm_lines(col = 'Nout', lwd = 2, title.col = 'N exported (kg/d), N load = 0.001') +
  tm_scale_bar()

tmap_save(n_exported_map, 'figs/FarmingtonHR_N_export_load_0.001.png', width = 6, height = 6)





End_time <- Sys.time()

print(End_time - Start_time)

# mapview::mapview(n_net_data_sf, zcol = 'Nout')