# Terrestrial Analyses
# based on Jasmine Lee's work and the conservation zones


# from Jasmine
#   High resolution vector polylines of the Antarctic coastline
#         https://data.bas.ac.uk/items/45c3cc90-098b-45e3-a809-16b80eed4ec2/

#           Coastline for Antarctica created from various mapping and remote sensing sources, consisting of the following coast types: ice coastline, rock coastline, grounding line, ice shelf and front, ice rumple, and rock against ice shelf. Covering all land and ice shelves south of 60°S. Suitable for topographic mapping and analysis. High resolution versions of ADD data are suitable for scales larger than 1:1,000,000. The largest suitable scale is changeable and dependent on the region.
#           Changes in v7.9 include updates to the South Orkney Islands, sections of the Filchner, West, Wilkins and George VI ice shelves, and ice fronts in eastern Dronning Maud Land and west of Law Dome.
#           Data compiled, managed and distributed by the Mapping and Geographic Information Centre and the UK Polar Data Centre, British Antarctic Survey on behalf of the Scientific Committee on Antarctic Research.
#           Citation:
#             Gerrish, L., Ireland, L., Fretwell, P., & Cooper, P. (2024). High resolution vector polylines of the Antarctic coastline (7.9) [Data set]. UK Polar Data Centre, Natural Environment Research Council, UK Research & Innovation. https://doi.org/10.5285/45c3cc90-098b-45e3-a809-16b80eed4ec2
#           If using for a graphic or if short on space, please cite as 'data from the SCAR Antarctic Digital Database, 2024'


#     Projections of Antarctic ice-free areas under two RCP climate change scenarios
#        AAS_4297_Future_Ice-free_Layers
#        https://data.aad.gov.au/metadata/AAS_4297_Future_Ice-free_Layers
#        Citation:
#            Lee, J., and Terauds, A. (2017) Projections of Antarctic ice-free areas under two RCP climate change scenarios, Ver. 1, Australian Antarctic Data Centre - doi:10.4225/15/585216f8703d0, Accessed: 2024-06-22
#        Distribution models of ice-free area expansion under various climate change scenarios will be generated using rock layers from the ADD (v6), the BAS bedmap ice-thickness model and various climate layers and projections. The data product will comprise ice-free area distribution and/or melt patterns in the form of pdf maps and/or ESRI shapefiles and GeoTiffs.
#
#        The current ice-free layer, that was used to project the future ice-free layers, is the 'medium resolution' rock outcrop layer from the Antarctic Digital Database (ADD) version 7 (add.scar.org). The projected melt was applied to the Bedmap2 'ice thickness' layer from the British Antarctic Survey (bas.ac.uk/project/bedmap-2/), in order to derive the future ice-free layers.
#
#        This data set conforms to the CCBY Attribution License (http://creativecommons.org/licenses/by/4.0/). Please follow instructions listed in the citation reference provided at http://data.aad.gov.au/aadc/metadata/citation.cfm?entry_id=AAS_4297_Future_Ice-free_Layers when using these data. These data have been used in the following publication: Lee, J.R., Raymond, B., Bracegirdle, T.J., Chades, I., Fuller, R.A., Shaw, J.D. and Terauds, A. 2017. Climate change drives expansion of Antarctic ice-free habitat. Nature

# Bioregions (from Dana)
#    
#      ASAC_3095_Antarctic_Biogeography_shapefiles
#         Citation:
#            Terauds, A. (2012) Conservation Biogeography of the Antarctic - Shapefiles, Ver. 1, Australian Antarctic Data Centre - doi:10.4225/15/54AF512520027, Accessed: 2024-06-22
#            doi:10.4225/15/54AF512520027

#         https://data.aad.gov.au/metadata/records/ASAC_3095_Antarctic_Biogeography_shapefiles
#
#    An update of the Antarctic Conservation Biogeographic Regions (ACBRs)
#         AAS_4296_Antarctic_Conservation_Biogeographic_Regions_v2
#         https://data.aad.gov.au/metadata/AAS_4296_Antarctic_Conservation_Biogeographic_Regions_v2?search=datasets
#         Citation:
#            Terauds, A. (2016) An update of the Antarctic Conservation Biogeographic Regions (ACBRs), Ver. 1, Australian Antarctic Data Centre - doi:10.4225/15/5729930925224, Accessed: 2024-06-22


library(terra)
library(ggplot2)
library(tidyterra)
library(patchwork)

# input shapefiles for developing buffer (SCAR ADD)
sCoast<-vect("/Users/acon/Desktop/_w/_d/Shapefiles/Antarctica ice free/add_coastline_high_res_line_v7_9/add_coastline_high_res_line_v7_9.shp",crs="EPSG:3031")
unique(sCoast$surface)
plot(sCoast[sCoast$surface=="rock coastline"])

SCoast_ice_rock<-sCoast[sCoast$surface=="rock coastline" | sCoast$surface=="ice coastline"]
plot(SCoast_ice_rock)
pCirBuffer<-aggregate(buffer(SCoast_ice_rock, width=5000))  # aggregate overlapping buffers
plot(pCirBuffer)

# input coast and ice polygons for plotting 

sCoast_poly<-vect("/Users/acon/Desktop/_w/_d/Shapefiles/Antarctica\ terrestrial/add_coastline_medium_res_polygon_v7_9.shp/add_coastline_medium_res_polygon_v7_9.shp",crs="EPSG:3031")
unique(sCoast_poly$surface)  # "land"       "ice shelf"  "ice tongue" "rumple" 

# input shapefiles of ice-free areas in biogeography
sTerrestrialBiogV2<-vect("/Users/acon/Desktop/_w/_d/Shapefiles/Antarctica ice free/AAS_4296_Antarctic_Conservation_Biogeographic_Regions_v2/ACBRs_v2_2016/ACBRs_v2_2016.shp")
plot(sTerrestrialBiogV2)
tbAP<-c("South Orkney Islands","North-west Antarctic Peninsula"
        ,"North-east Antarctic Peninsula","Central south Antarctic Peninsula"
        ,"South Antarctic Peninsula")
unique(sTerrestrialBiogV2$ACBR_Name)
plot(sTerrestrialBiogV2[!(sTerrestrialBiogV2$ACBR_Name %in% tbAP)])

lines(SCoast_ice_rock,col="red")
pTB<-buffer(sTerrestrialBiogV2,width=0)

tbCoast<-crop(pTB,pCirBuffer) # keep areas within the coastal buffer
# crop to each MEASO sector and sum the area
# use MEASO shapes extended inland to -80S and separating Antarctic Peninsula

MEASOinland<-vect("/Users/acon/Desktop/_w/_d/Shapefiles/MEASO\ polygons\ Inland/MEASO\ areas\ inland\ polar.shp")
MEASOinland<-project(MEASOinland,"EPSG:3031")
MEASOinland
plot(MEASOinland)
lines(SCoast_ice_rock,col="red")

plot(crop(project(tbCoast,"EPSG:4326"),ext(-65,-55,-66,-61)))
lines(project(MEASOinland,"EPSG:4326"),col="red")

MEASOareas<-c("AOA","CIA","EIA","WPA","EPA") # note there are 2 polygons for "WPA"
pM<-MEASOinland[MEASOinland$name %in% MEASOareas]
plot(pM)

RockArea<-do.call(c,lapply(MEASOareas,function(a,pM,tbC){
  v<-crop(tbC,pM[pM$name %in% a])
  sum(expanse(v))
  },MEASOinland,tbCoast))*1E-6  # km2
RockArea<-data.frame(Period=rep("Current",length(MEASOareas)),MEASO=MEASOareas,Area=RockArea)
RockArea


# change according to RCP8.5 (use best estimate rather than bounds)
# use same buffers to crop shapes

sFutureIceFree<-list(RCP45 = vect("/Users/acon/Desktop/_w/_d/Shapefiles/Antarctica ice free/AAS_4297_Future_Ice-free_Layers/AAS_4297_Ice_Free_Shapefiles/PS_RCP45_Best_Future_IceFree.shp")
                     ,RCP85 = vect("/Users/acon/Desktop/_w/_d/Shapefiles/Antarctica ice free/AAS_4297_Future_Ice-free_Layers/AAS_4297_Ice_Free_Shapefiles/PS_RCP85_Best_Future_IceFree.shp"))

plot(sFutureIceFree$RCP85)
#pRCP<-buffer(sFutureIceFree$RCP85,width=0)

# plot(crop(pCirBuffer,ext(21E5,22E5,12E5,13E5)))

rcp585Coast<-crop(pRCP,pCirBuffer) # keep areas within the coastal buffer
plot(rcp585Coast)

rcp585Area<-do.call(c,lapply(MEASOareas,function(a,pM,rcpC){
  v<-crop(rcpC,pM[pM$name %in% a])
  sum(expanse(v))
},MEASOinland,rcp585Coast))*1E-6  # km2
rcp585Area<-data.frame(Period=rep("2100",length(MEASOareas)),MEASO=MEASOareas,Area=rcp585Area)
rcp585Area

CoastalAreas<-rbind(RockArea,rcp585Area)
CoastalAreas$Period<-factor(CoastalAreas$Period, levels=c("Current","2100"))


# now plot biogeographic map with circles showing areas now and expected in future
sMEASOcoast<-vect("/Users/acon/Desktop/_w/_d/Shapefiles/MEASO\ polygons\ coast/MEASO_polygons.shp")
#sMEASOcoast<-sMEASOcoast[sMEASOcoast$MEASO_area %in% MEASOareas]
sMEASOcoast<-project(sMEASOcoast,"EPSG:3031")

MEASOsectorRGBA<-matrix(c(1,    249,188,188, 1  # alpha level does not work at present
                         ,4,    249,215,194, 1
                         ,7,    209,211,233, 1
                         ,10,   160,197,199, 1
                         ,13,   225,184,216, 1),ncol=5,byrow=TRUE)
dimnames(MEASOsectorRGBA)[[2]]<-c("Area","R","G","B", "Alpha")
mc<-MEASOsectorRGBA
MEASOsectorCols<-as.vector(rgb(red=mc[,"R"],green=mc[,"G"],blue=mc[,"B"],maxColorValue = 255))
MEASOcols<-as.vector(matrix(rep(MEASOsectorCols,3),nrow=3,byrow=TRUE))
names(MEASOcols)<-c("AOA", "AOS", "AON", "CIA", "CIS", "CIN", "EIA", "EIS", "EIN", "EPA", "EPS", "EPN", "WPA", "WPS", "WPN")


sMEASOcoast_crop<-crop(sMEASOcoast,ext(-2.9E6,3E6,-2.6E6,2.6E6))

p<-ggplot(sMEASOcoast_crop)+geom_spatvector(aes(fill=MEASO_area,colour=MEASO_area))+scale_fill_manual(values=MEASOcols)+scale_colour_manual(values=MEASOcols)
p<- p+ guides(fill=FALSE,colour=FALSE)

p<-ggplot(sMEASOcoast_crop)
p<-p+geom_spatvector(data=sCoast_poly[sCoast_poly$surface=="land"],fill="white",colour="black",size=0.1)
p<-p+geom_spatvector(data=sCoast_poly[sCoast_poly$surface!="land"],fill="lightblue",colour="black",size=0.1)
p<-p+geom_spatvector(data=sTerrestrialBiogV2[!(sTerrestrialBiogV2$ACBR_Name %in% tbAP)],fill="chartreuse1",colour="chartreuse1")
p<-p+geom_spatvector(data=sTerrestrialBiogV2[(sTerrestrialBiogV2$ACBR_Name %in% tbAP)],fill="darkgreen",colour="darkgreen")
p <- p+ theme_void()
pTB<-p

#################################
# Bubble plots on maps

# bubble locations

ca<-data.frame(MEASO = c("AOA", "CIA", "EIA", "WPA", "EPA")
             , x =     c(0.20,  0.95,  0.90,  0.45,  0.03)
             , y =     c(0.95,  0.79,  0.40,  0.38,  0.7)
             )
ca<-merge(ca,CoastalAreas,by="MEASO")

ca<-cbind(ca,Fill=ca$MEASO)
ca$Fill[ca$Period=="2100"]<-"transparent"
ca<-ca[order(ca$Period,decreasing=TRUE),]
ca

MEASOcols<-c(MEASOcols,NA)
MEASOcols[16]<-"grey90"
names(MEASOcols)<-c("AOA", "AOS", "AON", "CIA", "CIS", "CIN", "EIA", "EIS", "EIN", "EPA", "EPS", "EPN", "WPA", "WPS", "WPN","transparent")

#
# do bubble plot
p<-ggplot(ca,aes(x=x,y=y,size=Area, fill=Fill))
p<-p+ geom_point(alpha=1,shape=21,colour="black")  # alpha is opacity
p<-p+geom_text(aes(label=MEASO),size=2.8)
p<-p+ scale_size_area(name="Area"  #,guide="legend" ,position="bottom" #,range  = c(0, 0.7)
                      ,max_size=30, breaks = c(500, 1000, 5000, 15000))
p<-p+  scale_fill_manual(values=MEASOcols,aesthetics="fill")

p<-p+ coord_fixed(ratio=1) + xlim(0,1) + ylim(0,1) + xlab("")+ylab("")+
  theme(
    axis.text.x = element_blank(),axis.text.y = element_blank()
    ,axis.ticks = element_blank()
    ,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ,panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA) #transparent plot bg
    ,legend.position = "bottom"
    # ,legend.position.inside = c(0.8, 0.1)
    ,legend.background = element_rect(fill='transparent', color=NA)
  ) # end theme
p<-p+guides(size=guide_legend(position="top",title=bquote(Coastal~ice-free~area~km^{2}))
            ,fill=FALSE)
p


g<-pTB
g<-g+inset_element( p
                    ,left   = -0.2 
                    ,bottom = -0.3
                    ,right  = 1.2 
                    ,top    = 1.1
                    ,on_top = TRUE)
g



