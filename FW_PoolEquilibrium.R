n<-10 # number of species pools
B<-rep(1,N)
P<-rep(0.1,N) # per biomass production
M<-rep(0.001,N) # per biomass metabolism
W<-rep(0.0001,N) # per biomass waste of food - could be a matrix giving waste from each prey item
# a    matrix of availability of one taxon to another
# B.D=sum across all nodes[ aB(P+M+W) ]

# N+sum(bBW)=sum across all nodes[ aB(P+M+W) ]  # N = nutrients replenished in a year, b is the proportion of waste recycled
# D for terminal nodes are fixed


solveB<-function(B    # biomass vector to estimate
                ,N    # nutrients (if a vector)
                ,BT=0 # terminal node biomass if needed
                ,p    # per bimoass productivity
                ,m    # per biomass metabolism
                ,w    # per biomass waste
                ,a    # matrix of availabilities (Ecotrophic Efficiencies in Ecopath)
                ,b    # vector of waste recycling
){
  # loop through Nodes not including nutrients
  # calculate N =  sum across all nodes[ aB(P+M+W) ] - sum(bBW) 
}


Include competition

Randomise over uncertainty including starting vector
Do PCA on results to show main forces and orientations.  Need a centroid to display the food web structure
In the function, can include part of starting vector as absolute or relative abundances with error
Can you do a likelihood of a vector based on what the linkages would need to be to make it work (departure of linkages from estimated linkages))


Notes on sea ice
Will Hobbs paper (to come)
Changes in season to season persistence -Libera etal 2022  (should lose memory from one summer to next because of depth of mixed layer
see Schroeter et al 2023 on spatial coherence in sea ice anomalies since abrupt critical transition (2006)                                                            
                                                            )
Alexander Haumann (AWI) abrupt sea ice ocean transition
following on Will Hobbs talk on role of ocean in sea ice transition
declines observed previouslyin spring season but now in autumn (growing) season
abrupt change in complex systems (Scheffer etal 2009)
Importance of subsurface ocean warming (Meehle)


