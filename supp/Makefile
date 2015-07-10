all: sim-results1.rds \
     sim-results30.rds 

sim-results1.rds: 
  R CMD BATCH --vanilla '--args i=1' sim-script.R 
sim-results30.rds: 
  R CMD BATCH --vanilla '--args i=30' sim-script.R 
