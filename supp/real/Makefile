RCMD = R CMD BATCH --vanilla

all: results/ji.rds results/normal.rds results/laplace.rds \
     results/normal-parallel.rds results/laplace-parallel.rds

results/%.rds: %.R
	$(RCMD) $< 

