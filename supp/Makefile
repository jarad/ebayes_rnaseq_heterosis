######################################################################
#     DO NOT TOUCH: Automatically created using create_Makefile.R
######################################################################

RCMD = R CMD BATCH --vanilla

all: laplace normal
laplace: laplace4 laplace8 laplace16
normal: normal4 normal8 normal16

clean: clean-laplace clean-normal
clean-laplace: clean-laplace4 clean-laplace8 clean-laplace16
clean-normal: clean-normal4 clean-normal8 clean-normal16

files = README.txt script.R figs.R get_hyperparameters.R single_gene_analysis.R model.stan
zip: $(files); zip supp.zip $(files)

######################################################################
# Targets for 4 reps per variety using method laplace
######################################################################
laplace4: results/laplace/results-4-1.rds results/laplace/results-4-2.rds results/laplace/results-4-3.rds results/laplace/results-4-4.rds results/laplace/results-4-5.rds results/laplace/results-4-6.rds results/laplace/results-4-7.rds results/laplace/results-4-8.rds results/laplace/results-4-9.rds results/laplace/results-4-10.rds

results/laplace/results-4-1.rds: sim-script.R sim-data/sim-4-1.rds
	$(RCMD) '--args r=4 i=1 m="laplace"' sim-script.R Rout/laplace/sim-4-1.Rout 

results/laplace/results-4-2.rds: sim-script.R sim-data/sim-4-2.rds
	$(RCMD) '--args r=4 i=2 m="laplace"' sim-script.R Rout/laplace/sim-4-2.Rout 

results/laplace/results-4-3.rds: sim-script.R sim-data/sim-4-3.rds
	$(RCMD) '--args r=4 i=3 m="laplace"' sim-script.R Rout/laplace/sim-4-3.Rout 

results/laplace/results-4-4.rds: sim-script.R sim-data/sim-4-4.rds
	$(RCMD) '--args r=4 i=4 m="laplace"' sim-script.R Rout/laplace/sim-4-4.Rout 

results/laplace/results-4-5.rds: sim-script.R sim-data/sim-4-5.rds
	$(RCMD) '--args r=4 i=5 m="laplace"' sim-script.R Rout/laplace/sim-4-5.Rout 

results/laplace/results-4-6.rds: sim-script.R sim-data/sim-4-6.rds
	$(RCMD) '--args r=4 i=6 m="laplace"' sim-script.R Rout/laplace/sim-4-6.Rout 

results/laplace/results-4-7.rds: sim-script.R sim-data/sim-4-7.rds
	$(RCMD) '--args r=4 i=7 m="laplace"' sim-script.R Rout/laplace/sim-4-7.Rout 

results/laplace/results-4-8.rds: sim-script.R sim-data/sim-4-8.rds
	$(RCMD) '--args r=4 i=8 m="laplace"' sim-script.R Rout/laplace/sim-4-8.Rout 

results/laplace/results-4-9.rds: sim-script.R sim-data/sim-4-9.rds
	$(RCMD) '--args r=4 i=9 m="laplace"' sim-script.R Rout/laplace/sim-4-9.Rout 

results/laplace/results-4-10.rds: sim-script.R sim-data/sim-4-10.rds
	$(RCMD) '--args r=4 i=10 m="laplace"' sim-script.R Rout/laplace/sim-4-10.Rout 

clean-laplace4:
	rm -fv results/laplace/results-4-1.rds Rout/laplace/sim-4-1.Rout results/laplace/results-4-2.rds Rout/laplace/sim-4-2.Rout results/laplace/results-4-3.rds Rout/laplace/sim-4-3.Rout results/laplace/results-4-4.rds Rout/laplace/sim-4-4.Rout results/laplace/results-4-5.rds Rout/laplace/sim-4-5.Rout results/laplace/results-4-6.rds Rout/laplace/sim-4-6.Rout results/laplace/results-4-7.rds Rout/laplace/sim-4-7.Rout results/laplace/results-4-8.rds Rout/laplace/sim-4-8.Rout results/laplace/results-4-9.rds Rout/laplace/sim-4-9.Rout results/laplace/results-4-10.rds Rout/laplace/sim-4-10.Rout

######################################################################
# Targets for 8 reps per variety using method laplace
######################################################################
laplace8: results/laplace/results-8-1.rds results/laplace/results-8-2.rds results/laplace/results-8-3.rds results/laplace/results-8-4.rds results/laplace/results-8-5.rds results/laplace/results-8-6.rds results/laplace/results-8-7.rds results/laplace/results-8-8.rds results/laplace/results-8-9.rds results/laplace/results-8-10.rds

results/laplace/results-8-1.rds: sim-script.R sim-data/sim-8-1.rds
	$(RCMD) '--args r=8 i=1 m="laplace"' sim-script.R Rout/laplace/sim-8-1.Rout 

results/laplace/results-8-2.rds: sim-script.R sim-data/sim-8-2.rds
	$(RCMD) '--args r=8 i=2 m="laplace"' sim-script.R Rout/laplace/sim-8-2.Rout 

results/laplace/results-8-3.rds: sim-script.R sim-data/sim-8-3.rds
	$(RCMD) '--args r=8 i=3 m="laplace"' sim-script.R Rout/laplace/sim-8-3.Rout 

results/laplace/results-8-4.rds: sim-script.R sim-data/sim-8-4.rds
	$(RCMD) '--args r=8 i=4 m="laplace"' sim-script.R Rout/laplace/sim-8-4.Rout 

results/laplace/results-8-5.rds: sim-script.R sim-data/sim-8-5.rds
	$(RCMD) '--args r=8 i=5 m="laplace"' sim-script.R Rout/laplace/sim-8-5.Rout 

results/laplace/results-8-6.rds: sim-script.R sim-data/sim-8-6.rds
	$(RCMD) '--args r=8 i=6 m="laplace"' sim-script.R Rout/laplace/sim-8-6.Rout 

results/laplace/results-8-7.rds: sim-script.R sim-data/sim-8-7.rds
	$(RCMD) '--args r=8 i=7 m="laplace"' sim-script.R Rout/laplace/sim-8-7.Rout 

results/laplace/results-8-8.rds: sim-script.R sim-data/sim-8-8.rds
	$(RCMD) '--args r=8 i=8 m="laplace"' sim-script.R Rout/laplace/sim-8-8.Rout 

results/laplace/results-8-9.rds: sim-script.R sim-data/sim-8-9.rds
	$(RCMD) '--args r=8 i=9 m="laplace"' sim-script.R Rout/laplace/sim-8-9.Rout 

results/laplace/results-8-10.rds: sim-script.R sim-data/sim-8-10.rds
	$(RCMD) '--args r=8 i=10 m="laplace"' sim-script.R Rout/laplace/sim-8-10.Rout 

clean-laplace8:
	rm -fv results/laplace/results-8-1.rds Rout/laplace/sim-8-1.Rout results/laplace/results-8-2.rds Rout/laplace/sim-8-2.Rout results/laplace/results-8-3.rds Rout/laplace/sim-8-3.Rout results/laplace/results-8-4.rds Rout/laplace/sim-8-4.Rout results/laplace/results-8-5.rds Rout/laplace/sim-8-5.Rout results/laplace/results-8-6.rds Rout/laplace/sim-8-6.Rout results/laplace/results-8-7.rds Rout/laplace/sim-8-7.Rout results/laplace/results-8-8.rds Rout/laplace/sim-8-8.Rout results/laplace/results-8-9.rds Rout/laplace/sim-8-9.Rout results/laplace/results-8-10.rds Rout/laplace/sim-8-10.Rout

######################################################################
# Targets for 16 reps per variety using method laplace
######################################################################
laplace16: results/laplace/results-16-1.rds results/laplace/results-16-2.rds results/laplace/results-16-3.rds results/laplace/results-16-4.rds results/laplace/results-16-5.rds results/laplace/results-16-6.rds results/laplace/results-16-7.rds results/laplace/results-16-8.rds results/laplace/results-16-9.rds results/laplace/results-16-10.rds

results/laplace/results-16-1.rds: sim-script.R sim-data/sim-16-1.rds
	$(RCMD) '--args r=16 i=1 m="laplace"' sim-script.R Rout/laplace/sim-16-1.Rout 

results/laplace/results-16-2.rds: sim-script.R sim-data/sim-16-2.rds
	$(RCMD) '--args r=16 i=2 m="laplace"' sim-script.R Rout/laplace/sim-16-2.Rout 

results/laplace/results-16-3.rds: sim-script.R sim-data/sim-16-3.rds
	$(RCMD) '--args r=16 i=3 m="laplace"' sim-script.R Rout/laplace/sim-16-3.Rout 

results/laplace/results-16-4.rds: sim-script.R sim-data/sim-16-4.rds
	$(RCMD) '--args r=16 i=4 m="laplace"' sim-script.R Rout/laplace/sim-16-4.Rout 

results/laplace/results-16-5.rds: sim-script.R sim-data/sim-16-5.rds
	$(RCMD) '--args r=16 i=5 m="laplace"' sim-script.R Rout/laplace/sim-16-5.Rout 

results/laplace/results-16-6.rds: sim-script.R sim-data/sim-16-6.rds
	$(RCMD) '--args r=16 i=6 m="laplace"' sim-script.R Rout/laplace/sim-16-6.Rout 

results/laplace/results-16-7.rds: sim-script.R sim-data/sim-16-7.rds
	$(RCMD) '--args r=16 i=7 m="laplace"' sim-script.R Rout/laplace/sim-16-7.Rout 

results/laplace/results-16-8.rds: sim-script.R sim-data/sim-16-8.rds
	$(RCMD) '--args r=16 i=8 m="laplace"' sim-script.R Rout/laplace/sim-16-8.Rout 

results/laplace/results-16-9.rds: sim-script.R sim-data/sim-16-9.rds
	$(RCMD) '--args r=16 i=9 m="laplace"' sim-script.R Rout/laplace/sim-16-9.Rout 

results/laplace/results-16-10.rds: sim-script.R sim-data/sim-16-10.rds
	$(RCMD) '--args r=16 i=10 m="laplace"' sim-script.R Rout/laplace/sim-16-10.Rout 

clean-laplace16:
	rm -fv results/laplace/results-16-1.rds Rout/laplace/sim-16-1.Rout results/laplace/results-16-2.rds Rout/laplace/sim-16-2.Rout results/laplace/results-16-3.rds Rout/laplace/sim-16-3.Rout results/laplace/results-16-4.rds Rout/laplace/sim-16-4.Rout results/laplace/results-16-5.rds Rout/laplace/sim-16-5.Rout results/laplace/results-16-6.rds Rout/laplace/sim-16-6.Rout results/laplace/results-16-7.rds Rout/laplace/sim-16-7.Rout results/laplace/results-16-8.rds Rout/laplace/sim-16-8.Rout results/laplace/results-16-9.rds Rout/laplace/sim-16-9.Rout results/laplace/results-16-10.rds Rout/laplace/sim-16-10.Rout

######################################################################
# Targets for 4 reps per variety using method normal
######################################################################
normal4: results/normal/results-4-1.rds results/normal/results-4-2.rds results/normal/results-4-3.rds results/normal/results-4-4.rds results/normal/results-4-5.rds results/normal/results-4-6.rds results/normal/results-4-7.rds results/normal/results-4-8.rds results/normal/results-4-9.rds results/normal/results-4-10.rds

results/normal/results-4-1.rds: sim-script.R sim-data/sim-4-1.rds
	$(RCMD) '--args r=4 i=1 m="normal"' sim-script.R Rout/normal/sim-4-1.Rout 

results/normal/results-4-2.rds: sim-script.R sim-data/sim-4-2.rds
	$(RCMD) '--args r=4 i=2 m="normal"' sim-script.R Rout/normal/sim-4-2.Rout 

results/normal/results-4-3.rds: sim-script.R sim-data/sim-4-3.rds
	$(RCMD) '--args r=4 i=3 m="normal"' sim-script.R Rout/normal/sim-4-3.Rout 

results/normal/results-4-4.rds: sim-script.R sim-data/sim-4-4.rds
	$(RCMD) '--args r=4 i=4 m="normal"' sim-script.R Rout/normal/sim-4-4.Rout 

results/normal/results-4-5.rds: sim-script.R sim-data/sim-4-5.rds
	$(RCMD) '--args r=4 i=5 m="normal"' sim-script.R Rout/normal/sim-4-5.Rout 

results/normal/results-4-6.rds: sim-script.R sim-data/sim-4-6.rds
	$(RCMD) '--args r=4 i=6 m="normal"' sim-script.R Rout/normal/sim-4-6.Rout 

results/normal/results-4-7.rds: sim-script.R sim-data/sim-4-7.rds
	$(RCMD) '--args r=4 i=7 m="normal"' sim-script.R Rout/normal/sim-4-7.Rout 

results/normal/results-4-8.rds: sim-script.R sim-data/sim-4-8.rds
	$(RCMD) '--args r=4 i=8 m="normal"' sim-script.R Rout/normal/sim-4-8.Rout 

results/normal/results-4-9.rds: sim-script.R sim-data/sim-4-9.rds
	$(RCMD) '--args r=4 i=9 m="normal"' sim-script.R Rout/normal/sim-4-9.Rout 

results/normal/results-4-10.rds: sim-script.R sim-data/sim-4-10.rds
	$(RCMD) '--args r=4 i=10 m="normal"' sim-script.R Rout/normal/sim-4-10.Rout 

clean-normal4:
	rm -fv results/normal/results-4-1.rds Rout/normal/sim-4-1.Rout results/normal/results-4-2.rds Rout/normal/sim-4-2.Rout results/normal/results-4-3.rds Rout/normal/sim-4-3.Rout results/normal/results-4-4.rds Rout/normal/sim-4-4.Rout results/normal/results-4-5.rds Rout/normal/sim-4-5.Rout results/normal/results-4-6.rds Rout/normal/sim-4-6.Rout results/normal/results-4-7.rds Rout/normal/sim-4-7.Rout results/normal/results-4-8.rds Rout/normal/sim-4-8.Rout results/normal/results-4-9.rds Rout/normal/sim-4-9.Rout results/normal/results-4-10.rds Rout/normal/sim-4-10.Rout

######################################################################
# Targets for 8 reps per variety using method normal
######################################################################
normal8: results/normal/results-8-1.rds results/normal/results-8-2.rds results/normal/results-8-3.rds results/normal/results-8-4.rds results/normal/results-8-5.rds results/normal/results-8-6.rds results/normal/results-8-7.rds results/normal/results-8-8.rds results/normal/results-8-9.rds results/normal/results-8-10.rds

results/normal/results-8-1.rds: sim-script.R sim-data/sim-8-1.rds
	$(RCMD) '--args r=8 i=1 m="normal"' sim-script.R Rout/normal/sim-8-1.Rout 

results/normal/results-8-2.rds: sim-script.R sim-data/sim-8-2.rds
	$(RCMD) '--args r=8 i=2 m="normal"' sim-script.R Rout/normal/sim-8-2.Rout 

results/normal/results-8-3.rds: sim-script.R sim-data/sim-8-3.rds
	$(RCMD) '--args r=8 i=3 m="normal"' sim-script.R Rout/normal/sim-8-3.Rout 

results/normal/results-8-4.rds: sim-script.R sim-data/sim-8-4.rds
	$(RCMD) '--args r=8 i=4 m="normal"' sim-script.R Rout/normal/sim-8-4.Rout 

results/normal/results-8-5.rds: sim-script.R sim-data/sim-8-5.rds
	$(RCMD) '--args r=8 i=5 m="normal"' sim-script.R Rout/normal/sim-8-5.Rout 

results/normal/results-8-6.rds: sim-script.R sim-data/sim-8-6.rds
	$(RCMD) '--args r=8 i=6 m="normal"' sim-script.R Rout/normal/sim-8-6.Rout 

results/normal/results-8-7.rds: sim-script.R sim-data/sim-8-7.rds
	$(RCMD) '--args r=8 i=7 m="normal"' sim-script.R Rout/normal/sim-8-7.Rout 

results/normal/results-8-8.rds: sim-script.R sim-data/sim-8-8.rds
	$(RCMD) '--args r=8 i=8 m="normal"' sim-script.R Rout/normal/sim-8-8.Rout 

results/normal/results-8-9.rds: sim-script.R sim-data/sim-8-9.rds
	$(RCMD) '--args r=8 i=9 m="normal"' sim-script.R Rout/normal/sim-8-9.Rout 

results/normal/results-8-10.rds: sim-script.R sim-data/sim-8-10.rds
	$(RCMD) '--args r=8 i=10 m="normal"' sim-script.R Rout/normal/sim-8-10.Rout 

clean-normal8:
	rm -fv results/normal/results-8-1.rds Rout/normal/sim-8-1.Rout results/normal/results-8-2.rds Rout/normal/sim-8-2.Rout results/normal/results-8-3.rds Rout/normal/sim-8-3.Rout results/normal/results-8-4.rds Rout/normal/sim-8-4.Rout results/normal/results-8-5.rds Rout/normal/sim-8-5.Rout results/normal/results-8-6.rds Rout/normal/sim-8-6.Rout results/normal/results-8-7.rds Rout/normal/sim-8-7.Rout results/normal/results-8-8.rds Rout/normal/sim-8-8.Rout results/normal/results-8-9.rds Rout/normal/sim-8-9.Rout results/normal/results-8-10.rds Rout/normal/sim-8-10.Rout

######################################################################
# Targets for 16 reps per variety using method normal
######################################################################
normal16: results/normal/results-16-1.rds results/normal/results-16-2.rds results/normal/results-16-3.rds results/normal/results-16-4.rds results/normal/results-16-5.rds results/normal/results-16-6.rds results/normal/results-16-7.rds results/normal/results-16-8.rds results/normal/results-16-9.rds results/normal/results-16-10.rds

results/normal/results-16-1.rds: sim-script.R sim-data/sim-16-1.rds
	$(RCMD) '--args r=16 i=1 m="normal"' sim-script.R Rout/normal/sim-16-1.Rout 

results/normal/results-16-2.rds: sim-script.R sim-data/sim-16-2.rds
	$(RCMD) '--args r=16 i=2 m="normal"' sim-script.R Rout/normal/sim-16-2.Rout 

results/normal/results-16-3.rds: sim-script.R sim-data/sim-16-3.rds
	$(RCMD) '--args r=16 i=3 m="normal"' sim-script.R Rout/normal/sim-16-3.Rout 

results/normal/results-16-4.rds: sim-script.R sim-data/sim-16-4.rds
	$(RCMD) '--args r=16 i=4 m="normal"' sim-script.R Rout/normal/sim-16-4.Rout 

results/normal/results-16-5.rds: sim-script.R sim-data/sim-16-5.rds
	$(RCMD) '--args r=16 i=5 m="normal"' sim-script.R Rout/normal/sim-16-5.Rout 

results/normal/results-16-6.rds: sim-script.R sim-data/sim-16-6.rds
	$(RCMD) '--args r=16 i=6 m="normal"' sim-script.R Rout/normal/sim-16-6.Rout 

results/normal/results-16-7.rds: sim-script.R sim-data/sim-16-7.rds
	$(RCMD) '--args r=16 i=7 m="normal"' sim-script.R Rout/normal/sim-16-7.Rout 

results/normal/results-16-8.rds: sim-script.R sim-data/sim-16-8.rds
	$(RCMD) '--args r=16 i=8 m="normal"' sim-script.R Rout/normal/sim-16-8.Rout 

results/normal/results-16-9.rds: sim-script.R sim-data/sim-16-9.rds
	$(RCMD) '--args r=16 i=9 m="normal"' sim-script.R Rout/normal/sim-16-9.Rout 

results/normal/results-16-10.rds: sim-script.R sim-data/sim-16-10.rds
	$(RCMD) '--args r=16 i=10 m="normal"' sim-script.R Rout/normal/sim-16-10.Rout 

clean-normal16:
	rm -fv results/normal/results-16-1.rds Rout/normal/sim-16-1.Rout results/normal/results-16-2.rds Rout/normal/sim-16-2.Rout results/normal/results-16-3.rds Rout/normal/sim-16-3.Rout results/normal/results-16-4.rds Rout/normal/sim-16-4.Rout results/normal/results-16-5.rds Rout/normal/sim-16-5.Rout results/normal/results-16-6.rds Rout/normal/sim-16-6.Rout results/normal/results-16-7.rds Rout/normal/sim-16-7.Rout results/normal/results-16-8.rds Rout/normal/sim-16-8.Rout results/normal/results-16-9.rds Rout/normal/sim-16-9.Rout results/normal/results-16-10.rds Rout/normal/sim-16-10.Rout

