# Need to pass as arguments (or assign)
# r: number of replicates per variety (4, 8, 16)
# i: simulation number (1-10)
# m: statistical method ('laplace' or 'normal')

# R CMD BATCH --vanilla '--args r=4 i=1 m ="laplace"' sim-script.R 
args=(commandArgs(TRUE))

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  i = 1
}else{
  for(ii in 1:length(args)){
    eval(parse(text=args[[ii]]))
  }
}
