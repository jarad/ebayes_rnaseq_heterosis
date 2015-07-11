n_reps = c(4,8,16)
n_sims = 10

makefile = 'Makefile'

r_command = 'R CMD BATCH --vanilla'

screen_width = 60

cat(rep('#', screen_width), '\n', sep='', file=makefile, append=FALSE)
cat('#     DO NOT TOUCH \n', file=makefile, append=TRUE)
cat(rep('#', screen_width), '\n', sep='', file=makefile, append=TRUE)
cat('# Automatically created using make_makefile.R\n\n', 
    file=makefile, 
    append=TRUE)

for (r in n_reps) {
  cat(rep("#",screen_width), "\n", sep='', file=makefile, append=TRUE)
  cat("# Targets for ", r, "reps per variety\n", file=makefile, append=TRUE)
  cat(rep("#",screen_width), "\n", sep='', file=makefile, append=TRUE)  
  
  # Vectors of files
  result_file = paste0("results/results-", r, "-", 1:n_sims,'.rds')
  data_file   = paste0("sim-data/sim-", r, "-", 1:n_sims,'.rds')
  Rout_file   = paste0("sim-", r, "-", 1:n_sims,'.Rout')
  
  cat("make",r, ': ', paste(result_file, collapse=' '),"\n\n" , 
      sep = '', 
    file=makefile, 
    append=TRUE)
  
  for (i in 1:n_sims) {
    cat(result_file[i], ": sim-script.R ", data_file[i], "\n\t", 
        r_command, " '--args r=", r, " i=", i, "' sim-script.R ",
        Rout_file[i], " \n\n",
        sep='',
        file=makefile,
        append=TRUE)
  }
  
  cat("clean",r,":\n\trm -fv ", 
      paste(result_file, collapse=' '),"\n\n",
      sep='',
    file=makefile, 
    append=TRUE)
}

