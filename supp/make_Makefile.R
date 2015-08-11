n_reps = c(4,8,16)
n_sims = 10

makefile = 'Makefile'

r_command = 'R CMD BATCH --vanilla'

screen_width = 60

catf = function(...) cat(..., sep='', file=makefile, append=TRUE)

cat(rep('#', screen_width), '\n', file=makefile, append=FALSE)
catf('#     DO NOT TOUCH \n')
catf(rep('#', screen_width), '\n')
catf('# Automatically created using make_Makefile.R\n\n')

# Main targets
methods = c('laplace','normal')

catf('all: ', paste(methods, collapse=' '), '\n\n')

catf('laplace: ', paste('laplace',n_reps, sep='', collapse=' '), '\n\n')

catf('files = README.txt script.R figs.R get_hyperparameters.R single_gene_analysis.R model.stan\n\nzip: $(files); zip supp.zip $(files)\n\n')


# Sub targets

for (m in methods) {
  for (r in n_reps) {
    catf(rep('#',screen_width), '\n')
    catf('# Targets for ', r, 'reps per variety\n')
    catf(rep('#',screen_width), '\n')  
    
    # Vectors of files
    result_file = paste0('results/',m,'/results-', r, '-', 1:n_sims,'.rds')
    data_file   = paste0('sim-data/sim-',          r, '-', 1:n_sims,'.rds')
    Rout_file   = paste0('Rout/',m,'/sim-',        r, '-', 1:n_sims,'.Rout')
    
    catf(m,r, ': ', paste(result_file, collapse=' '),'\n\n')
    
    for (i in 1:n_sims) {
      catf(result_file[i], ": sim-script.R ", data_file[i], "\n\t", 
           r_command, " '--args r=", r, " i=", i, "' sim-script.R ",
           Rout_file[i], " \n\n")
    }
    
    catf("clean",r,":\n\trm -fv ", 
         paste(result_file, collapse=' '),"\n\n")
  }
}

