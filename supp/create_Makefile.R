n_reps = c(4,8,16)
n_sims = 10

makefile = 'Makefile'

screen_width = 70

catf = function(...) cat(..., '\n', sep='', file=makefile, append=TRUE)

cat(rep('#', screen_width), '\n', sep='', file=makefile, append=FALSE)
catf('#     DO NOT TOUCH: Automatically created using create_Makefile.R')
catf(rep('#', screen_width),'\n')

catf('RCMD = R CMD BATCH --vanilla\n')

# Main targets
methods = c('laplace','normal')

catf('all: ', paste(methods, collapse=' '))

for (m in methods) {
  catf(m, ': ', paste(m, n_reps, sep='', collapse=' '))
}


catf('\nclean: ', paste('clean-', methods, sep='', collapse=' '))

for (m in methods) {
  catf('clean-', m, ': ', paste('clean-', m, n_reps, sep='', collapse=' '))
}


catf('\nfiles = README.txt script.R figs.R get_hyperparameters.R single_gene_analysis.R model.stan')
catf('zip: $(files); zip supp.zip $(files)\n')


# Sub targets

for (m in methods) {
  for (r in n_reps) {
    catf(rep('#',screen_width))
    catf('# Targets for ', r, ' reps per variety using method ', m)
    catf(rep('#',screen_width))  
    
    # Vectors of files
    result_file = paste0('results/',m,'/results-', r, '-', 1:n_sims,'.rds')
    data_file   = paste0('sim-data/sim-',          r, '-', 1:n_sims,'.rds')
    Rout_file   = paste0('Rout/',m,'/sim-',        r, '-', 1:n_sims,'.Rout')
    
    catf(m, r, ': ', paste(result_file, collapse=' '),'\n')
    
    for (i in 1:n_sims) {
      catf(result_file[i], ": sim-script.R ", data_file[i], "\n\t", 
           "$(RCMD) '--args r=", r, " i=", i,' m=\"', m, "\"' sim-script.R ",
           Rout_file[i], " \n")
    }
    
    catf("clean-", m, r,":\n\trm -fv ", 
         paste(result_file, Rout_file, collapse=' '),"\n")
  }
}

