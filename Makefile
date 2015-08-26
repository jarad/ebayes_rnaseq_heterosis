depends = \
  ebayes_rnaseq_heterosis.pdf \
  response.pdf \
  ebayes_rnaseq_heterosis.tex \
  response.tex \
  jarad.bib biom.cls biom.bst \
  include/estimates.pdf \
  include/volcano.pdf \
  include/exampleROC0_1.pdf \
  include/auc-facet-TRUE.pdf \
  include/gene_specific_estimates.pdf

zip: $(depends); zip ebayes_rnaseq_heterosis.zip $(depends)
