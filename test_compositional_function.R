comp_transform <- function (x, transform = "compositional", target = "OTU", shift = 0, 
          scale = 1) 
{
  y <- NULL
  xorig <- x
  if (class(x) == "phyloseq") {
    x <- abundances(x)
  }
  if (!all(sample(round(prod(dim(abundances(x)))/10)))%%1 == 
      0) {
    warning("The OTU abundances are not integers. \n        Check that the OTU input data is given as original counts \n        to avoid transformation errors!")
  }
  if (min(sample_sums(otu_table(ps))) == 0) {
    warning("At least one sample has no reads...cannot divide by zero!")
  }
  if (transform == "compositional") {
    if (target == "OTU") {
      xt <- apply(x, 2, function(x) {
        log(x/sum(x) + (1/min(sample_sums(otu_table(ps)))))
      })
    }
  }
}
  xret <- xt
  if (class(xorig) == "phyloseq") {
    if (taxa_are_rows(xorig)) {
      otu_table(xorig)@.Data <- xret
    }
    else {
      otu_table(xorig)@.Data <- t(xret)
    }
    xret <- xorig
  }
  xret

#<environment: namespace:microbiome>