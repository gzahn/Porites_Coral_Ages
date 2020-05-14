comp_transform <- function (x, transform = "compositional", target = "OTU", shift = 0, 
    scale = 1) 
{
    y <- NULL
    xorig <- x
    if (class(x) == "phyloseq") {
        x <- abundances(x)
    }
    if (transform == "relative.abundance") {
        transform <- "compositional"
    }
    if (min(sample_sums(otu_table(ps))) == 0) {
    warning("At least one sample has no reads...cannot divide by zero!")
    }
    if (!all(sample(round(prod(dim(abundances(x)))/10)))%%1 == 
        0) {
        warning("The OTU abundances are not integers. \n        Check that the OTU input data is given as original counts \n        to avoid transformation errors!")
    }
    if (transform == "compositional") {
        if (target == "OTU") {
            xt <- apply(x, 2, function(x) {
                abs(log(x/sum(x) + (1/min(sample_sums(otu_table(ps))))))
            })
        }
        else if (target == "sample") {
            xt <- apply(x, 1, function(x) {
                x/max(sum(x), 1e-32)
            })
        }
    }
    else if (transform == "Z") {
        xt <- ztransform(x, target)
    }
    else if (transform == "clr") {
        if (any(abundances(x) < 0)) {
            stop("Non-negative data matrix required for the \n            clr transformation. Exiting.")
        }
        xt <- x
        xt <- transform(xt, "compositional")
        colnames(xt) <- colnames(x)
        if (any(xt == 0)) {
            v <- as.vector(xt)
            minval <- min(v[v > 0])/2
            xt <- xt + minval
        }
        d <- t(apply(xt, 2, function(x) {
            log(x) - mean(log(x))
        }))
        if (nrow(d) == ncol(xt)) {
            rownames(d) <- colnames(xt)
            colnames(d) <- rownames(xt)
        }
        else {
            colnames(d) <- colnames(xt)
            rownames(d) <- rownames(xt)
        }
        xt <- t(d)
    }
    else if (transform == "log10") {
        if (min(x) == 0) {
            warning("OTU table contains zeroes. Using log10(1 + x) transform.")
            xt <- log10(1 + x)
        }
        else {
            xt <- log10(x)
        }
    }
    else if (transform == "log10p") {
        xt <- log10(1 + x)
    }
    else if (transform == "identity") {
        xt <- x
    }
    else if (transform == "shift") {
        xt <- x + shift
    }
    else if (transform == "scale") {
        xt <- scale * x
    }
    else {
        a <- try(xt <- decostand(x, method = transform, MARGIN = 2))
        if (class(a) == "try-error") {
            xt <- NULL
            stop(paste("Transformation", transform, "not defined."))
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
}
