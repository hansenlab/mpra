# mpralm adapted for allowing individual barcode input

This is a fork from [the original mpralm package](https://github.com/hansenlab/mpra/tree/master), an R package that provides tools for differential analysis in MPRA studies.  
The code has been modified to allow for using individual barcodes as input to the model.  
In order to use individual barcodes, call the ``mpralm`` function with ``aggregate = "none"`` and indicate in the passed ``block`` argument which replicate each barcode belongs to.   
Apart from this, the package works exactly the same. For an example on how to run, see the [original vignette](https://rdrr.io/bioc/mpra/f/vignettes/mpra.Rmd).

Original publication:  
Myint, Leslie, Dimitrios G. Avramopoulos, Loyal A. Goff, and Kasper D. Hansen. *Linear models enable powerful differential activity analysis in massively parallel reporter assays.* BMC Genomics 2019, 209. doi: 10.1186/s12864-019-5556-x.
