# mpralm adapted for allowing individual barcode input

This is a fork from [the original mpralm package](https://github.com/hansenlab/mpra/tree/master), an R package that provides tools for differential analysis in MPRA studies.  
The code has been modified to allow for using individual barcodes as input to the model.  
In order to use individual barcodes, call the ``mpralm`` function with ``aggregate = "none"`` and indicate in the passed ``block`` argument which replicate each barcode belongs to. 
