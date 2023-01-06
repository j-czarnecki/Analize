# Dispersion relation anlysis

This program changes numeration of dispersion relation bands from lowest-to-highest order to real, physical numeration, consistent with quantum numbers.
It is used to analise dispersion relations and other quantities derived from them, obtained from a self-consistent kp calcuation of semiconducting nanowire. 

**Find_bands.py** contains procedures taking numpy arrays as input and numerating bands in a proper way. 
*  *find_bands()* finds pairs of bands split by zeeman splitting at k=0 
*  *numerate_bands()* changes found pairs of Zeeman-split bands to band numbers
*  *crossing()* uses polynomial extrapolation to predict band crossing and change numeration when bands cross

**Crossing_bands.ipynb** is a demostration of program workflow and output, comparing lowest-to-highest numeration and the one consistent with quantum numbers on two plots.
  
