## Context specific growth of the SARS-CoV-2 Delta variant
 Code and data for all the statistical analysis done for estimating the context specific growth of the SARS-CoV-2 Delta variant across England and US.

### The packages used in the R scripts are: 
dplyr, lubridate, ggplot2, rstan, reshape2, loo, posterior, bayesplot, gridExtra, data.table, cowplot, boot

### How to run the code:
1. Clone the repository to your local machine.<br>
2. Open terminal/bash and change the directory to where the file is cloned, using ```cd```. You should be able to see folders for data, outputs and src in this directory.<br>
3. Type ```make``` and press enter.<br>
4. You can now check the pdf files of the outputs in the ```outputs``` folder.

### Note:
Due to data sharing agreements, the data on relative self mobility has not been made available, and instead a column with random numbers has been used. Hence the results will be different to those reported in the paper.
