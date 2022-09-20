## Context specific growth of the SARS-CoV-2 Delta variant
 Code and data for all the statistical analysis done for estimating the context specific growth of the SARS-CoV-2 Delta variant across England.

### The packages used in the R scripts are: 
dplyr, lubridate, ggplot2, rstan, reshape2, loo, posterior, bayesplot, gridExtra, data.table, cowplot, boot

### How to run the code:
1. Clone the repository to your local machine.<br>
2. Open terminal/bash and change the directory to where the file is cloned, using ```cd```. You should be able to see folders for data, outputs and src in this directory.<br>
3. Type ```make``` and press enter.<br>
4. You can now check the pdf files of the outputs in the ```outputs``` folder.

### Note:
Due to data sharing agreements, the data on relative self mobility has not been made available, and instead a column with random numbers has been used. Hence the results will be different to those reported in the paper.

### Session Info
```
R version 4.0.4 (2021-02-15)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] cowplot_1.1.1        loo_2.4.1            reshape2_1.4.4       boot_1.3-26         
 [5] posterior_1.1.0      rstan_2.21.2         StanHeaders_2.21.0-7 forcats_0.5.1       
 [9] stringr_1.4.0        purrr_0.3.4          readr_2.1.1          tidyr_1.1.4         
[13] tibble_3.1.6         ggplot2_3.3.5        tidyverse_1.3.1      lubridate_1.8.0     
[17] dplyr_1.0.7         

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.6           prettyunits_1.1.1    ps_1.6.0             assertthat_0.2.1    
 [5] utf8_1.2.2           V8_3.4.2             plyr_1.8.6           R6_2.5.1            
 [9] cellranger_1.1.0     backports_1.4.1      reprex_2.0.1         stats4_4.0.4        
[13] httr_1.4.2           pillar_1.6.4         rlang_0.4.12         curl_4.3.2          
[17] readxl_1.3.1         rstudioapi_0.13      callr_3.7.0          checkmate_2.0.0     
[21] munsell_0.5.0        tinytex_0.36         broom_0.7.11         compiler_4.0.4      
[25] modelr_0.1.8         xfun_0.29            pkgconfig_2.0.3      pkgbuild_1.2.0      
[29] tidyselect_1.1.1     tensorA_0.36.2       gridExtra_2.3        codetools_0.2-18    
[33] matrixStats_0.58.0   fansi_1.0.2          crayon_1.4.2         tzdb_0.2.0          
[37] dbplyr_2.1.1         withr_2.4.3          distributional_0.2.2 grid_4.0.4          
[41] jsonlite_1.7.3       gtable_0.3.0         lifecycle_1.0.1      DBI_1.1.2           
[45] magrittr_2.0.1       scales_1.1.1         RcppParallel_5.1.4   cli_3.1.0           
[49] stringi_1.7.6        farver_2.1.0         fs_1.5.2             xml2_1.3.3          
[53] ellipsis_0.3.2       generics_0.1.1       vctrs_0.3.8          tools_4.0.4         
[57] glue_1.6.0           hms_1.1.1            abind_1.4-5          parallel_4.0.4      
[61] processx_3.5.2       inline_0.3.18        colorspace_2.0-2     rvest_1.0.2         
[65] haven_2.4.3 ```
