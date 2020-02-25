## FPWM

# Installation
```r
install.packages("devtools")
library("devtools")
devtools::install_github("https://github.com/rtmag/FPWM")
```

# Quick usage
MiniCofactor Report
```r
miniCofactorReport( TF = "CEBPB", cell = "K562",cobinding_threshold=.25)
```
FPWM Creation and plot

```r
fpwm <- createFPWM(mainTF ="CEBPB",
                        partners = c("ATF4","ATF7","ATF3","JUND","FOS","CEBPD"),
                        cell = "K562", 
                        forkPosition = 5,
                        flipMatrix = FALSE)

plotFPWM(fpwm,pdfName="fpwm_plot.pdf")
```
