# FPWM

### Installation
```r
install.packages("devtools")
library("devtools")
devtools::install_github("https://github.com/benoukraflab/forkedTF")
```

### Quick usage
- MiniCofactor Report
```r
library("forkedTF")
miniCofactorReport( TF = "CEBPB", cell = "K562",cobinding_threshold=.06)
```
- FPWM creation and plot
```r
fpwm <- createFPWM(mainTF ="CEBPB",
                        partners = c("ATF4","ATF7","ATF3","JUND","FOS","CEBPD"),
                        cell = "K562", 
                        forkPosition = 5,
                        flipMatrix = FALSE)

plotFPWM(fpwm,pdfName="fpwm_plot.pdf")
```
- Writing FPWM
```r
write.FPWM(FPWM = fpwm, format = "transfac", fileName = "FPWM.transfact" )
write.FPWM(FPWM = fpwm, format = "FPWMtransfac", fileName = "FPWM.FPWMtransfac" )
```
