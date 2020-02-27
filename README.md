# FPWM

### Installation
```r
install.packages("devtools")
library("devtools")
devtools::install_github("https://github.com/rtmag/FPWM")
```

### Quick usage
- MiniCofactor Report
```r
library("FPWM")
miniCofactorReport( TF = "CEBPB", cell = "K562",cobinding_threshold=.06)
```
<div align="center">
<a name="miniCofactorReport"/>
<img src="./inst/MM1_HSA_K562_CEBPB_cofactor_minireport.png" alt="miniCofactorReport" width="490" height="630" ></img>
</a>
</div>

- FPWM creation and plot
```r
fpwm <- createFPWM(mainTF ="CEBPB",
                        partners = c("ATF4","ATF7","ATF3","JUND","FOS","CEBPD"),
                        cell = "K562", 
                        forkPosition = 5,
                        flipMatrix = FALSE)

plotFPWM(fpwm,pdfName="fpwm_plot.pdf")
```
<div align="center">
<a name="fpwm_plot"/>
<img src="./inst/fpwm_plot.png" alt="fpwm_plot" width="400" height="630" ></img>
</a>
</div>

- Writing FPWM
```r
write.FPWM(FPWM = fpwm, format = "transfac", fileName = "FPWM.transfact" )
write.FPWM(FPWM = fpwm, format = "FPWMtransfac", fileName = "FPWM.FPWMtransfac" )
```
