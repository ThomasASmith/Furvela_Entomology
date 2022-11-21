# Effect of growing (or shrinking) population on parous rate

P_0 =0.6 # true survival per cycle
R = 2


P_0_1 = (-(R + P_0) + sqrt((R + P_0)^2 - 4*R * P_0))/(-2*R)

steps =1000
Rvector = rep(seq(1:steps), each=steps)/(steps/4)
Pvector = rep(seq(1:steps), times=steps)/steps

parous_rate = function(R,P_0){
  value = ifelse(R<1,
                 (-(R + P_0) - sqrt((R + P_0)^2 - 4*R * P_0))/(-2*R),
  (-(R + P_0) + sqrt((R + P_0)^2 - 4*R * P_0))/(-2*R))
  value =   (-(R + P_0) + sqrt((R + P_0)^2 - 4*R * P_0))/(-2*R)
  value= P_0/R
  return(value)}

bias = parous_rate(R=Rvector,P_0=Pvector) - Pvector

toPlot = data.frame(Rvector,Pvector,bias)

#
library(ggplot2)
library(grid)
library(scales)
# define colourBlind safe colors
# slightly more blue pink/purple #C879C8 to replace pink#CC79A7
mp_recomCol = c("#D55E00", "#009E73", "#0072A7","#C879C8")
fmt_dcimals <- function(decimals=0){function(x) format(x,nsmall = decimals,scientific = FALSE)}
toPlot$bias[toPlot$Pvector > toPlot$Rvector] <- NA
toPlot$bias <- cut(toPlot$bias,breaks = c(-Inf,-0.75,-0.45,-0.15,0.15,0.45,Inf))
ggplot(data=toPlot,aes(x=Rvector, y=Pvector, fill=bias)) + theme_bw() +
  theme(text = element_text(size = 14)) +
  geom_tile() +
  scale_fill_manual(values = c("black","darkblue", "blue",
             "green","yellow","red","lightgrey"),na.translate = FALSE) +
            labs(x = "R", y=expression("P"[0])) 
        geom_polygon(data=triangle,mapping=aes(x=x,y=y, fill='lightgrey'), show.legend = FALSE)
