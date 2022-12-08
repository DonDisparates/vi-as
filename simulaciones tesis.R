########### Simulaciones en la tesis #############

library(copula)
library(gridExtra)
library(lattice)
library(scatterplot3d)
library(ggplot2)
library(grid)

####### Copulas Fundamentales ##############################

u <- seq(0,1, length.out = 40)
u12 <- expand.grid("u[1]"=u, "u[2]"=u)


W <- pmax(u12[,1]+u12[,2]-1,0)
M <- pmin(u12[,1],u12[,2])

val.W <- cbind(u12, "W(u[1],u[2])"=W)
val.M <- cbind(u12, "M(u[1],u[2])"=M)

##Cópulas de comonotonicidas y contracomonotnicidad

w1 <- wireframe2(val.W, main = "Contracomonotonicidad", zlab=NULL, shade =T)
w2 <- wireframe2(val.M, main = "Comonotonicidad", zlab=NULL)
c1 <- contourplot2(val.W,xlim= 0:1,ylim =0:1)
c2 <- contourplot2(val.M,xlim= 0:1,ylim =0:1)


#Cópula de independencia

d <- 2
ic <- indepCopula(dim = d)
ui <- runif(d) # puntos aleatorios 
(Pi <- pCopula(ui, copula = ic)) # simulamos la cópula de independencia
w3 <-wireframe2  (ic, FUN = pCopula, main = "Independencia")
c3 <-contourplot2(ic, FUN = pCopula) 

#Gráficamos
grid.arrange(w1,c1,ncol=1)
grid.arrange(w2,c2,ncol=1)
grid.arrange(w3,c3,ncol=1)


###Copula Gaussiana #############

normal  <- normalCopula ( param  =  0.7 , dim  =  2 )
gaus <- rCopula(2000, normal)
pdf_n <- dCopula(gaus, normal)

par(mfrow = c(1, 3))
scatterplot3d(gaus[,1], gaus[,2], pdf_n, color="red",xlab ="a)", ylab="", zlab="", pch=".")
persp(normal, dCopula, col = "red", xlab="", ylab="",zlab="")
contour(normal, dCopula, xlim = c(0, 1), ylim=c(0, 1), xlab="c)", ylab="")
mtext(expression(paste("Cópula Gaussiana con ",rho,"=0.7",sep="")),line = -3,side = 3,outer = TRUE, cex =2)
mtext("b)",line = -25,side = 3,outer = TRUE, cex =1)


##### Copula t #####

t1  <- tCopula ( param  =  0.8 , dim  =  2 , df  =  2 )
t <- rCopula (2000, t1)
pdf_t <- dCopula(t, t1)
par(mfrow = c(1, 3))
scatterplot3d(t[,1], t[,2], pdf_t, color="red", xlab ="a)", ylab="", zlab="", pch=".")
persp(t1, dCopula, col = "red", xlab="", ylab="",zlab="")
contour(t1, dCopula, xlim = c(0, 1), ylim=c(0, 1),xlab="c)", ylab="")
mtext(expression(paste("Cópula t-Student con ",rho,"=0.8 y ",italic(v)," =2",sep="")),line = -3,side = 3,outer = TRUE, cex =2)
mtext("b)",line = -25,side = 3,outer = TRUE, cex =1)


##### Copulas elipticas comparacion ####

##Primero ejecutar codigo 1 y 2

p0 <-  qplot(gaus[,1], gaus[,2],col = gaus[,1] ,main="", xlab = "a)", ylab = "")
p4 <- qplot(t[,1], t[,2], col = t[,1],main="", xlab = "b)", ylab = "")
grid.arrange(p0,p4 ,ncol=2, top = textGrob("Comparación de cópulas elípticas",gp=gpar(fontsize=20,font=3)))



### Copula Clayton ####

clayton  <- claytonCopula ( dim  =  2 , param  =  19 )
cl <- rCopula(2000, clayton)
pdf_cl <- dCopula(cl, clayton)

par(mfrow = c(1, 3))
scatterplot3d(cl[,1], cl[,2], pdf_cl, color="blue",xlab="a)", ylab="",zlab="", pch=".")
persp(clayton, dCopula, col = "blue",xlab="", ylab="",zlab="")
contour(clayton, dCopula, xlim = c(0, 1), ylim=c(0, 1),xlab="c)", ylab="")
mtext(expression(paste("Cópula Clayton con ",theta,"=19",sep="")),line = -3,side = 3,outer = TRUE, cex =2)
mtext("b)",line = -25,side = 3,outer = TRUE, cex =1)

### Copula Gumbel ####


gumbel  <- gumbelCopula ( dim  =  2 , param  =  5.6 )
gu <- rCopula(2000, gumbel)
pdf_gu <- dCopula(gu,gumbel)

par(mfrow = c(1, 3))
scatterplot3d(gu[,1], gu[,2], pdf_gu, color="blue",xlab ="a)", ylab="", zlab="", pch=".")
persp(gumbel, dCopula, col = "blue",xlab="", ylab="",zlab="")
contour(gumbel, dCopula, xlim = c(0, 1), ylim=c(0, 1),xlab="c)", ylab="")
mtext(expression(paste("Cópula Gumbel con ",theta,"= 5.6",sep="")),line = -3,side = 3,outer = TRUE, cex =2)
mtext("b)",line = -25,side = 3,outer = TRUE, cex =1)

###copula frank ###

frank  <- frankCopula ( dim  =  2 , param  =  8 )
fr <- rCopula(2000, frank)
pdf_fr <- dCopula(fr,frank)


par(mfrow = c(1, 3))
scatterplot3d(fr[,1], fr[,2], pdf_fr, color="blue",xlab ="a)", ylab="", zlab="", pch=".")
persp(frank, dCopula, col = "blue",xlab ="", ylab="", zlab="")
contour(frank, dCopula, xlim = c(0, 1), ylim=c(0, 1),xlab ="c)", ylab="")
mtext(expression(paste("Cópula Frank con ",theta,"= 8",sep="")),line = -3,side = 3,outer = TRUE, cex =2)
mtext("b)",line = -25,side = 3,outer = TRUE, cex =1)

##comparacion de copulas aqrquimedianas###

#primero ejecutar los codigos 5,6,7

p1 <- qplot(cl[,1], cl[,2], colour = cl[,1], main="", xlab = "a)", ylab = "")
p2 <- qplot(gu[,1], gu[,2], colour = gu[,1], main="", xlab = "b)", ylab = "") 
p3 <- qplot(fr[,1], fr[,2], colour = fr[,1], main="", xlab = "c)", ylab = "")

grid.arrange(p1,p2,p3 ,ncol=3, top = textGrob("Comparación de Cópulas Arquimedianas",gp=gpar(fontsize=20,font=3)))




################ Cópulas rotadas ##############



############Clayton Rotada #############
rc0  <- claytonCopula ( dim  =  2 , param  =  2)
rc_0 <- rCopula(2000, rc0)

#0 grados
pr00 <-qplot(rc_0[,1], rc_0[,2], main="0 grados", xlab = "", ylab = "")
pr0<-contourplot2(rc0,FUN=dCopula, main="0 grados", xlab ="", ylab="")

#90 grados

rc90 <-rotCopula(claytonCopula(dim = 2,param = 2),flip = c(TRUE,FALSE))
rc_90 <- rCopula(2000,rc90)
pr1<- qplot(rc_90[,1], rc_90[,2], main="90 grados", xlab ="", ylab="")
pr2 <-contourplot2(rc90,FUN=dCopula,main="90 grados",xlab ="", ylab="")

#180 grados

rc180 <-rotCopula(claytonCopula(dim = 2,param = 2),flip = c(TRUE,TRUE))
rc_180 <- rCopula(2000,rc180)
pr3 <-qplot(rc_180[,1], rc_180[,2], main="180 grados", xlab ="", ylab="")
pr4 <-contourplot2(rc180,FUN=dCopula,main="180 grados",xlab ="", ylab="")

#270 grados

rc270 <-rotCopula(claytonCopula(dim = 2,param = 2),flip = c(FALSE,TRUE))
rc_270 <- rCopula(2000,rc270)
pr5 <-qplot(rc_270[,1], rc_270[,2], main="270 grados", xlab ="", ylab="")
pr6 <- contourplot2(rc270,FUN=dCopula,main="270 grados",xlab ="", ylab="")

grid.arrange(pr00, pr1,pr3,pr5,pr0, pr2,pr4, pr6, nrow =2, ncol=4)


############ Gumbel Rotada #############
rg0  <- gumbelCopula ( dim  =  2 , param  =  6)
rg_0 <- rCopula(2000, rc0)

#0 grados
pr001 <-qplot(rg_0[,1], rg_0[,2], main="0 grados", xlab = "", ylab = "")
pr01<-contourplot2(rg0,FUN=dCopula, main="0 grados", xlab ="", ylab="")

#90 grados

rg90 <-rotCopula(gumbelCopula(dim = 2,param = 6),flip = c(TRUE,FALSE))
rg_90 <- rCopula(2000,rg90)
pr11<- qplot(rc_90[,1], rc_90[,2], main="90 grados", xlab ="", ylab="")
pr21 <-contourplot2(rg90,FUN=dCopula,main="90 grados",xlab ="", ylab="")

#180 grados

rg180 <-rotCopula(gumbelCopula(dim = 2,param = 6),flip = c(TRUE,TRUE))
rg_180 <- rCopula(2000,rg180)
pr31 <-qplot(rg_180[,1], rg_180[,2], main="180 grados", xlab ="", ylab="")
pr41 <-contourplot2(rg180,FUN=dCopula,main="180 grados",xlab ="", ylab="")

#270 grados

rg270 <-rotCopula(gumbelCopula(dim = 2,param = 6),flip = c(FALSE,TRUE))
rg_270 <- rCopula(2000,rg270)
pr51 <-qplot(rg_270[,1], rg_270[,2], main="270 grados", xlab ="", ylab="")
pr61 <- contourplot2(rg270,FUN=dCopula,main="270 grados",xlab ="", ylab="")

grid.arrange(pr001, pr11,pr31,pr51,pr01, pr21,pr41, pr61, nrow =2, ncol=4)




## correlació de rango####



#### tau ##

## Cuatro copulas con marginales N(0,1), y tau =0.7
tau <- 0.7
th.n <- iTau(normalCopula(),  tau = tau)
th.t <- iTau(tCopula(df = 3), tau = tau)
th.c <- iTau(claytonCopula(), tau = tau)
th.g <- iTau(gumbelCopula(),  tau = tau)
set.seed(271)
n <- 10000
N01m <- list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)) # margins
X.n <- rMvdc(n, mvdc = mvdc(normalCopula(th.n),    c("norm", "norm"), N01m))
X.t <- rMvdc(n, mvdc = mvdc(tCopula(th.t, df = 3), c("norm", "norm"), N01m))
X.c <- rMvdc(n, mvdc = mvdc(claytonCopula(th.c),   c("norm", "norm"), N01m))
X.g <- rMvdc(n, mvdc = mvdc(gumbelCopula(th.g),    c("norm", "norm"), N01m))

plotCorners <- function(X, qu, lim, smooth = FALSE, main, ...)
{
  plot(X, xlim = lim, ylim = lim, xlab = main, ylab = "",
       col = adjustcolor("black", 0.5), ...) # or pch = 16
  abline(h = qu, v = qu, lty = 2, col = adjustcolor("black", 0.6))
  ll <- sum(apply(X <= qu[1], 1, all)) * 100 / n
  ur <- sum(apply(X >= qu[2], 1, all)) * 100 / n
  invisible()
}
## Plots
a. <- 0.005
q <- qnorm(c(a., 1 - a.)) # a- and (1-a)-quantiles of N(0,1)
lim <- range(q, X.n, X.t, X.c, X.g)
lim <- c(floor(lim[1]), ceiling(lim[2]))

par(mfrow= c(2,2))
plotCorners(X.n, qu = q, lim = lim, cex = 0.4, main="a)")
plotCorners(X.t, qu = q, lim = lim, cex = 0.4, main="b)")
plotCorners(X.c, qu = q, lim = lim, cex = 0.4, main="c)")
plotCorners(X.g, qu = q, lim = lim, cex = 0.4, main="d)")
mtext(expression(paste("Copulas con marginales N(0,1) y ",tau,"=0.7",sep="")),line = -3,side = 3,outer = TRUE, cex =2)







