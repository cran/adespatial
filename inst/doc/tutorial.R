## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(fig.align = "center")

## ------------------------------------------------------------------------
library(adespatial)
library(ade4)
library(adegraphics)
library(spdep)
library(maptools)

## ------------------------------------------------------------------------
data(mafragh)
class(mafragh$Spatial)
nb.maf <- poly2nb(mafragh$Spatial)
s.Spatial(mafragh$Spatial, nb = nb.maf, plabel.cex = 0, pnb.edge.col = 'red')

## ------------------------------------------------------------------------
xygrid <- expand.grid(x = 1:10, y = 1:8)
s.label(xygrid, plabel.cex = 0)

## ------------------------------------------------------------------------
nb2.q <- cell2nb(8, 10, type = "queen")
nb2.r <- cell2nb(8, 10, type = "rook")
s.label(xygrid, nb = nb2.q, plabel.cex = 0, main = "Queen neighborhood")
s.label(xygrid, nb = nb2.r, plabel.cex = 0, main = "Rook neighborhood")

## ------------------------------------------------------------------------
xytransect <- expand.grid(1:20, 1)
nb3 <- cell2nb(20, 1)

summary(nb3)

## ------------------------------------------------------------------------
set.seed(3)
xyir <- matrix(runif(40), 20, 2)
s.label(xyir, main = "Irregular sampling with 10 sites")

## ---- fig.width = 5------------------------------------------------------
nbnear1 <- dnearneigh(xyir, 0, 0.2)
nbnear2 <- dnearneigh(xyir, 0, 1.5)

g1 <- s.label(xyir, nb = nbnear1, pnb.edge.col = "red", main = "neighbors if 0<d<0.2", plot = FALSE)
g2 <- s.label(xyir, nb = nbnear2, pnb.edge.col = "red", main = "neighbors if 0<d<1.5", plot = FALSE)
cbindADEg(g1, g2, plot = TRUE)

## ------------------------------------------------------------------------
nbnear1

## ------------------------------------------------------------------------
nbnear2

## ---- fig.width = 5------------------------------------------------------
knn1 <- knearneigh(xyir, k = 1)
nbknn1 <- knn2nb(knn1, sym = TRUE)
knn2 <- knearneigh(xyir, k = 2)
nbknn2 <- knn2nb(knn2, sym = TRUE)

g1 <- s.label(xyir, nb = nbknn1, pnb.edge.col = "red", main = "Nearest neighbors (k=1)", plot = FALSE)
g2 <- s.label(xyir, nb = nbknn2, pnb.edge.col = "red", main = "Nearest neighbors (k=2)", plot = FALSE)
cbindADEg(g1, g2, plot = TRUE)

## ------------------------------------------------------------------------
n.comp.nb(nbknn1)

## ------------------------------------------------------------------------
nbtri <- tri2nb(xyir)
nbgab <- graph2nb(gabrielneigh(xyir), sym = TRUE)
nbrel <- graph2nb(relativeneigh(xyir), sym = TRUE)
nbsoi <- graph2nb(soi.graph(nbtri, xyir), sym = TRUE)

g1 <- s.label(xyir, nb = nbtri, pnb.edge.col = "red", main = "Delaunay", plot = FALSE)
g2 <- s.label(xyir, nb = nbgab, pnb.edge.col = "red", main = "Gabriel", plot = FALSE)
g3 <- s.label(xyir, nb = nbrel, pnb.edge.col = "red", main = "Relative", plot = FALSE)
g4 <- s.label(xyir, nb = nbsoi, pnb.edge.col = "red", main = "Sphere of influence", plot = FALSE)

ADEgS(list(g1, g2, g3, g4))

## ------------------------------------------------------------------------
nbgab[[1]]

## ------------------------------------------------------------------------
diffnb(nbsoi, nbrel)

## ------------------------------------------------------------------------
nb2listw(nbgab)

## ------------------------------------------------------------------------
distgab <- nbdists(nbgab, xyir)
nbgab[[1]]
distgab[[1]]

## ------------------------------------------------------------------------
fdist <- lapply(distgab, function(x) 1 - x/max(dist(xyir)))

## ------------------------------------------------------------------------
listwgab <- nb2listw(nbgab, glist = fdist, style = "B")
listwgab
names(listwgab)
listwgab$neighbours[[1]]
listwgab$weights[[1]]

## ------------------------------------------------------------------------
print(listw2mat(listwgab)[1:10, 1:10], digits = 3)

## ------------------------------------------------------------------------
mem.gab <- mem(listwgab)
mem.gab

## ------------------------------------------------------------------------
str(mem.gab)

## ---- echo = -1----------------------------------------------------------
par(mar = c(0, 2, 3, 0))
    barplot(attr(mem.gab, "values"), 
        main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7)

## ------------------------------------------------------------------------
plot(mem.gab)

## ---- fig.width = 5, fig.height = 5--------------------------------------
plot(mem.gab, SpORcoords = xyir, nb = nbgab)

## ------------------------------------------------------------------------
moranI <- moran.randtest(mem.gab, listwgab, 99)
moranI

## ------------------------------------------------------------------------
attr(mem.gab, "values") / moranI$obs

## ---- fig.width = 6, fig.height = 4--------------------------------------
signi <- which(moranI$pvalue < 0.05)
signi
plot(mem.gab[,signi], SpORcoords = xyir, nb = nbgab)

## ------------------------------------------------------------------------
data("mafragh")
class(mafragh)
names(mafragh)
dim(mafragh$flo)

## ------------------------------------------------------------------------
str(mafragh$env)

## ---- fig.height = 4, fig.width = 4--------------------------------------
s.label(mafragh$xy, plabel.cex = 0, ppoint.pch = 15, ppoint.col = "darkseagreen4", Sp = mafragh$Spatial)

## ------------------------------------------------------------------------
mafragh$spenames[c(1, 11), ]

## ---- fig.height=3, fig.width=6------------------------------------------
fpalette <- colorRampPalette(c("white", "darkseagreen2", "darkseagreen3", "palegreen4"))
sp.flo <- SpatialPolygonsDataFrame(Sr = mafragh$Spatial, data = mafragh$flo, match.ID = FALSE)
s.Spatial(sp.flo[,c(1, 11)], col = fpalette(3), nclass = 3)

## ------------------------------------------------------------------------
lw <- nb2listw(mafragh$nb)

## ---- fig.height=6, fig.width=8, out.width="80%"-------------------------
sp.env <- SpatialPolygonsDataFrame(Sr = mafragh$Spatial, data = mafragh$env, match.ID = FALSE)
maps.env <- s.Spatial(sp.env, col = fpalette(4), nclass = 4)
MC.env <- moran.randtest(mafragh$env, lw, nrepet = 999)
MC.env

## ------------------------------------------------------------------------
mc.bounds <- moran.bounds(lw)
mc.bounds

## ---- fig = TRUE---------------------------------------------------------
env.maps <- s1d.barchart(MC.env$obs, labels = MC.env$names, plot = FALSE, ylim = 1.1 * mc.bounds, p1d.hori = FALSE)
addline(env.maps, h = mc.bounds, plot = TRUE, pline.col = 'red', pline.lty = 3)

