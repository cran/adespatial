## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.align = "center")

## -----------------------------------------------------------------------------
library(ade4)
library(adespatial)
library(adegraphics)
library(spdep)
library(maptools)

## -----------------------------------------------------------------------------
data("mafragh")
class(mafragh)
names(mafragh)
dim(mafragh$flo)

## -----------------------------------------------------------------------------
str(mafragh$env)

## ---- fig.height = 4, fig.width = 4-------------------------------------------
mxy <- as.matrix(mafragh$xy)
rownames(mxy) <- NULL
s.label(mxy, ppoint.pch = 15, ppoint.col = "darkseagreen4", Sp = mafragh$Spatial.contour)

## -----------------------------------------------------------------------------
mafragh$spenames[c(1, 11), ]

## ---- fig.height=3, fig.width=6-----------------------------------------------
fpalette <- colorRampPalette(c("white", "darkseagreen2", "darkseagreen3", "palegreen4"))
sp.flo <- SpatialPolygonsDataFrame(Sr = mafragh$Spatial, data = mafragh$flo, match.ID = FALSE)
s.Spatial(sp.flo[,c(1, 11)], col = fpalette(3), nclass = 3)

## -----------------------------------------------------------------------------
data(mafragh)
class(mafragh$Spatial)
nb.maf <- poly2nb(mafragh$Spatial)
s.Spatial(mafragh$Spatial, nb = nb.maf, plabel.cex = 0, pnb.edge.col = 'red')

## -----------------------------------------------------------------------------
xygrid <- expand.grid(x = 1:10, y = 1:8)
s.label(xygrid, plabel.cex = 0)

## -----------------------------------------------------------------------------
nb2.q <- cell2nb(8, 10, type = "queen")
nb2.r <- cell2nb(8, 10, type = "rook")
s.label(xygrid, nb = nb2.q, plabel.cex = 0, main = "Queen neighborhood")
s.label(xygrid, nb = nb2.r, plabel.cex = 0, main = "Rook neighborhood")

## -----------------------------------------------------------------------------
xytransect <- expand.grid(1:20, 1)
nb3 <- cell2nb(20, 1)

summary(nb3)

## -----------------------------------------------------------------------------
set.seed(3)
xyir <- mxy[sample(1:nrow(mafragh$xy), 20),]
s.label(xyir, main = "Irregular sampling with 20 sites")

## ---- fig.width = 5-----------------------------------------------------------
nbnear1 <- dnearneigh(xyir, 0, 50)
nbnear2 <- dnearneigh(xyir, 0, 305)

g1 <- s.label(xyir, nb = nbnear1, pnb.edge.col = "red", main = "neighbors if 0<d<50", plot = FALSE)
g2 <- s.label(xyir, nb = nbnear2, pnb.edge.col = "red", main = "neighbors if 0<d<305", plot = FALSE)
cbindADEg(g1, g2, plot = TRUE)

## -----------------------------------------------------------------------------
nbnear1

## -----------------------------------------------------------------------------
nbnear2

## ---- fig.width = 5-----------------------------------------------------------
knn1 <- knearneigh(xyir, k = 1)
nbknn1 <- knn2nb(knn1, sym = TRUE)
knn2 <- knearneigh(xyir, k = 2)
nbknn2 <- knn2nb(knn2, sym = TRUE)

g1 <- s.label(xyir, nb = nbknn1, pnb.edge.col = "red", main = "Nearest neighbors (k=1)", plot = FALSE)
g2 <- s.label(xyir, nb = nbknn2, pnb.edge.col = "red", main = "Nearest neighbors (k=2)", plot = FALSE)
cbindADEg(g1, g2, plot = TRUE)

## -----------------------------------------------------------------------------
n.comp.nb(nbknn1)

## -----------------------------------------------------------------------------
nbtri <- tri2nb(xyir)
nbgab <- graph2nb(gabrielneigh(xyir), sym = TRUE)
nbrel <- graph2nb(relativeneigh(xyir), sym = TRUE)
nbsoi <- graph2nb(soi.graph(nbtri, xyir), sym = TRUE)

g1 <- s.label(xyir, nb = nbtri, pnb.edge.col = "red", main = "Delaunay", plot = FALSE)
g2 <- s.label(xyir, nb = nbgab, pnb.edge.col = "red", main = "Gabriel", plot = FALSE)
g3 <- s.label(xyir, nb = nbrel, pnb.edge.col = "red", main = "Relative", plot = FALSE)
g4 <- s.label(xyir, nb = nbsoi, pnb.edge.col = "red", main = "Sphere of influence", plot = FALSE)

ADEgS(list(g1, g2, g3, g4))

## -----------------------------------------------------------------------------
nbgab[[1]]

## -----------------------------------------------------------------------------
diffnb(nbsoi, nbrel)

## -----------------------------------------------------------------------------
nbgab <- graph2nb(gabrielneigh(mxy), sym = TRUE)

## -----------------------------------------------------------------------------
nb2listw(nbgab)

## -----------------------------------------------------------------------------
distgab <- nbdists(nbgab, mxy)
nbgab[[1]]
distgab[[1]]

## -----------------------------------------------------------------------------
fdist <- lapply(distgab, function(x) 1 - x/max(dist(mxy)))

## -----------------------------------------------------------------------------
listwgab <- nb2listw(nbgab, glist = fdist)
listwgab
names(listwgab)
listwgab$neighbours[[1]]
listwgab$weights[[1]]

## -----------------------------------------------------------------------------
print(listw2mat(listwgab)[1:10, 1:10], digits = 3)

## -----------------------------------------------------------------------------
mem.gab <- mem(listwgab)
mem.gab

## -----------------------------------------------------------------------------
class(mem.gab)
names(attributes(mem.gab))

## ---- echo = -1---------------------------------------------------------------
par(mar = c(0, 2, 3, 0))
    barplot(attr(mem.gab, "values"), 
        main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7)

## ---- eval = FALSE------------------------------------------------------------
#  plot(mem.gab[,c(1, 5, 10, 20, 30, 40, 50, 60, 70)], SpORcoords = mxy)

## ---- fig.width = 5, fig.height = 5-------------------------------------------
s.value(mxy, mem.gab[,c(1, 5, 10, 20, 30, 40, 50, 60, 70)], symbol = "circle", ppoint.cex = 0.6)

## -----------------------------------------------------------------------------
moranI <- moran.randtest(mem.gab, listwgab, 99)

## -----------------------------------------------------------------------------
head(attr(mem.gab, "values") / moranI$obs)

## ---- fig.height=6, fig.width=8, out.width="80%"------------------------------
sp.env <- SpatialPolygonsDataFrame(Sr = mafragh$Spatial, data = mafragh$env, match.ID = FALSE)
maps.env <- s.Spatial(sp.env, col = fpalette(6), nclass = 6)
MC.env <- moran.randtest(mafragh$env, listwgab, nrepet = 999)
MC.env

## -----------------------------------------------------------------------------
mc.bounds <- moran.bounds(listwgab)
mc.bounds

## ---- fig = TRUE--------------------------------------------------------------
env.maps <- s1d.barchart(MC.env$obs, labels = MC.env$names, plot = FALSE, xlim = 1.1 * mc.bounds, paxes.draw = TRUE, pgrid.draw = FALSE)
addline(env.maps, v = mc.bounds, plot = TRUE, pline.col = 'red', pline.lty = 3)

## -----------------------------------------------------------------------------
NP.Mg <- moranNP.randtest(mafragh$env[,5], listwgab, nrepet = 999, alter = "two-sided") 
NP.Mg
plot(NP.Mg)
sum(NP.Mg$obs)
MC.env$obs[5]

## -----------------------------------------------------------------------------
pca.hell <- dudi.pca(mafragh$flo, scale = FALSE, scannf = FALSE, nf = 2)

## ---- fig.height=3, fig.width=6-----------------------------------------------
moran.randtest(pca.hell$li, listw = listwgab)
s.value(mxy, pca.hell$li, Sp = mafragh$Spatial.contour, symbol = "circle", col = c("white", "palegreen4"), ppoint.cex = 0.6)

## -----------------------------------------------------------------------------
ms.hell <- multispati(pca.hell, listw = listwgab, scannf = F)

## -----------------------------------------------------------------------------
summary(ms.hell)

## ---- fig.height=3, fig.width= 6----------------------------------------------
g.ms.maps <- s.value(mafragh$xy, ms.hell$li, Sp = mafragh$Spatial.contour, symbol = "circle", col = c("white", "palegreen4"), ppoint.cex = 0.6)

## ---- fig.width = 5, fig.height = 5-------------------------------------------
g.ms.spe <- s.arrow(ms.hell$c1, plot = FALSE)
g.abund <- s.value(mxy, mafragh$flo[, c(12,11,31,16)],
    Sp = mafragh$Spatial.contour, symbol = "circle", col = c("black", "palegreen4"), plegend.drawKey = FALSE, ppoint.cex = 0.4, plot = FALSE)
p1 <- list(c(0.05, 0.65), c(0.01, 0.25), c(0.74, 0.58), c(0.55, 0.05))
for (i in 1:4)
g.ms.spe <- insert(g.abund[[i]], g.ms.spe, posi = p1[[i]], ratio = 0.25, plot = FALSE)
g.ms.spe

## -----------------------------------------------------------------------------
scalo <- scalogram(mafragh$flo[,11], mem.gab)
plot(scalo)

## -----------------------------------------------------------------------------
sum(scalo$obs)

## -----------------------------------------------------------------------------
plot(scalogram(mafragh$flo[,11], mem(listwgab), nblocks = 20))

## ---- fig.height=5, fig.width=5-----------------------------------------------
mspa.hell <- mspa(pca.hell, listwgab, scannf = FALSE, nf = 2)

g.mspa <- scatter(mspa.hell, posieig = "topright", plot = FALSE)
g.mem <- s.value(mafragh$xy, mem.gab[, c(1, 2, 6, 3)], Sp = mafragh$Spatial.contour, ppoints.cex = 0.4, plegend.drawKey = FALSE, plot = FALSE)
g.abund <- s.value(mafragh$xy, mafragh$flo[, c(31,54,25)], Sp = mafragh$Spatial.contour, symbol = "circle", col = c("black", "palegreen4"), plegend.drawKey = FALSE, ppoint.cex = 0.4, plot = FALSE)

p1 <- list(c(0.01, 0.44), c(0.64, 0.15), c(0.35, 0.01), c(0.15, 0.78))
for (i in 1:4)
g.mspa <- insert(g.mem[[i]], g.mspa, posi = p1[[i]], plot = FALSE)

p2 <- list(c(0.27, 0.54), c(0.35, 0.35), c(0.75, 0.31))
for (i in 1:3)
g.mspa <- insert(g.abund[[i]], g.mspa, posi = p2[[i]], plot = FALSE)

g.mspa

## -----------------------------------------------------------------------------
mem.gab.sel <- mem.select(pca.hell$tab, listw = listwgab)
mem.gab.sel$global.test
mem.gab.sel$summary

## -----------------------------------------------------------------------------
class(mem.gab.sel$MEM.select)
dim(mem.gab.sel$MEM.select)

## -----------------------------------------------------------------------------
cand.lw <- listw.candidates(mxy, nb = c("gab", "rel"), weights = c("bin", "flin"))

## -----------------------------------------------------------------------------
sel.lw <- listw.select(pca.hell$tab, candidates = cand.lw, nperm = 99)

## -----------------------------------------------------------------------------
sel.lw$candidates
sel.lw$best.id

## -----------------------------------------------------------------------------
lw.best <- cand.lw[[sel.lw$best.id]]

## -----------------------------------------------------------------------------
sel.lw$best

## -----------------------------------------------------------------------------
rda.hell <- pcaiv(pca.hell, sel.lw$best$MEM.select, scannf = FALSE)

## -----------------------------------------------------------------------------
test.rda <- randtest(rda.hell)
test.rda
plot(test.rda)

## ---- fig.height=3, fig.width= 6----------------------------------------------
s.value(mxy, rda.hell$li, Sp = mafragh$Spatial.contour, symbol = "circle", col = c("white", "palegreen4"), ppoint.cex = 0.6)

## ---- fig.height = 5, fig.width = 5-------------------------------------------
library(vegan)
vp1 <- varpart(pca.hell$tab, mafragh$env, sel.lw$best$MEM.select)
vp1
plot(vp1, bg = c(3, 5), Xnames = c("environment", "spatial"))

## -----------------------------------------------------------------------------
vp2 <- varipart(pca.hell$tab, mafragh$env, sel.lw$best$MEM.select)
vp2

## -----------------------------------------------------------------------------
cor(mafragh$env[,10], mafragh$env[,11])

## -----------------------------------------------------------------------------
cor.test(mafragh$env[,10], mafragh$env[,11])

## -----------------------------------------------------------------------------
moran.randtest(mafragh$env[,10], lw.best)

## -----------------------------------------------------------------------------
msr1 <- msr(mafragh$env[,10], lw.best)
summary(moran.randtest(msr1, lw.best, nrepet = 2)$obs)

msr2 <- msr(mafragh$env[,11], lw.best)

## -----------------------------------------------------------------------------
obs <- cor(mafragh$env[,10], mafragh$env[,11])
sim <- sapply(1:ncol(msr1), function(i) cor(msr1[,i], msr2[,i]))
testmsr <- as.randtest(obs = obs, sim = sim, alter = "two-sided")
testmsr

## -----------------------------------------------------------------------------
msr(vp2, listwORorthobasis = lw.best)

