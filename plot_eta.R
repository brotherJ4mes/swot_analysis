library(rasterVis)



layout(matrix(c(1,1,2,3),2,2), widths=c(.7,.4))
plot(r_eta)
par(mar=c(2,0,1,0))
xlim = c(-2,2)*100
ubar = round(global(abs(u), 'mean', na.rm=T)*100, digits=1)
vbar = round(global(abs(v), 'mean', na.rm=T)*100, digits=1)
hist(u*100, main='u (cm/s)', yaxt='n', xlab=NA, xlim=xlim, ylab=NA)
mtext(sprintf('mean: %4.1f', ubar), adj=.95, line=-2)
hist(v*100, main='v (cm/s)', yaxt='n', xlab=NA, xlim=xlim, ylab=NA)
mtext(sprintf('mean: %4.1f', vbar), adj=.95, line=-2)

vectorplot(raster::stack(c(u,v)), isField = "dXY", par.settings = viridisTheme())

