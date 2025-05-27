library(animint2)
library(data.table)

## Example: 4 plots, 2 selectors.
data(intreg)
signal.colors <- c(estimate="#07c107", latent="#0098ef")
breakpoint.colors <- c("1breakpoint"="#ff7d7d", "0breakpoints"='#f6f4bf')
model.linetypes <- c(margin="dotted",limit="dashed",regression="solid")
intreg$annotations$logratio <- max(intreg$sig$log)
## To get the bottom 3 plots to line up properly, we need to plot some
## geom_blanks bigger than the x range, so we calculate that here.
blank.items <- with(intreg,{
  list(segments=list(data=selection,x="min.L",y="segments"),
       error=list(data=selection,x="max.L",y="cost"),
       regression=list(data=model,x=c("min.L","max.L"),
                       y=c("min.feature","max.feature")),
       intervals=list(data=intervals,x=c("min.L","max.L"),y="feature"))
})
Lrange <- c()
for(N in names(blank.items)){
  L <- blank.items[[N]]
  Lrange <- range(c(Lrange,unlist(L$data[,L$x])),finite=TRUE)
  blank.items[[N]]$yrange <- range(unlist(L$data[,L$y]))
}
Lrange[1] <- Lrange[1]-1
Lrange[2] <- Lrange[2]+1
for(N in names(blank.items)){
  L <- blank.items[[N]]
  blank.items[[N]]$blank <- data.frame(x=Lrange, y=L$yrange)
}

## Regions with linetype indicating errors.
breaks.by.signal <- split(intreg$breaks, intreg$breaks$signal)
anns.by.signal <- split(intreg$ann, intreg$ann$signal)
error.regions.list <- list()
for(signal in names(breaks.by.signal)){
  signal.breaks <- breaks.by.signal[[signal]]
  signal.anns <- anns.by.signal[[signal]]
  signal.anns$target.breaks <-
    ifelse(signal.anns$annotation=="1breakpoint", 1, 0)
  for(model.i in 1:20){
    model.breaks <- subset(signal.breaks, segments==model.i)
    signal.anns$breaks <- NA
    for(region.i in 1:nrow(signal.anns)){
      region <- signal.anns[region.i, ]
      after.start <- region$first.base < model.breaks$base
      before.end <- model.breaks$base < region$last.base
      signal.anns$breaks[region.i] <- sum(after.start & before.end)
    }
    signal.anns$error.type <- with(signal.anns, {
      ifelse(breaks < target.breaks, "false negative",
             ifelse(target.breaks < breaks, "false positive", "correct"))
    })
    error.regions.list[[paste(model.i, signal)]] <-
      data.frame(segments=model.i, signal.anns)
  }
}
error.regions <- do.call(rbind, error.regions.list)

reg <- subset(intreg$model, line=="regression")
slope <- with(reg, (min.L-max.L)/(min.feature-max.feature))
intreg$intervals$pred.L <-
  slope * (intreg$intervals$feature - reg$min.feature) + reg$min.L

## TODO: plot reconstruction error vs model complexity! (instead of
## penalty)
library(data.table)
library(dplyr)
segs <- data.table(intreg$segments)
sigs <- data.table(intreg$signals) %>%
  mutate(base.after=base+1)
setkey(segs, signal, first.base, last.base)
setkey(sigs, signal, base, base.after)
ov <- foverlaps(sigs, segs)
model.selection <- ov %>%
  group_by(signal, segments) %>%
  summarise(error=sum((logratio-mean)^2),
            data=n()) %>%
  mutate(log.data=log(data),
         penalized.error=error + log.data * segments)
sig.labels <- model.selection %>%
  group_by(signal) %>%
  filter(segments==1)
sig.seg.names <- ov %>%
  group_by(signal, segments) %>%
  summarise(min.base=min(base),
            max.base=max(base)) %>%
  mutate(base=(min.base+max.base)/2)
sig.names <- sig.seg.names %>%
  group_by(signal, base) %>%
  summarise(logratio=max(sigs$logratio))
seg.names <- sig.seg.names %>%
  group_by(signal, segments, base) %>%
  summarise(logratio=min(sigs$logratio)-0.2)
# tallrects <- make_tallrect(model.selection, "segments")
tallrects <- geom_tallrect(aes(xmin=segments-0.5, xmax=segments+0.5
                          ),clickSelects="segments",
                      data=model.selection, fill = signal.colors[["estimate"]],
                      alpha=0.3, size=0.5)
tallrects$geom_params$colour <- signal.colors[["estimate"]]

library(reshape2)
model.tall <- melt(model.selection, measure.vars=c("error", "penalized.error"))

## Plot error AND penalized error versus number of segments.
penalized <- model.selection %>%
  group_by(signal) %>%
  filter(penalized.error < penalized.error[1]*2)
mmir.BIC <- 
  list(
    title = "Segmentation model selection with BIC",
    source = "https://github.com/suhaani-agarwal/animint/blob/master/inst/examples/intreg_BIC.R",
    segments=ggplot()+
       theme_bw()+
       ggtitle("Select profile and number of segments")+
       tallrects+
       scale_x_continuous("segments", breaks=c(1, 5, 10, 20),
                          limits=c(-2, 21))+
       xlab("")+
       facet_grid(variable ~ ., scales="free_y")+
       geom_line(aes(segments, error,
                     group=signal
                     ),clickSelects="signal",
                 data=data.frame(model.selection,
                   variable="un-penalized error"),
                 alpha=0.6, size=8)+
       geom_line(aes(segments, penalized.error,
                     group=signal
                     ),clickSelects="signal",
                 data=data.frame(penalized, variable="penalized error (BIC)"),
                 alpha=0.6, size=8),

       signal=ggplot()+
       theme_bw()+
       theme_animint(width=800)+       
       scale_x_continuous("position on chromosome (mega base pairs)",
                          breaks=c(100,200))+
       scale_fill_manual(values=breakpoint.colors,guide="none")+
       geom_blank(aes(first.base/1e6, logratio+2/8), data=intreg$ann)+
       ggtitle("Copy number profile and maximum likelihood segmentation")+
       ylab("logratio")+
       geom_point(aes(base/1e6, logratio
                      ),showSelected="signal",
                  data=intreg$sig)+
       geom_segment(aes(first.base/1e6, mean, xend=last.base/1e6, yend=mean
                        ),showSelected = c("signal", "segments"),
                    data=intreg$seg, colour=signal.colors[["estimate"]])+
       geom_segment(aes(base/1e6, min(sigs$logratio),
                        xend=base/1e6, yend=max(sigs$logratio)
                        ),showSelected = c("signal", "segments"),
                  colour=signal.colors[["estimate"]],
                  linetype="dashed",
                  data=intreg$breaks)+
       geom_text(aes(base/1e6, logratio, label=paste("signal", signal)
                     ),showSelected="signal", size = 15,
                 data=sig.names)+
       geom_text(aes(base/1e6, logratio,
                     label=sprintf("%d segment%s", segments,
                       ifelse(segments==1, "", "s"))
                     ),showSelected=c("signal", "segments"),
                 data=seg.names, color=signal.colors[["estimate"]]),
        

       first=list(signal="4.2", segments=4),
       duration=list(signal=2000, segments=2000),
       time=list(variable="segments", ms=3000),
       selector.types=list(signal="single", segments="single")
)
animint2pages(mmir.BIC, "intreg-BIC")