DEBUG = FALSE
if (!DEBUG) {
    options(warn=-1)
}
suppressMessages(suppressWarnings(library(corrplot)))
d = read.csv('cor.dat', sep=' ', header=TRUE, row.names = 1)
if (DEBUG) {
    sapply(d, typeof)
}
omit = c()
for (i in names(d)) {
    columnData = as.numeric(unlist(d[i])) # yeah, R is "fun"
    if (sd(columnData, na.rm = TRUE) == 0) {
        cat(i, ' ', summary(columnData), '\n')
        omit = c(omit, i)
    }
}
cat('omitted', length(omit), '\n')
variable = d[, -which(names(d) %in% omit)]
cat('kept', dim(variable)[1], '\n')
cm = cor(variable, method = 'pearson', use = 'pairwise.complete.obs')
groups = hclust(dist(abs(cm))) # groups from corr_single_scalar.png

suppressMessages(suppressWarnings(library(dendextend)))
d = as.dendrogram(groups)
db = color_branches(d, k = 4)
dl = color_labels(db, k = 4)
# bottom left top right
png('clust_single_scalar.png', width = 800, height = 1000)
par(mar = c(0,18,0,0))
plot_horiz.dendrogram(dl, side = TRUE, sub="", main="", axes=F)
# plot(groups, xlab="", sub="", main="", axes=F, ylab="")
invisible(dev.off())
png('corr_single_scalar.png', width = 1000, height = 1000)
par(mar = c(0,0,0,0))
co = rev(groups$order)
ocm = cm[co, co]
corrplot(ocm, type = 'upper', sig.level = 0.01, tl.cex = 0.9, tl.col = "black")
invisible(dev.off())

suppressMessages(suppressWarnings(require(ggplot2)))
suppressMessages(suppressWarnings(require(scales)))
suppressMessages(suppressWarnings(require(varhandle)))
suppressMessages(suppressWarnings(library(tidyverse))) # incl. dplyr
suppressMessages(suppressWarnings(require(ggrepel)))
suppressMessages(suppressWarnings(require(gridExtra)))
suppressMessages(suppressWarnings(require(cowplot)))

args = commandArgs(trailingOnly=TRUE)
dur = as.numeric(args[1])
cat(sprintf("Analyzing experiments with %d second time limit\n", dur))
ss = read.csv(paste('single_scalar_', dur, 'sec.dat', sep=""), sep=' ', header=FALSE)
names(ss) = c('generator', 'order', 'measure', 'value', 'runtime', 'replica')
ss$generator = as.factor(ss$generator)
ss$measure = as.factor(ss$measure)
ss$n = as.factor(ss$order)
ss$t = ss$runtime * 1000
ss$value = as.numeric(ss$value)

clusters = cutree(groups, k = 4)
A = clusters["vertexLoad"]
B = clusters["robustnessIndex"]
C = clusters["fragility"]
group1 = which(clusters == A)
group2 = which(clusters == B)
group3 = which(clusters == C)

cat('Median runtime', median(ss$t), '\n')
rf = ss[ss$measure == 'resilienceFactor',]
mrg = median(rf$t)
summary(rf$t)
suppressMessages(suppressWarnings(require(FSA)))
cat(perc(ss$t, mrg, 'lt'), 'are less\n')

g1 = ss[ss$measure %in% names(group1),]
la = length(levels(droplevels(g1)$measure))
cat('A',  la, dim(g1)[1], '\n')

g2 = ss[ss$measure %in% names(group2),]
lb = length(levels(droplevels(g2)$measure))
cat('B', lb, dim(g2)[1], '\n')


g3 = ss[ss$measure %in% names(group3),]
lc = length(levels(droplevels(g3)$measure))
cat('C', lc, dim(g3)[1], '\n')


cat(c('A', names(group1)), sep = '\n')
cat(c('B', names(group2)), sep = '\n')
cat(c('C', names(group3)), sep = '\n')


fontsize = 15
yrange = c(0.08, 1200)
ybreaks = c(0.1, 1, 10, 100, 1000)
ylabels = c('0.1 ms', '1 ms', '10 ms', '100 ms', '1 s')
p = ggplot(g1, aes(x = measure, y = t, fill = measure)) +
    geom_boxplot(width=0.5) +
    scale_y_continuous(trans='log2',
                       limits = yrange, breaks = ybreaks,
                       labels = ylabels, name = 'Runtime') +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
          legend.position = "none") +
    scale_fill_discrete(name = "Measure")
                                        #    guides(fill = guide_legend(ncol = 2))
                                        #           axis.ticks.x=element_blank(),
                                        #           axis.text.x=element_blank(),
ggsave('poscor_g1.png', unit='cm', width=la, height=16)

fontsize = 15
yrange = c(0.08, 1200)
ybreaks = c(0.1, 1, 10, 100, 1000)
ylabels = c('0.1 ms', '1 ms', '10 ms', '100 ms', '1 s')
p = ggplot(g2, aes(x = measure, y = t, fill = measure)) +
    geom_boxplot(width=0.5) +
    scale_y_continuous(trans='log2',
                       limits = yrange, breaks = ybreaks,
                       labels = ylabels, name = 'Runtime') +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
          legend.position = "none") +
    scale_fill_discrete(name = "Measure")
ggsave('poscor_g2.png', unit='cm', width=lb, height=14)


fontsize = 15
yrange = c(0.08, 1200)
ybreaks = c(0.1, 1, 10, 100, 1000)
ylabels = c('0.1 ms', '1 ms', '10 ms', '100 ms', '1 s')
p = ggplot(g3, aes(x = measure, y = t, fill = measure)) +
    geom_boxplot(width=0.5) +
    scale_y_continuous(trans='log2',
                       limits = yrange, breaks = ybreaks,
                       labels = ylabels, name = 'Runtime') +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
          legend.position = "none") +
    scale_fill_discrete(name = "Measure")
ggsave('poscor_g3.png', unit='cm', width=lc, height=14)


palette = c("#8bbd8b","#c1cc99","#f5a65b","#5b8266","#b0daf1")
fontsize = 17
offset = 0.3
yrange = c(0.01, 1200)
ybreaks = c(0.1, 1, 10, 100, 1000)
ylabels = c('0.1 ms', '1 ms', '10 ms', '100 ms', '1 s')
p = ggplot(ss, aes(x = n, y = t, fill = n)) +
    geom_violin(trim = TRUE) + geom_boxplot(width=0.1) + scale_fill_manual(values = palette) +
    scale_y_continuous(trans='log2',
                       limits = yrange, breaks = ybreaks,
                       labels = ylabels, name = 'Runtime') +
    scale_x_discrete(name = 'Graph order') +
    theme_classic(base_size = fontsize) + theme(legend.position="none")
counts = as.data.frame(table(ss$n))
names(counts) = c('n', 'actual')
k = length(levels(ss$generator))
h = length(levels(ss$measure))
counts$x = unfactor(counts$n)
counts$replicas = 10 - ceiling(log(counts$x, 2) / 2) + 1
counts$expected = counts$replicas * k * h
counts$perc = 100 * counts$actual / counts$expected
for (row in 1:nrow(counts)) {
    p = p + annotate(geom = "label", label=sprintf("%.0f %%", counts[row, "perc"]),
                     x = row + offset, y = 1100, label.size = 1, color = "black")
}
ggsave(paste('single_scalar_', dur, 'sec.png', sep=""), width = 7, height = 7)

sa = read.csv(paste('single_avg_', dur, 'sec.dat', sep=""), sep=' ', header=FALSE)
names(sa) = c('generator', 'order', 'measure', 'value', 'runtime')
sa$generator = as.factor(sa$generator)
sa$measure = as.factor(sa$measure)
sa$n = as.factor(sa$order)
sa$t = sa$runtime * 1000
p = ggplot(sa, aes(x = n, y = t, fill = n)) +
    geom_violin(trim = TRUE) + geom_boxplot(width=0.1) + scale_fill_manual(values = palette) +
    scale_y_continuous(trans='log2', labels = ylabels, limits = yrange, breaks = ybreaks, name = 'Runtime') +
    scale_x_discrete(name = 'Graph order') +
    theme_classic(base_size = fontsize) + theme(legend.position="none")
counts = as.data.frame(table(sa$n))
names(counts) = c('n', 'actual')
k = length(levels(sa$generator))
h = length(levels(sa$measure))
counts$x = unfactor(counts$n)
counts$replicas = 10 - ceiling(log(counts$x, 2) / 2) + 1
counts$expected = counts$replicas * k * h
counts$perc = 100 * counts$actual / counts$expected
for (row in 1:nrow(counts)) {
    p = p + annotate(geom = "label", label=sprintf("%.0f %%", counts[row, "perc"]),
                     x = row + offset, y = 1100, label.size = 1, color = "black")
}
ggsave(paste('single_avg_', dur, 'sec.png', sep=""), width = 7, height = 7)

ds = read.csv(paste('double_scalar_', dur, 'sec.dat', sep=""), sep=' ', header=FALSE)
names(ds) = c('first', 'second', 'order', 'measure', 'value', 'runtime')
ds$measure = as.factor(ds$measure)
ds$n = as.factor(ds$order)
ds$t = ds$runtime * 1000
p = ggplot(ds, aes(x = n, y = t, fill = n)) +
    geom_violin(trim = TRUE) + geom_boxplot(width=0.1) + scale_fill_manual(values = palette) +
    scale_y_continuous(trans='log2', labels = ylabels, limits = yrange, breaks = ybreaks, name = 'Runtime') +
    scale_x_discrete(name = 'Graph order') +
     theme_classic(base_size = fontsize) + theme(legend.position="none")
counts = as.data.frame(table(ds$n))
names(counts) = c('n', 'actual')
k1 = length(levels(ds$first))
k2 = length(levels(ds$second))
h = length(levels(ds$measure))
counts$x = unfactor(counts$n)
counts$replicas = 10 - ceiling(log(counts$x, 2) / 2) + 1
counts$expected = counts$replicas * k1 * k2 * h
counts$perc = 100 * counts$actual / counts$expected
for (row in 1:nrow(counts)) {
    p = p + annotate(geom = "label", label=sprintf("%.0f %%", counts[row, "perc"]),
                     x = row + offset, y = 1100, label.size = 1, color = "black")
}
ggsave(paste('double_scalar_',  dur, 'sec.png', sep=""), width = 7, height = 7)

da = read.csv(paste('double_avg_', dur, 'sec.dat', sep=""), sep=' ', header=FALSE)
names(da) = c('first', 'second', 'order', 'measure', 'value', 'runtime')
da$measure = as.factor(da$measure)
da$n = as.factor(da$order)
da$t = da$runtime * 1000
p = ggplot(da, aes(x = n, y = t, fill = n)) +
    geom_violin(trim = TRUE) + geom_boxplot(width=0.1) + scale_fill_manual(values = palette) +
    scale_y_continuous(trans='log2', labels = ylabels, limits = yrange, breaks = ybreaks, name = 'Runtime') +
    scale_x_discrete(name = 'Graph order') +
    theme_classic(base_size = fontsize) + theme(legend.position="none")
counts = as.data.frame(table(da$n))
names(counts) = c('n', 'actual')
k1 = length(levels(da$first))
k1 = length(levels(da$second))
h = length(levels(da$measure))
counts$x = unfactor(counts$n)
counts$replicas = 10 - ceiling(log(counts$x, 2) / 2) + 1
counts$expected = counts$replicas * k1 * k2 * h
counts$perc = 100 * counts$actual / counts$expected
for (row in 1:nrow(counts)) {
    p = p + annotate(geom = "label", label=sprintf("%.0f %%", counts[row, "perc"]),
                     x = row + offset, y = 1100, label.size = 1, color = "black")
}
ggsave(paste('double_avg_', dur, 'sec.png', sep = ''), width = 7, height = 7)

eff = c('closedWalkNumber',
        'normalizedSubgraphCentrality',
        'redundancyOfAlternativePaths',
        'fragmentation',
        'relativeEntropy',
        'pathDiversity',
        'percolatedPath',
        'perturbationScore')

fontsize = 20
xlab = c('BA', 'BB', 'CL', 'CG', 'ER', 'SL', 'HL', 'LG', 'PG', 'RG', 'HG', 'TL', 'WG')
for (char in eff) {
    ms = ss[ss$measure == char,]
#    print(levels(ms$generator))
    p = ggplot(ms, aes(x = generator, y = value, fill = generator)) +
        geom_boxplot(width=0.5, lwd=2) +
        scale_y_continuous(name = 'Value') +
        scale_x_discrete(name = 'Generation model', labels  = xlab) +
        theme_classic(base_size = fontsize) +
        theme(legend.position="none") +
        theme(axis.text.x = element_text(angle=90, hjust=1))
    ggsave(paste('values_', char, '_', dur, 'sec.png', sep = ''), width = 12, height = 5)
}
