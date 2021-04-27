# some of the less common packages to install: corrplot, dendextend, tidyverse, FSA, varhandle
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
groups = hclust(dist(abs(cm))) # groups from corr_single_scalar
ordering = rev(groups$order)
cm = cm[ordering, ordering] # reorder

suppressMessages(suppressWarnings(library(dendextend)))
d = as.dendrogram(groups)
db = color_branches(d, k = 4)
dl = color_labels(db, k = 4)
# bottom left top right
postscript('clust_single_scalar.eps', horizontal = TRUE, onefile = FALSE, paper = "special", height = 15, width = 20, colormodel="rgb")
par(mar = c(0,18,0,0))
g = plot_horiz.dendrogram(rev(dl), side = TRUE, sub="", main="", axes=F) # reverse order
# plot(groups, xlab="", sub="", main="", axes=F, ylab="")
invisible(dev.off())
postscript('corr_single_scalar.eps', horizontal = FALSE, onefile = FALSE, paper = "special", height = 15, width = 15, colormodel="rgb")
par(mar = c(0,0,0,0))
corrplot(cm, type = 'upper', sig.level = 0.01, tl.cex = 0.9, tl.col = "black")
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
A = clusters["normalizedSubgraphCentrality"]
B = clusters["randomRobustnessIndex"]
C = clusters["RCB"]
D = clusters["fragmentation"]
group1 = which(clusters == A)
group2 = which(clusters == B)
group3 = which(clusters == C)
group4 = which(clusters == D)

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

g4 = ss[ss$measure %in% names(group4),]
lc = length(levels(droplevels(g4)$measure))
cat('DX', lc, dim(g4)[1], '\n')

cat(c('A', names(group1)), sep = '\n')
cat(c('B', names(group2)), sep = '\n')
cat(c('C', names(group3)), sep = '\n')
cat(c('D', names(group4)), sep = '\n')

fontsize = 15
yrange = c(0.08, 1200)
ybreaks = c(0.1, 1, 10, 100, 1000, 10000)
ylabels = c('0.1 ms', '1 ms', '10 ms', '100 ms', '1 s', '10 s')
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
ggsave(plot = p, device='eps', filename = 'poscor_g1.eps', units='cm', width = la, height = 16)
pc1 = p

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
ggsave(plot = p, filename = 'poscor_g2.eps', device='eps', unit='cm', width=lb, height=14)
pc2 = p

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
ggsave(plot = p, filename = 'poscor_g3.eps', device='eps', unit='cm', width=lc, height=14)
pc3 = p

fontsize = 15
yrange = c(0.08, 1200)
ybreaks = c(0.1, 1, 10, 100, 1000)
ylabels = c('0.1 ms', '1 ms', '10 ms', '100 ms', '1 s')
p = ggplot(g4, aes(x = measure, y = t, fill = measure)) +
    geom_boxplot(width=0.5) +
    scale_y_continuous(trans='log2',
                       limits = yrange, breaks = ybreaks,
                       labels = ylabels, name = 'Runtime') +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1),
          legend.position = "none") +
    scale_fill_discrete(name = "Measure")
ggsave(plot = p, filename = 'poscor_g4.eps', device='eps', unit='cm', width=lc, height=14)
pc4 = p

fig5a = pc1 + ggtitle("A") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 20))
fig5b = pc2 + ggtitle("B") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 20))
fig5c = pc3 + ggtitle("C") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 20))
fig5d = pc4 + ggtitle("D") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 20))
fig5 = grid.arrange(fig5a, fig5b, fig5c, fig5d, nrow = 2, ncol = 2)
ggsave(plot = fig5, filename = "fig5.eps", device = 'eps', width = 12, height = 12)

eff = c('splittingNumber',
        'RCB',
        'hubDensity',
        'connectivityRobustnessFunction',
        'electricalNodalRobustness',
        'percolatedPath')

fontsize = 26

v1 = 1
v2 = 1
v3 = 1
v4 = 1
v5 = 1
v6 = 1

genorder = c('complete_graph',
             'triangular_lattice_graph',
             'lattice.grid_2d_graph',
             'hexagonal_lattice_graph',
             'connected_watts_strogatz_graph',
             'random_regular_graph',
             'gnm_random_graph',
             'circular_ladder_graph',
             'ladder_graph',
             'barabasi_albert_graph',
             'random_powerlaw_tree',
             'star_graph',
             'barbell_graph',
             'wheel_graph',
             'path_graph')

xlab = c('CG',
         'TL',
         'SL',
         'HL',
         'WS',
         'RR',
         'ER',
         'LL',
         'CL',
         'BA',
         'PT',
         'HG',
         'BB',
         'WG',
         'PG')
for (char in eff) {
    ms = ss[ss$measure == char,]
    ms$gen = factor(ms$generator, levels = genorder)
#    print(levels(ms$gen))
    p = ggplot(ms, aes(x = gen, y = value, fill = gen)) +
        geom_boxplot(width=0.5, lwd=1.5) +
        scale_y_continuous(name = 'Reported value') +
        scale_x_discrete(name = 'Generation model', labels  = xlab) +
        theme_classic(base_size = fontsize) +
        theme(legend.position="none") +
        theme(axis.text.x = element_text(angle=90, hjust=1))
    ggsave(plot = p, filename = paste('values_', char, '_', dur, 'sec.eps', sep = ''), device = 'eps', width = 12, height = 5)
    if (char == "splittingNumber") {
        v1 = p
    } else if (char == "connectivityRobustnessFunction") {
        v2 = p
    } else if (char == "electricalNodalRobustness") {
        v3 = p
    } else if (char == "percolatedPath") {
        v4 = p
    } else if (char == "RCB") {
        v5 = p
    } else if (char == "hubDensity") {
        v6 = p
    }
}

fig6a = v1 + ggtitle("A") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 30))
fig6b = v2 + ggtitle("B") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 30))
fig6c = v3 + ggtitle("C") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 30))
fig6d = v4 + ggtitle("D") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 30))
fig6e = v5 + ggtitle("E") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 30))
fig6f = v6 + ggtitle("F") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 30))
fig6 = grid.arrange(fig6a, fig6b, fig6c, fig6d, fig6e, fig6f, nrow = 3, ncol = 2)
ggsave(plot = fig6, filename = "fig6.eps", device = 'eps', width = 20, height = 30)


palette = c("#8bbd8b","#c1cc99","#f5a65b","#5b8266","#b0daf1")
fontsize = 12
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
ggsave(plot = p, filename = paste('single_scalar_', dur, 'sec.eps', sep=""), device = 'eps', width = 7, height = 7)
fig2a = p

print('single-graph multi-valued characteristics, if any')
filename = paste('single_avg_', dur, 'sec.dat', sep="")
if (file.exists(filename)) {
    sa = read.csv(filename, sep=' ', header=FALSE)
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
    ggsave(plot = p, filename = paste('single_avg_', dur, 'sec.eps', sep=""), device = eps, width = 7, height = 7)
}

print('two-graph single-valued characteristics, if any')
filename = paste('double_scalar_', dur, 'sec.dat', sep="")
if (file.exists(filename)) {
    ds = read.csv(filename, sep=' ', header=FALSE)
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
        if (is.finite(counts[row, "perc"])) {
            p = p + annotate(geom = "label", label=sprintf("%.0f %%", counts[row, "perc"]),
                             x = row + offset, y = 1100, label.size = 1, color = "black")
        }
    }
    ggsave(plot = p, filename = paste('double_scalar_',  dur, 'sec.eps', sep=""), device = 'eps', width = 7, height = 7)
    fig2b = p
}

#library(cowplot)
#fig2a = fig2a + draw_text("A", x = 1, y = 80, hjust = 0, vjust = 1, size = 25)
                                        #fig2b = fig2b + draw_text("B", x = 1, y = 80, hjust = 0, vjust = 1, size = 25)
fig2a = fig2a + ggtitle("A") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 20))
fig2b = fig2b + ggtitle("B") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(size = 20))
fig2 = grid.arrange(fig2a, fig2b, nrow = 1)
ggsave(plot = fig2, filename = "fig2.eps", device = 'eps', width = 14, height = 7)

print('two-graph double-valued characteristics, if any')
filename = paste('double_avg_', dur, 'sec.dat', sep="")
if (file.exists(filename)) {
    da = read.csv(filename, sep=' ', header=FALSE)
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
    ggsave(plot = p, filename = paste('double_avg_', dur, 'sec.eps', sep = ''), device = 'eps', width = 7, height = 7)
}
