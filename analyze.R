args = commandArgs(trailingOnly=TRUE)
dur = as.numeric(args[1])
cat(sprintf("Analyzing experiments with %d second time limit\n", dur))

require(ggplot2)
require(scales)
require(varhandle)

palette = c("#8bbd8b","#c1cc99","#f5a65b","#5b8266","#b0daf1")
fontsize = 17
offset = 0.3
yrange = c(0.01, 1200)
ybreaks = c(0.01, 0.1, 1, 10, 100, 1000)

ss = read.csv(paste('single_scalar_', dur, 'sec.dat', sep=""), sep=' ', header=FALSE)
names(ss) = c('generator', 'order', 'measure', 'value', 'runtime')
ss$generator = as.factor(ss$generator)
ss$measure = as.factor(ss$measure)
ss$n = as.factor(ss$order)
ss$t = ss$runtime * 1000

p = ggplot(ss, aes(x = n, y = t, fill = n)) +
    geom_violin(trim = TRUE) + geom_boxplot(width=0.1) + scale_fill_manual(values = palette) +
    scale_y_continuous(trans='log2', labels = scientific, limits = yrange, breaks = ybreaks, name = 'Runtime (ms)') +
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
ggsave(paste('single_scalar_runtime_', dur, 'sec.png', sep=""))

sa = read.csv(paste('single_avg_', dur, 'sec.dat', sep=""), sep=' ', header=FALSE)
names(sa) = c('generator', 'order', 'measure', 'value', 'runtime')
sa$generator = as.factor(sa$generator)
sa$measure = as.factor(sa$measure)
sa$n = as.factor(sa$order)
sa$t = sa$runtime * 1000
p = ggplot(sa, aes(x = n, y = t, fill = n)) +
    geom_violin(trim = TRUE) + geom_boxplot(width=0.1) + scale_fill_manual(values = palette) +
    scale_y_continuous(trans='log2', labels = scientific, limits = yrange, breaks = ybreaks, name = 'Runtime (ms)') +
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
ggsave(paste('single_avg_runtime_', dur, 'sec.png', sep=""))

ds = read.csv(paste('double_scalar_', dur, 'sec.dat', sep=""), sep=' ', header=FALSE)
names(ds) = c('first', 'second', 'order', 'measure', 'value', 'runtime')
ds$measure = as.factor(ds$measure)
ds$n = as.factor(ds$order)
ds$t = ds$runtime * 1000
p = ggplot(ds, aes(x = n, y = t, fill = n)) +
    geom_violin(trim = TRUE) + geom_boxplot(width=0.1) + scale_fill_manual(values = palette) +
    scale_y_continuous(trans='log2', labels = scientific, limits = yrange, breaks = ybreaks, name = 'Runtime (ms)') +
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
ggsave(paste('double_scalar_runtime_',  dur, '.png', sep=""))

da = read.csv(paste('double_avg_', dur, 'sec.dat', sep=""), sep=' ', header=FALSE)
names(da) = c('first', 'second', 'order', 'measure', 'value', 'runtime')
da$measure = as.factor(da$measure)
da$n = as.factor(da$order)
da$t = da$runtime * 1000
p = ggplot(da, aes(x = n, y = t, fill = n)) +
    geom_violin(trim = TRUE) + geom_boxplot(width=0.1) + scale_fill_manual(values = palette) +
    scale_y_continuous(trans='log2', labels = scientific, limits = yrange, breaks = ybreaks, name = 'Runtime (ms)') +
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
ggsave(paste('double_avg_runtime_', dur, 'sec.png', sep = ''))

