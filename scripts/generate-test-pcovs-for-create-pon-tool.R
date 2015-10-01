# script to generate test combined pcovs for CreatePanelOfNormals tool.
#
# arguments:

output.file = "all.pcovs"
sample.n = 100
target.n = 10000

# Some parameters to control the distribution of sample, target and noise factor in the output counts.

sample.scale = 150
sample.dispersion = .1

target.scale = 50
target.dispersion = .1

noise.scale = 10
noise.dispersion = 2.

# Then you can create the PoN like this:
#
# java -Xmx12g -Djava.library.path=/whereever-hdf5-libs-are-installed -jar hellbender-all-lib.jar CreatePanelOfNormals -I all.pcovs -o pon.hd5

sample.names = paste("sample_", 1:sample.n, sep="")
target.names = paste("target_", 1:target.n, sep="")

sample.factor = rnorm(sample.n)
target.factor = rnorm(target.n)

noise.factor = rnorm(sample.n * target.n)

counts = noise.factor * noise.dispersion * noise.scale + noise.scale +
      rep.int(sample.factor * sample.scale * sample.dispersion + sample.scale, target.n) +
  unlist(lapply(1:target.n, function(i) rep.int(target.factor[i] * target.scale * target.dispersion + target.scale , sample.n)))

# Normalize by sample totals
counts.mat = matrix(counts, nrow=target.n, ncol=sample.n, byrow=T)
col.sums = colSums(counts.mat)
col.sums.mat = matrix(rep.int(col.sums, nrow(counts.mat)), byrow=T, nrow=nrow(counts.mat))
counts.mat = counts.mat / col.sums.mat
counts.df = as.data.frame(counts.mat)
names(counts.df) = sample.names
row.names(counts.df) = target.names

# the target coord and name appear in the last columns.
counts.df$contig = rep("1", nrow(counts.df))
counts.df$start = as.integer(1:target.n * 1000)
counts.df$end = as.integer(1:target.n * 1000 + 200)
counts.df$name = target.names

write.table(counts.df, sep="\t", file=output.file, quote=F, row.names=F)
