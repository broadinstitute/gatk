library(HMM)
outfile = "hmm-unit-test-truth.tsv"
sequences = strsplit(readLines("hmm-unit-test-sequences.txt"), ',')
sequences.n = length(sequences)

models = list(
	initHMM(States=c("A", "B", "C"), Symbols=c("X", "Y", "Z"), startProbs=c(0.009804, 0.9804, 0.009804), transProbs=matrix(c(0.9804, 0.009804, 0.009804, 0.009804, 0.9804, 0.009804, 0.009804, 0.009804, 0.9804), nrow=3, byrow=T), emissionProbs=matrix(c(0.9804, 0.009804, 0.009804, 0.009804, 0.9804, 0.009804, 0.009804, 0.009804, 0.9804), nrow=3, byrow=T)),
	initHMM(States=c("A", "B", "C"), Symbols=c("X", "Y", "Z"), startProbs=c(0.3333, 0.3333, 0.3333), transProbs=matrix(c(0.3436, 0.3282, 0.3282, 0.3282, 0.3436, 0.3282, 0.3282, 0.3282, 0.3436), nrow=3, byrow=T), emissionProbs=matrix(c(0.3436, 0.3282, 0.3282, 0.3282, 0.3436, 0.3282, 0.3282, 0.3282, 0.3436), nrow=3, byrow=T)),
	initHMM(States=c("A", "B", "C"), Symbols=c("X", "Y", "Z"), startProbs=c(0.3333, 0.3333, 0.3333), transProbs=matrix(c(0.3333, 0.3333, 0.3333, 0.3333, 0.3333, 0.3333, 0.3333, 0.3333, 0.3333), nrow=3, byrow=T), emissionProbs=matrix(c(0.3333, 0.3333, 0.3333, 0.3333, 0.3333, 0.3333, 0.3333, 0.3333, 0.3333), nrow=3, byrow=T)),
	initHMM(States=c("A", "B", "C"), Symbols=c("X", "Y", "Z"), startProbs=c(0.001010, 0.9003, 0.09871), transProbs=matrix(c(0.9023, 0.0009894, 0.09669, 0.08895, 0.9102, 0.0008692, 0.0008727, 0.08528, 0.9138), nrow=3, byrow=T), emissionProbs=matrix(c(0.9111, 0.07935, 0.009540, 9.109e-05, 0.9109, 0.08902, 0.1821, 0.1864, 0.6315), nrow=3, byrow=T))
)
models.n = length(models)
write(x = paste('SEQUENCE', 'MODEL', 'POSITION', 'BEST_PATH', 'FW_A', 'FW_B', 'FW_C', 'BW_A', 'BW_B', 'BW_C', 'PP_A', 'PP_B', 'PP_C', sep ='\t'), file = outfile, append = F)
for (i in 1:models.n) {
    for (j in 1:sequences.n) {
        best_path = viterbi(models[[i]], sequences[[j]])
        fw = forward(models[[i]], sequences[[j]])
        bw = backward(models[[i]], sequences[[j]])
        pp = log(posterior(models[[i]], sequences[[j]]))
        for (k in 1:length(best_path)) {
            write(x = paste(j - 1, i - 1, k - 1, best_path[k], fw['A', k], fw['B', k], fw['C', k], bw['A', k], bw['B', k], bw['C', k], pp['A', k], pp['B', k], pp['C', k], sep = '\t'), file = outfile, append = T)
}}}
