package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;
import java.nio.file.Path;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.IOException;
import java.util.List;
import java.util.Map;

/*
A class for logging likelihood matrices, in a more consise way - we just show
for read the best haplotype and the differences in the scores from the second best and the reference haplotype
 */

public class ConciseAlleleLikelihoodWriter extends AlleleLikelihoodWriter{
    public ConciseAlleleLikelihoodWriter(final Path _outputPath) {
        super(_outputPath, null);
    }
    /**
     * Write read x haplotype likelihood matrix as a matrix of
     * @param likelihoods - matrix to add
     */
    public void writeAlleleLikelihoodsAsMatrix(final AlleleLikelihoods<GATKRead, Haplotype> likelihoods,
                                               Map<String, String> haplotypeToNameMap,
                                               boolean writeHeader, int readCount){
        final List<String> samples = likelihoods.samples();
        final List<Haplotype> haplotypes = likelihoods.alleles();
        try {
            for (int s = 0 ; s < samples.size(); s++) {
                if (writeHeader){
                    output.write("Read\t");
                    output.write("Best_hap\t");
                    output.write("Best_score\t");
                    output.write("Diff_from_second\t");
                    output.write("Diff_from_ref\n");
                }

                List<GATKRead> reads = likelihoods.sampleEvidence(s);
                for (int read = 0; read < likelihoods.sampleMatrix(s).evidenceCount(); read++) {
                    output.write(String.format("%s\t", reads.get(read).getName()));
                    String bestHap = "";
                    double bestScore = Double.NEGATIVE_INFINITY;
                    double secondBestScore = Double.NEGATIVE_INFINITY;
                    double refScore = Double.NEGATIVE_INFINITY;
                    for (int allele = 0; allele < likelihoods.sampleMatrix(s).numberOfAlleles(); allele++) {
                        double newScore = likelihoods.sampleMatrix(s).get(allele, read);
                        if (newScore > bestScore){
                            secondBestScore = bestScore;
                            bestScore = newScore;
                            bestHap = haplotypeToNameMap.get(haplotypes.get(allele).toString());
                            if (haplotypes.get(allele).isReference()){
                                refScore = newScore;
                            }
                        } else if (newScore > secondBestScore){
                            secondBestScore = newScore;
                        }
                    }
                    output.write(bestHap + "\t");
                    output.write(String.format("%.03f", bestScore) + "\t");
                    output.write(String.format("%.03f", bestScore-secondBestScore) + "\t");
                    output.write(String.format("%.03f", bestScore-refScore) + "\n");
                }
            }
            output.flush();
        } catch (IOException err) {
            throw new RuntimeException("Unable to write matrix to file");
        }
    }
}
