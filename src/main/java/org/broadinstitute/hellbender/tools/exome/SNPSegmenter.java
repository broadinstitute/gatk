package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Finds segments based on allelic counts at SNP sites.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SNPSegmenter {
    private SNPSegmenter() {
    }

    /**
     * Write segment file based on maximum-likelihood estimates of the minor allele fraction at SNP sites,
     * assuming no allelic bias.  These estimates are converted to target coverages,
     * which are written to a temporary file and then passed to {@link RCBSSegmenter}.
     * @param snps                  TargetCollection of allelic counts at SNP sites
     * @param sampleName            sample name
     * @param outputFile            segment file to write to and return
     */
    public static void writeSegmentFile(final TargetCollection<AllelicCount> snps,
                                        final String sampleName, final File outputFile) {
        writeSegmentFile(snps, sampleName, outputFile, 1.);
    }

    /**
     * Write segment file based on maximum-likelihood estimates of the minor allele fraction at SNP sites,
     * assuming the specified allelic bias.  These estimates are converted to target coverages,
     * which are written to a temporary file and then passed to {@link RCBSSegmenter}.
     * @param snps                  TargetCollection of allelic counts at SNP sites
     * @param sampleName            sample name
     * @param outputFile            segment file to write to and return
     * @param allelicBias           allelic bias to use in estimate of minor allele fraction
     */
    public static void writeSegmentFile(final TargetCollection<AllelicCount> snps, final String sampleName,
                                        final File outputFile, final double allelicBias) {
        try {
            final File targetsFromSNPCountsFile = File.createTempFile("targets-from-snps", ".tsv");

            List<TargetCoverage> targetsFromSNPCounts = snps.targets().stream()
                    .map(count -> count.toMinorAlleleFractionTargetCoverage(
                            "snp-target" + count.getContig() + ":" + count.getStart() + "-" + count.getEnd(),
                            allelicBias))
                    .collect(Collectors.toList());

            TargetCoverageUtils.writeTargetsWithCoverage(targetsFromSNPCountsFile, sampleName, targetsFromSNPCounts);

            //segment SNPs based on observed log_2 minor allele fraction (log_2 is applied in CBS.R)
            RCBSSegmenter.writeSegmentFile(sampleName, targetsFromSNPCountsFile.getAbsolutePath(),
                    outputFile.getAbsolutePath(), false);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not create temporary output file during " +
                    "SNP segmentation.", e);
        }
    }
}
