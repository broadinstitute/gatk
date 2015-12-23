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
    /**
     * Write segment file based on allelic counts at SNP sites.  Converts allelic counts to target coverages,
     * which are written to a temporary file and then passed to {@link RCBSSegmenter}.
     * @param snps                  TargetCollection of allelic counts at SNP sites
     * @param sampleName            sample name
     * @param outputFile            segment file to write to and return
     */
    public static void writeSegmentFile(final TargetCollection<AllelicCount> snps,
                                        final String sampleName, final File outputFile) {
        try {
            final File targetsFromSNPCountsFile = File.createTempFile("targets-from-snps", ".tsv");

            List<TargetCoverage> targetsFromSNPCounts = snps.targets().stream()
                    .map(count -> count.toMinorAlleleFractionTargetCoverage("snp-target" + count.getContig() + ":" +
                            count.getStart() + "-" + count.getEnd())).collect(Collectors.toList());

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
