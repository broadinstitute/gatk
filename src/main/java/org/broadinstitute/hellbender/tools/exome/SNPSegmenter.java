package org.broadinstitute.hellbender.tools.exome;

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
     * @param snpCounts             list of allelic counts at SNP sites
     * @param sampleName            sample name
     * @param outputFile            segment file to write to and return
     * @throws IOException          if temporary target file cannot be created or written to
     */
    public static void writeSegmentFile(final List<AllelicCount> snpCounts,
                                        final String sampleName, final File outputFile)
            throws IOException {
        final File targetsFromSNPCountsFile = File.createTempFile("targets-from-snps", ".tsv");

        List<TargetCoverage> targetsFromSNPCounts = snpCounts.stream()
                .map(count -> count.toMinorAlleleFractionTargetCoverage("snp-target" + count.getContig() + ":" +
                        count.getStart() + "_" + count.getEnd())).collect(Collectors.toList());

        TargetCoverageUtils.writeTargetsWithCoverage(targetsFromSNPCountsFile, sampleName, targetsFromSNPCounts);

        RCBSSegmenter.writeSegmentFile(sampleName, targetsFromSNPCountsFile.getAbsolutePath(),
                outputFile.getAbsolutePath(), false);
    }
}
