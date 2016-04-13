package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.segmenter.RCBSSegmenter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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

            final List<Target> targets = snps.targets().stream()
                    .map(ac -> new Target(name(ac), ac.getInterval())).collect(Collectors.toList());

            final RealMatrix minorAlleleFractions = new Array2DRowRealMatrix(snps.targetCount(), 1);
            minorAlleleFractions.setColumn(0, snps.targets().stream()
                    .mapToDouble(ac -> ac.estimateMinorAlleleFraction(allelicBias)).toArray());

            ReadCountCollectionUtils.write(targetsFromSNPCountsFile, new ReadCountCollection(targets, Arrays.asList(sampleName), minorAlleleFractions));

            //segment SNPs based on observed log_2 minor allele fraction (log_2 is applied in CBS.R)
            RCBSSegmenter.writeSegmentFile(sampleName, targetsFromSNPCountsFile.getAbsolutePath(),
                    outputFile.getAbsolutePath(), false);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not create temporary output file during " +
                    "SNP segmentation.", e);
        }
    }

    private static String name(final AllelicCount ac) {
        return "snp-target" + ac.getContig() + ":" + ac.getStart() + "-" + ac.getEnd();
    }
}
