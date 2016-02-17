package org.broadinstitute.hellbender.tools.exome.samplenamefinder;

import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.exome.TargetCoverageUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.List;

/**
 * Static class for finding the sample name from the various tsv file formats.
 */
public class SampleNameFinder {
    private SampleNameFinder() {}

    /**
     *  Useful for determining a sample name to give the segmenter.
     *
     * @param targetCoverageFile a target coverage file (e.g. output of tangent normalization), where each target has a raw copy ratio estimate.  Never {@code null}
     * @return list of sample names, each which will be one of the column headers in the input file.  Never {@code null}
     */
    public static List<String> determineSampleNamesFromTargetCoverageFile(final File targetCoverageFile) {
        Utils.nonNull(targetCoverageFile);
        Utils.regularReadableUserFile(targetCoverageFile);

        return TargetCoverageUtils.retrieveSampleNamesFromTargetCoverageFile(targetCoverageFile);
    }

    public static List<String> determineSampleNamesFromSegmentFile(final File segmentFile) {
        Utils.nonNull(segmentFile);
        Utils.regularReadableUserFile(segmentFile);

        return SegmentUtils.readSampleNamesFromSegmentFile(segmentFile);
    }
}
