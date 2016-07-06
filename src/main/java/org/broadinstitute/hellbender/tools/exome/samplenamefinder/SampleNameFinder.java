package org.broadinstitute.hellbender.tools.exome.samplenamefinder;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
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
     * @param readCountsFile a target coverage file (e.g. output of tangent normalization), where each target has a raw copy ratio estimate.  Never {@code null}
     * @return list of sample names, each which will be one of the column headers in the input file.  Never {@code null}
     */
    public static List<String> determineSampleNamesFromReadCountsFile(final File readCountsFile) {
        Utils.nonNull(readCountsFile);
        Utils.regularReadableUserFile(readCountsFile);
        final List<String> result = ReadCountCollectionUtils.retrieveSampleNamesFromReadCountsFile(readCountsFile);
        if (result.isEmpty()) {
            throw new UserException.BadInput("Input read counts file contains no sample columns");
        } else {
            return result;
        }
    }

    public static List<String> determineSampleNamesFromSegmentFile(final File segmentFile) {
        Utils.nonNull(segmentFile);
        Utils.regularReadableUserFile(segmentFile);

        return SegmentUtils.readSampleNamesFromSegmentFile(segmentFile);
    }
}
