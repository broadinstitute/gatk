package org.broadinstitute.hellbender.tools.exome.samplenamefinder;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.ReadCountCollectionUtils;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.Reader;
import java.util.List;

/**
 * Static class for finding the sample name from the various tsv file formats.
 */
public class SampleNameFinder {
    private SampleNameFinder() {}

    /**
     * Useful for determining a sample name to give the segmenter.
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

    /**
     * Useful for determining a sample name to give the segmenter.
     *
     * @param readCountsReader a target coverage reader (e.g. output of tangent normalization), where each target has a raw copy ratio estimate.  Never {@code null}
     * @param readCountsSourceName a string identifier for {@param readCountsReader} to use in error messages.  Never {@code null}
     * @return list of sample names, each which will be one of the column headers in the input file.  Never {@code null}
     */
    public static List<String> determineSampleNamesFromReadCountsReader(final Reader readCountsReader, final String readCountsSourceName) {
        Utils.nonNull(readCountsReader);
        Utils.nonNull(readCountsSourceName);
        final List<String> result = ReadCountCollectionUtils.retrieveSampleNamesFromReadCountsReader(readCountsReader, readCountsSourceName);
        if (result.isEmpty()) {
            throw new UserException.BadInput("Input read counts file contains no sample columns");
        } else {
            return result;
        }
    }

    /**
     * Read only the sample names from a segment file.
     *
     * @param segmentFile Never {@code null}
     * @return a list of the unique sample names in the seg file.
     */    public static List<String> determineSampleNamesFromSegmentFile(final File segmentFile) {
        return SegmentUtils.readSampleNamesFromSegmentFile(segmentFile);
    }

    /**
     * Read only the sample names from a segment reader.
     *
     * @param segmentReader Never {@code null}; an instance of {@link Reader} for input segments
     * @param segmentSourceName Never {@code null}; a string identifier for the input segments reader (used in error messages)
     * @return a list of the unique sample names in the seg file.
     */
    public static List<String> determineSampleNamesFromSegmentReader(final Reader segmentReader, final String segmentSourceName) {
        return SegmentUtils.readSampleNamesFromSegmentReader(segmentReader, segmentSourceName);
    }

}
