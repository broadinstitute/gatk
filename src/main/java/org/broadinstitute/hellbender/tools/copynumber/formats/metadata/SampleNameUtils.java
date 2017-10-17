package org.broadinstitute.hellbender.tools.copynumber.formats.metadata;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SampleNameUtils {
    public static final String SAMPLE_NAME_COMMENT_TAG = "SAMPLE_NAME=";
    public static final String SAMPLE_NAME_COMMENT_PREFIX = TableUtils.COMMENT_PREFIX + SAMPLE_NAME_COMMENT_TAG;

    private SampleNameUtils() {}

    public static String readSampleName(final File file) {
        IOUtils.canReadFile(file);
        final List<String> sampleNameCommentLines = new ArrayList<>();
        try (final XReadLines reader = new XReadLines(file)) {
            for (final String line : reader) {
                if (!line.startsWith(SAMPLE_NAME_COMMENT_PREFIX)) {
                    break;
                }
                sampleNameCommentLines.add(line);
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(file);
        }
        if (sampleNameCommentLines.size() != 1) {
            throw new UserException.BadInput(String.format("File does not contain exactly one sample name specified by %s.",
                    SAMPLE_NAME_COMMENT_PREFIX));
        }
        return sampleNameCommentLines.get(0).replace(SAMPLE_NAME_COMMENT_PREFIX, "");
    }

    public static String readSampleName(final SAMFileHeader header) {
        if (header == null || header.getReadGroups() == null) {
            throw new UserException.BadInput("The input BAM has no header or no read groups. Cannot determine a sample name.");
        }
        final List<String> sampleNames = header.getReadGroups().stream().map(SAMReadGroupRecord::getSample).distinct().collect(Collectors.toList());
        if (sampleNames.size() > 1) {
            throw new UserException.BadInput(String.format("The input BAM has more than one unique sample name: %s", StringUtils.join(sampleNames, ", ")));
        }
        if (sampleNames.isEmpty()) {
            throw new UserException.BadInput("The input BAM has no sample names.");
        }
        return sampleNames.get(0);
    }
}
