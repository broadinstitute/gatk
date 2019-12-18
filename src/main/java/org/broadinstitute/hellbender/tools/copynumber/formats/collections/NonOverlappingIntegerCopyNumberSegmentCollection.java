package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Represents a collection of {@link IntegerCopyNumberSegment} for a sample, guaranteed to be non-overlapping.
 */
public final class NonOverlappingIntegerCopyNumberSegmentCollection extends IntegerCopyNumberSegmentCollection {
    public NonOverlappingIntegerCopyNumberSegmentCollection(final File inputFile) {
        super(inputFile);
        CopyNumberArgumentValidationUtils.validateNonOverlappingIntervals(getRecords(), getMetadata().getSequenceDictionary());
    }

    NonOverlappingIntegerCopyNumberSegmentCollection(final SampleLocatableMetadata metadata,
                                                            final List<IntegerCopyNumberSegment> integerCopyNumberSegmentList) {
        super(metadata, integerCopyNumberSegmentList);
        CopyNumberArgumentValidationUtils.validateNonOverlappingIntervals(getRecords(), getMetadata().getSequenceDictionary());
    }
}
