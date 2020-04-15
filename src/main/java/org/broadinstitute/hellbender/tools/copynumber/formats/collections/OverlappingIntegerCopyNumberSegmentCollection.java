package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;

import java.io.File;

public class OverlappingIntegerCopyNumberSegmentCollection extends AbstractRecordCollection<SampleLocatableMetadata, IntegerCopyNumberSegment> {
    public OverlappingIntegerCopyNumberSegmentCollection(final File inputFile) {
        super(inputFile, IntegerCopyNumberSegmentCollection.IntegerCopyNumberSegmentTableColumn.COLUMNS,
                IntegerCopyNumberSegmentCollection.INTEGER_COPY_NUMBER_SEGMENT_RECORD_DECODER,
                IntegerCopyNumberSegmentCollection.INTEGER_COPY_NUMBER_SEGMENT_RECORD_ENCODER);
    }

    @Override
    Metadata.Type getMetadataType() {
        return Metadata.Type.SAMPLE_LOCATABLE;
    }
}
