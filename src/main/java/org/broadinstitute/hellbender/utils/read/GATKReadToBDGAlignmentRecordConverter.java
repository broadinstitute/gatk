package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.bdgenomics.formats.avro.AlignmentRecord;
import org.bdgenomics.adam.converters.SAMRecordConverter;
import org.bdgenomics.adam.models.SequenceDictionary;
import org.bdgenomics.adam.models.RecordGroupDictionary;

/**
 * Converts a GATKRead to a BDG AlignmentRecord
 */
public class GATKReadToBDGAlignmentRecordConverter {
    private static final SAMRecordConverter converter = new SAMRecordConverter();

    private SAMFileHeader header;
    private SequenceDictionary dict;
    private RecordGroupDictionary readGroups;

    public GATKReadToBDGAlignmentRecordConverter(SAMFileHeader header) {
        this.header = header;
        this.dict = SequenceDictionary.fromSAMSequenceDictionary(header.getSequenceDictionary());
        this.readGroups = RecordGroupDictionary.fromSAMHeader(header);
    }

    public static AlignmentRecord convert( final GATKRead gatkRead, final SAMFileHeader header ) {
        SequenceDictionary dict = SequenceDictionary.fromSAMSequenceDictionary(header.getSequenceDictionary());
        RecordGroupDictionary readGroups = RecordGroupDictionary.fromSAMHeader(header);
        return GATKReadToBDGAlignmentRecordConverter.convert(gatkRead, header, dict, readGroups);
    }

    public static AlignmentRecord convert(
            final GATKRead gatkRead, final SAMFileHeader header, final SequenceDictionary dict, final RecordGroupDictionary readGroups ) {
        return converter.convert(gatkRead.convertToSAMRecord(header));
    }

    public static AlignmentRecord convert(
            final SAMRecord sam, final SequenceDictionary dict, final RecordGroupDictionary readGroups ) {
        return converter.convert(sam);
    }
}
