package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;

public final class TestingReferenceDataSource implements ReferenceDataSource{

    private final byte[] bytes;
    private final SAMSequenceDictionary dict;
    private final String contigName;

    public TestingReferenceDataSource(final String contigName, final byte[] bytes){
        Utils.nonNull(contigName, "contig name cannot be null");
        Utils.nonNull(bytes, "bytes cannot be null");
        this.bytes= bytes.clone();
        this.contigName= contigName;
        this.dict = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord(contigName, bytes.length)));
    }

    @Override
    public ReferenceSequence queryAndPrefetch(final String contig, final long start, final long stop) {
        if (contig == null || !contig.equals(contigName)){
            throw new IllegalArgumentException(String.format("contig %s is not equal to the given %s", contig, contigName));
        }
        if (start > Integer.MAX_VALUE || stop > Integer.MAX_VALUE){
            throw new IllegalArgumentException("invalid index, too high");
        }
        Utils.validIndex((int)start, bytes.length);
        Utils.validIndex((int)stop-1, bytes.length);
        if (start > stop){
            throw new IllegalArgumentException("start needs to be lower than stop");
        }

        final byte[] newbytes= Arrays.copyOfRange(bytes, (int)start, (int)stop);
        return new ReferenceSequence(contig, 0, newbytes);
    }

    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return dict;
    }

}
