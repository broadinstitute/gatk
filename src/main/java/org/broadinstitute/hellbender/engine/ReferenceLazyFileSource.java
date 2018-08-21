package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.Serializable;
import java.util.Iterator;

public class ReferenceLazyFileSource implements ReferenceDataSource, Serializable {
    private static final long serialVersionUID = 1L;

    private String path;
    private transient ReferenceFileSource referenceFileSource;

    public ReferenceLazyFileSource(String path) {
        this.path = path;
        this.referenceFileSource = new ReferenceFileSource(IOUtils.getPath(path));
    }

    public ReferenceFileSource getReferenceFileSource() {
        if (referenceFileSource == null) {
            referenceFileSource = new ReferenceFileSource(IOUtils.getPath(path));
        }
        return referenceFileSource;
    }

    @Override
    public ReferenceSequence queryAndPrefetch(String contig, long start, long stop) {
        return getReferenceFileSource().queryAndPrefetch(contig, start, stop);
    }

    @Override
    public SAMSequenceDictionary getSequenceDictionary() {
        return getReferenceFileSource().getSequenceDictionary();
    }

    @Override
    public Iterator<Byte> iterator() {
        return getReferenceFileSource().iterator();
    }

    @Override
    public void close() {
        getReferenceFileSource().close();
    }

    public ReferenceSequenceFile getReferenceSequenceFile() {
        return getReferenceFileSource().getReferenceSequenceFile();
    }
}
