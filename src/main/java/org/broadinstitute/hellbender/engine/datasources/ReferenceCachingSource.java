package org.broadinstitute.hellbender.engine.datasources;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class ReferenceCachingSource implements ReferenceSource, Serializable {
    private static final long serialVersionUID = 1L;

    private ReferenceSequenceFile referenceSequenceFile;
    private Map<String, ReferenceSequence> contigs = new HashMap<>();

    public ReferenceCachingSource(String referenceURL) {
        this.referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(IOUtils.getPath(referenceURL));
    }

    @Override
    public synchronized ReferenceBases getReferenceBases(SimpleInterval interval) {
        String contig = interval.getContig();
        ReferenceSequence referenceSequence = contigs.get(contig);
        if (referenceSequence == null) {
            referenceSequence = referenceSequenceFile.getSequence(contig);
            contigs.put(contig, referenceSequence);
        }
        byte[] bytes = Arrays.copyOfRange(referenceSequence.getBases(), interval.getStart(), interval.getEnd() + 1);
        return new ReferenceBases(bytes, interval);
    }

    @Override
    public SAMSequenceDictionary getReferenceSequenceDictionary(SAMSequenceDictionary optReadSequenceDictionaryToMatch) {
        return referenceSequenceFile.getSequenceDictionary();
    }

    @Override
    public boolean isCompatibleWithSparkBroadcast() {
        return true;
    }
}
