package org.broadinstitute.hellbender.tools.spark.linkedreads;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerShort;
import org.broadinstitute.hellbender.tools.spark.sv.SVKmerizer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

public class UpstreamKmerReadFilter extends ReadFilter {
    private static final long serialVersionUID = 1L;

    private final ReferenceMultiSource reference;
    public Set<String> kmers = new HashSet<>();



    public UpstreamKmerReadFilter(final Set<String> kmers, final ReferenceMultiSource referenceMultiSource) {
        this.kmers = kmers;
        this.reference = referenceMultiSource;
    }

    @Override
    public boolean test(final GATKRead r) {
        final SAMSequenceDictionary referenceSequenceDictionary = reference.getReferenceSequenceDictionary(null);
        final SAMSequenceRecord sequence = referenceSequenceDictionary.getSequence(r.getContig());
        if (sequence == null) {
            return false;
        }

        if (r.isReverseStrand() && (r.getUnclippedEnd() < 1 || r.getUnclippedEnd() > sequence.getSequenceLength() - 7)) {
            return false;
        }

        if (!r.isReverseStrand() && (r.getUnclippedStart() < 8 || r.getUnclippedStart() > sequence.getSequenceLength() - 7)) {
            return false;
        }

        final SimpleInterval motifInterval = new SimpleInterval(r.getContig(), r.isReverseStrand() ? r.getUnclippedEnd() + 1 : r.getUnclippedStart() - 7,
                r.isReverseStrand() ? r.getUnclippedEnd() + 7 : r.getUnclippedStart() - 1 );
        if (motifInterval.getStart() < 1) return false;

        if (motifInterval.getEnd() >= sequence.getSequenceLength()) {
            return false;
        }

        final ReferenceBases referenceBases;
        try {
            referenceBases = reference.getReferenceBases(null, motifInterval);
        } catch (IOException e) {
            throw new GATKException("Could not get reference bases for " + r);
        }
        final String baseString = new String(referenceBases.getBases());
        if (baseString.contains("N")) return false;

        final String kmer = SVKmerizer.toKmer(baseString, new SVKmerShort(7)).canonical(7).toString(7);

        return kmers.contains(kmer);
    }
}
