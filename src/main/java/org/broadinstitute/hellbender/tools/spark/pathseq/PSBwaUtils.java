package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;

/**
 * Utility functions for PathSeq Bwa tool
 */
public final class PSBwaUtils {

    static void addReferenceSequencesToHeader(final SAMFileHeader header,
                                              final String referencePath,
                                              final SerializableFunction<GATKRead, SimpleInterval> windowFunction,
                                              final PipelineOptions options) {
        final List<SAMSequenceRecord> refSeqs = getReferenceSequences(referencePath, windowFunction, options);
        for (final SAMSequenceRecord rec : refSeqs) {
            if (header.getSequence(rec.getSequenceName()) == null) {
                header.addSequence(rec);
            }
        }
    }

    private static List<SAMSequenceRecord> getReferenceSequences(final String referencePath,
                                                                 final SerializableFunction<GATKRead, SimpleInterval> windowFunction,
                                                                 final PipelineOptions options) {
        final ReferenceMultiSource referenceSource = new ReferenceMultiSource(options, referencePath, windowFunction);
        final SAMSequenceDictionary referenceDictionary = referenceSource.getReferenceSequenceDictionary(null);
        if (referenceDictionary == null) {
            throw new UserException.MissingReferenceDictFile(referencePath);
        }
        return referenceDictionary.getSequences();
    }

    /**
     * Returns list of sequence names that were aligned to at least once in the reads
     */
    static List<String> getAlignedSequenceNames(final JavaRDD<GATKRead> reads) {
        return reads.flatMap(read -> {
            if (read.isUnmapped() || read.getAssignedContig().equals("*")) Collections.emptySet().iterator();
            final Set<String> sequenceNames;
            if (read.hasAttribute("SA")) {
                final String[] saTokens = read.getAttributeAsString("SA").split(";");
                sequenceNames = new HashSet<>(1 + saTokens.length);
                sequenceNames.add(read.getAssignedContig());
                for (final String token : saTokens) {
                    final String[] alignmentTokens = token.split(",", 1);
                    sequenceNames.add(alignmentTokens[0]);
                }
            } else {
                sequenceNames = Collections.singleton(read.getAssignedContig());
            }
            return sequenceNames.iterator();
        }).filter(Objects::nonNull).distinct().collect();
    }

}