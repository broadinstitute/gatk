package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.spark.datasources.ReferenceMultiSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.apache.logging.log4j.Logger;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Utility functions for PathSeq Bwa tool
 */
public final class PSBwaUtils {

    static void addReferenceSequencesToHeader(final SAMFileHeader header,
                                              final String referencePath,
                                              final SerializableFunction<GATKRead, SimpleInterval> windowFunction) {
        final List<SAMSequenceRecord> refSeqs = getReferenceSequences(referencePath, windowFunction);
        for (final SAMSequenceRecord rec : refSeqs) {
            if (header.getSequence(rec.getSequenceName()) == null) {
                header.addSequence(rec);
            }
        }
    }

    private static List<SAMSequenceRecord> getReferenceSequences(final String referencePath,
                                                                 final SerializableFunction<GATKRead, SimpleInterval> windowFunction) {
        final ReferenceMultiSparkSource referenceSource = new ReferenceMultiSparkSource(referencePath, windowFunction);
        final SAMSequenceDictionary referenceDictionary = referenceSource.getReferenceSequenceDictionary(null);
        if (referenceDictionary == null) {
            throw new UserException.MissingReferenceDictFile(referencePath);
        }
        return referenceDictionary.getSequences();
    }

    /**
     * Returns header with sequences that were aligned to at least once in reads
     */
    static SAMFileHeader removeUnmappedHeaderSequences(final SAMFileHeader header,
                                                       final JavaRDD<GATKRead> reads,
                                                       final Logger logger) {
        final List<String> usedSequences = PSBwaUtils.getAlignedSequenceNames(reads);
        final List<SAMSequenceRecord> usedSequenceRecords = usedSequences.stream()
                .map(seqName -> header.getSequence(seqName))
                .filter(seq -> {
                    if (seq != null) return true;
                    logger.warn("One or more reads are aligned to sequence " + seq + " but it is not in the header");
                    return false;
                })
                .collect(Collectors.toList());
        header.setSequenceDictionary(new SAMSequenceDictionary(usedSequenceRecords));
        return header;
    }

    /**
     * Returns list of sequence names that were aligned to at least once in the reads
     */
    static List<String> getAlignedSequenceNames(final JavaRDD<GATKRead> reads) {
        return reads.flatMap(PSBwaUtils::getSequenceNames).distinct().collect();
    }

    /**
     * Returns set of sequence names of the read
     */
    private static Iterator<String> getSequenceNames(final GATKRead read) {
        if (read.isUnmapped() || read.getAssignedContig().equals("*")) return Collections.emptyIterator();
        if (!read.hasAttribute("SA")) return Collections.singleton(read.getAssignedContig()).iterator();
        final String[] saTokens = read.getAttributeAsString("SA").split(";");
        final Set<String> sequenceNames = new HashSet<>(SVUtils.hashMapCapacity(1 + saTokens.length));
        sequenceNames.add(read.getAssignedContig());
        for (final String token : saTokens) {
            final String[] alignmentTokens = token.split(",", 1);
            sequenceNames.add(alignmentTokens[0]);
        }
        return sequenceNames.iterator();
    }

}