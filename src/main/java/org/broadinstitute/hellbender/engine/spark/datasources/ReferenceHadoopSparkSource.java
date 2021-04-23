package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.Serializable;
import java.net.URI;

/**
 * Class to load a reference sequence from a fasta file on HDFS.
 */
public class ReferenceHadoopSparkSource implements ReferenceSparkSource, Serializable {
    private static final long serialVersionUID = 1L;

    private final URI referenceURI;

    /**
     * @param referencePathSpecifier the path to the reference file on HDFS
     */
    public ReferenceHadoopSparkSource( final GATKPath referencePathSpecifier) {
        // It would simplify this class if we could cache the GATKPath, but ReferenceFileSparkSource
        // objects are used as Spark broadcast variables, and caching GATKPath here triggers a known
        // issue during broadcast with the Java 11 GATK build. See https://issues.apache.org/jira/browse/SPARK-26963.
        this.referenceURI = referencePathSpecifier.getURI();
    }

    @Override
    public ReferenceBases getReferenceBases(final SimpleInterval interval) {
        ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new GATKPath(referenceURI.toString()).toPath());
        ReferenceSequence sequence = referenceSequenceFile.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd());
        return new ReferenceBases(sequence.getBases(), interval);
    }

    @Override
    public SAMSequenceDictionary getReferenceSequenceDictionary(final SAMSequenceDictionary optReadSequenceDictionaryToMatch) {
        ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new GATKPath(referenceURI.toString()).toPath());
        return referenceSequenceFile.getSequenceDictionary();
    }

}
