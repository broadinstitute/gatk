package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Class to load a reference sequence from a fasta file on Spark.
 */
public class ReferenceFileSparkSource implements ReferenceSparkSource, Serializable {
    private static final long serialVersionUID = 1L;

    private final URI referenceUri;
    // ReferenceFileSource is serializable (for running across a cluster), so we re-make the path
    // if necessary. Of course if you give us a Path to a RAM filesystem, don't expect
    // it to survive serialization.
    private transient Path referencePath;

    /**
     * @param referenceSpecifier the path to the reference file
     */
    public ReferenceFileSparkSource( final GATKPath referenceSpecifier) {
        // It would simplify this class if we could cache the GATKPath, but ReferenceFileSparkSource
        // objects are used as Spark broadcast variables, and caching GATKPath here triggers a known
        // issue during broadcast with the Java 11 GATK build. See https://issues.apache.org/jira/browse/SPARK-26963.
        referencePath = referenceSpecifier.toPath();
        referenceUri = referencePath.toUri();
        if (!Files.exists(referencePath)) {
            throw new UserException.MissingReference("The specified fasta file (" + referenceSpecifier.getRawInputString() + ") does not exist.");
        }
    }

    private synchronized Path getReferencePath() {
        if (null == referencePath) {
            this.referencePath = (new GATKPath(referenceUri.toString()).toPath());
        }
        return referencePath;
    }

    @Override
    public ReferenceBases getReferenceBases(final SimpleInterval interval) throws IOException {
        try ( ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(getReferencePath()) ) {
            ReferenceSequence sequence = referenceSequenceFile.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd());
            return new ReferenceBases(sequence.getBases(), interval);
        }
    }

    public Map<String, ReferenceBases> getAllReferenceBases() throws IOException {
        try ( ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(getReferencePath()) ) {
            Map<String, ReferenceBases> bases = new LinkedHashMap<>();
            ReferenceSequence seq;
            while ( (seq = referenceSequenceFile.nextSequence()) != null ) {
                String name = seq.getName();
                bases.put(name, new ReferenceBases(seq.getBases(), new SimpleInterval(name, 1, seq.length())));
            }
            return bases;
        }
    }

    @Override
    public SAMSequenceDictionary getReferenceSequenceDictionary(final SAMSequenceDictionary optReadSequenceDictionaryToMatch) throws IOException {
        try ( ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(getReferencePath()) ) {
            return referenceSequenceFile.getSequenceDictionary();
        }
    }
}
