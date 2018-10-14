package org.broadinstitute.hellbender.engine.spark.datasources;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
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
     * @param referenceUri the path to the reference file
     */
    public ReferenceFileSparkSource( final String referenceUri) {
        this(IOUtils.getPath(referenceUri));
    }

    /**
     * @param referencePath the path to the reference file
     */
    public ReferenceFileSparkSource( final Path referencePath) {
        this.referencePath = referencePath;
        this.referenceUri = referencePath.toUri();
        if (!Files.exists(this.referencePath)) {
            throw new UserException.MissingReference("The specified fasta file (" + referencePath.toAbsolutePath().toUri().toString() + ") does not exist.");
        }
    }

    private synchronized Path getReferencePath() {
        if (null == referencePath) {
            this.referencePath = Paths.get(referenceUri);
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
