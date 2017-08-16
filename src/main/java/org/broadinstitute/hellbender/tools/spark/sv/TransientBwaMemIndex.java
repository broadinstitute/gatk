package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;

/**
 * Wraps a {@link BwaMemIndex} that keeps track of the actual file containg the index in disk and
 * deletes it when closed.
 */
public final class TransientBwaMemIndex implements AutoCloseable {

    private final BwaMemIndex memIndex;
    private final File fileName;
    private boolean closed = false;

    public TransientBwaMemIndex(final String name, final byte[] bases) {
        try {
            final File fasta = File.createTempFile(name, ".fasta");
            fileName = File.createTempFile(name, ".img");
            FastaReferenceWriter.writeSingleSequenceReference(fasta.toPath(), false, false, "seq1", "", bases);
            BwaMemIndex.createIndexImageFromFastaFile(fasta.toString(), fileName.toString());
            if (!fasta.delete()) {
                LogManager.getLogger(GenotypeStructuralVariantsSpark.class)
                        .warn("Could not remove temporary fasta file: " + fasta.toString());
            }
            fileName.deleteOnExit();
            memIndex = new BwaMemIndex(fileName.toString());
        } catch (final IOException ex) {
            throw new UncheckedIOException(ex);
        }
    }

    public BwaMemIndex get() {
        if (closed) {
            throw new IllegalStateException("Cannot get a closed index");
        }
        return memIndex;
    }

    @Override
    public void close()  {
        if (!closed) {
            memIndex.close();
            closed = true;
            if (!fileName.delete() && fileName.exists()) {
                    throw new IllegalStateException("deleting the index file didn't work for " + fileName);
            }
        }
    }
}
