package org.broadinstitute.hellbender.utils.diffengine;

import java.io.File;
import java.io.IOException;

/**
 * Interface for readers creating diffable objects from a file
 */
public interface DiffableReader {
    /**
     * Return the name of this DiffableReader type.  For example, the VCF reader returns 'VCF' and the
     * bam reader 'BAM'
     */
    public String getName();

    /**
     * Read up to maxElementsToRead DiffElements from file, and return them.
     * Note: If value -1 is used for maxElementsToRead, then all elements are to be read from the file.
     * @throws IOException if there's a problem reading from the file.
     */
    public DiffElement readFromFile(File file, int maxElementsToRead) throws IOException;

    /**
     * Return true if the file can be read into DiffElement objects with this reader. This should
     * be uniquely true/false for all readers, as the system will use the first reader that can read the
     * file.  This routine should never throw an exception.  The VCF reader, for example, looks at the
     * first line of the file for the ##format=VCF4.1 header, and the BAM reader for the BAM_MAGIC value
     * @param file
     * @return
     */
    public boolean canRead(File file);
}
