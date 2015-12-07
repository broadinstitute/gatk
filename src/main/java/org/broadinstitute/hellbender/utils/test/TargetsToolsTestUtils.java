package org.broadinstitute.hellbender.utils.test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;

/**
 * Some common elements for target collection analysis tool unit and integration tests.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TargetsToolsTestUtils {

    private TargetsToolsTestUtils(){} //don't instantiate this

    /**
     * Returns a {@link File} pointing to the directory that contains the test data.
     * @return never {@code null}.
     */
    public static File getTestDataDir(){
        return new File(CommandLineProgramTest.getTestDataDir(),"exome");
    }

    /**
     * {@link File} pointing to the test toy reference used in targets analysis tool tests.
     */
    public static final File REFERENCE_FILE = new File(getTestDataDir(),"test_reference.fasta");

    /**
     * Sequence dictionary extracted from {@link #REFERENCE_FILE}.
     */
    public static final SAMSequenceDictionary REFERENCE_DICTIONARY = SAMSequenceDictionaryExtractor.extractDictionary(REFERENCE_FILE);

    /**
     * Creates a {@link SimpleInterval} instance given its contig and base range.
     * @param contig the new location contig name.
     * @param start  the new location start base index.
     * @param stop the new location stop base index.
     * @return never {@code null}.
     * @throws UserException if there was some problem when creating the location.
     */
    public static SimpleInterval createInterval(final String contig, final int start, final int stop) {
        return new SimpleInterval(REFERENCE_DICTIONARY.getSequence(contig).getSequenceName(),start,stop);
    }


    /**
     * Creates a {@link SimpleInterval} instance on an entire contig.
     * @param contigIndex the new location contig index.
     * @return never {@code null}.
     * @throws UserException if there was some problem when creating the location.
     */
    public static SimpleInterval createOverEntireContig(final int contigIndex) {
        final int contigLength = REFERENCE_DICTIONARY.getSequence(contigIndex).getSequenceLength();
        return new SimpleInterval(REFERENCE_DICTIONARY.getSequence(contigIndex).getSequenceName(),1,contigLength);
    }

    /**
     * Creates a {@link SimpleInterval} instance on an entire contig.
     * @param contig the new location contig.
     * @return never {@code null}.
     * @throws UserException if there was some problem when creating the location.
     */
    public static SimpleInterval createOverEntireContig(final String contig) {
        final int contigLength = REFERENCE_DICTIONARY.getSequence(contig).getSequenceLength();
        return new SimpleInterval(REFERENCE_DICTIONARY.getSequence(contig).getSequenceName(),1,contigLength);
    }

    /**
     * Creates a {@link SimpleInterval} at a give contig and position.
     * @param contig the contig name.
     * @param start the start and stop position.
     * @return never {@code null}.
     * @throws UserException if there was some problem when creating the location.
     */
    public static SimpleInterval createInterval(final String contig, final int start) {
        return new SimpleInterval(contig,start,start);
    }
}
