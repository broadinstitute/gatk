package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;

/**
 * Some common elements for exome analysis tool unit and integration tests.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ExomeToolsTestUtils {

    /**
     * Returns a {@link File} pointing to the directory that contains the test data.
     * @return never {@code null}.
     */
    protected static File getTestDataDir(){
        return new File(CommandLineProgramTest.getTestDataDir(),"exome");
    }

    /**
     * {@link File} pointing to the test toy reference used in exome analysis tool tests.
     */
    protected final static File REFERENCE_FILE = new File(getTestDataDir(),"test_reference.fasta");

    /**
     * Sequence dictionary extracted from {@link #REFERENCE_FILE}.
     */
    protected final static SAMSequenceDictionary REFERENCE_DICTIONARY = SAMSequenceDictionaryExtractor.extractDictionary(REFERENCE_FILE);

    /**
     * {@link GenomeLocParser} instance to use for creating genome locations.
     */
    protected final static GenomeLocParser GENOME_LOC_FACTORY = new GenomeLocParser(REFERENCE_DICTIONARY);


    /**
     * Creates a {@link GenomeLoc} instance given its contig and base range.
     * @param contig the new location contig.
     * @param start  the new location start base index.
     * @param stop the new location stop base index.
     * @return never {@code null}.
     * @throws UserException if there was some problem when creating the location.
     */
    protected static GenomeLoc createGenomeLoc(final String contig, final int start, final int stop) {
        return GENOME_LOC_FACTORY.createGenomeLoc(contig,start,stop);
    }

    /**
     * Creates a {@link GenomeLoc} instance on an entire contig.
     * @param contig the new location contig.
     * @return never {@code null}.
     * @throws UserException if there was some problem when creating the location.
     */
    protected static GenomeLoc createOverEntireContig(final String contig) {
        return GENOME_LOC_FACTORY.createOverEntireContig(contig);
    }

    /**
     * Creates a {@link GenomeLoc} at a give contig and position.
     * @param contig the contig name.
     * @param start the start and stop position.
     * @return never {@code null}.
     * @throws UserException if there was some problem when creating the location.
     */
    protected static GenomeLoc createGenomeLoc(final String contig, final int start) {
        return GENOME_LOC_FACTORY.createGenomeLoc(contig,start);
    }
}
