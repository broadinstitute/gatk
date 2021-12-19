package org.broadinstitute.hellbender.testutils;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;

import java.nio.file.Path;

/**
 * Utilities for comparing Fasta files.
 */
public final class FastaTestUtils {

    private FastaTestUtils() {}

    /**
     * Assert that two fasta files contain the same sequences, with the same names, in the same order.
     * This is different from exact file equality because they may have different sequence descriptions and line lengths.
     */
    public static void assertFastaFilesContainTheSameSequence(Path actualFasta, Path expectedFasta){
        try(final CachingIndexedFastaSequenceFile actual = new CachingIndexedFastaSequenceFile(actualFasta, 100000, true, true);
            final CachingIndexedFastaSequenceFile expected = new CachingIndexedFastaSequenceFile(expectedFasta,10000, true, true))
        {
            actual.getSequenceDictionary().assertSameDictionary(expected.getSequenceDictionary());
            for( SAMSequenceRecord record : actual.getSequenceDictionary().getSequences()){
                final byte[] actualBases = actual.getSequence(record.getSequenceName()).getBases();
                final byte[] expectedBases = expected.getSequence(record.getSequenceName()).getBases();
                Assert.assertEquals(actualBases, expectedBases, "Bases in reference were different.\n" +
                       new String(actualBases) + "\n" +new String(expectedBases));
            }
        }
    }

    public static void assertFastaFilesContainTheSameSequenceCaseInsensitive(Path actualFasta, Path expectedFasta){
        try(final CachingIndexedFastaSequenceFile actual = new CachingIndexedFastaSequenceFile(actualFasta, 100000, false, true);
            final CachingIndexedFastaSequenceFile expected = new CachingIndexedFastaSequenceFile(expectedFasta,10000, false, true))
        {
            actual.getSequenceDictionary().assertSameDictionary(expected.getSequenceDictionary());
            for( SAMSequenceRecord record : actual.getSequenceDictionary().getSequences()){
                final byte[] actualBases = actual.getSequence(record.getSequenceName()).getBases();
                final byte[] expectedBases = expected.getSequence(record.getSequenceName()).getBases();
                Assert.assertEquals(actualBases, expectedBases, "Bases in reference were different.\n" +
                        new String(actualBases) + "\n" +new String(expectedBases));
            }
        }
    }
}
