package org.broadinstitute.hellbender.testutils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.reference.FastaReferenceWriter;
import org.testng.Assert;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Map;

public class FastaTestUtils {

    public static void assertOutput(final Path path, final boolean mustHaveIndex, final boolean mustHaveDictionary,
                                    final boolean withDescriptions, final SAMSequenceDictionary dictionary, final int defaultBpl,
                                    final Map<String, byte[]> bases, final Map<String, Integer> basesPerLine)
            throws IOException {
        assertFastaContent(path, withDescriptions, dictionary, defaultBpl, bases, basesPerLine);
        if (mustHaveDictionary) {
            assertFastaDictionaryContent(ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(path), dictionary);
        }
        if (mustHaveIndex) {
            assertFastaIndexContent(path, ReferenceSequenceFileFactory.getFastaIndexFileName(path), dictionary, bases);
        }
    }

    public static void assertFastaContent(final Path path, final boolean withDescriptions, final SAMSequenceDictionary dictionary, final int defaultBpl,
                                          final Map<String, byte[]> bases, final Map<String, Integer> basesPerLine)
            throws IOException {
        try (final BufferedReader reader = new BufferedReader(new InputStreamReader(path.getFileSystem().provider().newInputStream(path)))) {
            for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
                final String description = String.format("index=%d\tlength=%d",
                        dictionary.getSequenceIndex(sequence.getSequenceName()), sequence.getSequenceLength());
                final String expectedHeader = FastaReferenceWriter.HEADER_START_CHAR + sequence.getSequenceName()
                                + ((withDescriptions) ? FastaReferenceWriter.HEADER_NAME_AND_DESCRIPTION_SEPARATOR + description : "");
                Assert.assertTrue(reader.readLine().startsWith(expectedHeader));
                final byte[] expectedBases = bases.get(sequence.getSequenceName());
                final int bpl_ = basesPerLine.get(sequence.getSequenceName());
                final int bpl = bpl_ < 0 ? (defaultBpl < 0 ? FastaReferenceWriter.DEFAULT_BASES_PER_LINE : defaultBpl) : bpl_;
                int offset = 0;
                while (offset < expectedBases.length) {
                    final int expectedLength = Math.min(expectedBases.length - offset, bpl);
                    final byte[] expectedBaseLine = SequenceUtil.upperCase(
                            Arrays.copyOfRange(expectedBases, offset, offset + expectedLength));
                    final byte[] actualBaseLine = SequenceUtil.upperCase(reader.readLine().getBytes());
                    Assert.assertEquals(actualBaseLine, expectedBaseLine);
                    offset += expectedLength;
                }
            }
        }
    }

    public static void assertFastaIndexContent(final Path path, final Path indexPath, final SAMSequenceDictionary dictionary,
                                               final Map<String, byte[]> bases) {
        final FastaSequenceIndex index = new FastaSequenceIndex(indexPath);
        final IndexedFastaSequenceFile indexedFasta = new IndexedFastaSequenceFile(path, index);
        for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
            final String name = sequence.getSequenceName();
            final int length = sequence.getSequenceLength();
            final ReferenceSequence start = indexedFasta.getSubsequenceAt(name, 1, Math.min(length, 30));
            final ReferenceSequence end = indexedFasta.getSubsequenceAt(name, Math.max(1, length - 29), length);
            final int middlePos = Math.max(1, Math.min(length, length / 2));
            final ReferenceSequence middle = indexedFasta.getSubsequenceAt(name, middlePos, Math.min(middlePos + 29, length));
            Assert.assertEquals(start.getBases(), Arrays.copyOfRange(bases.get(name), 0, start.length()));
            Assert.assertEquals(end.getBases(), Arrays.copyOfRange(bases.get(name), Math.max(0, length - 30), length));
            Assert.assertEquals(middle.getBases(), Arrays.copyOfRange(bases.get(name), middlePos - 1, middlePos - 1 + middle.length()));
        }
    }

    public static void assertFastaDictionaryContent(final Path dictPath, final SAMSequenceDictionary dictionary) {
        final ReadsDataSource readsDataSource = new ReadsDataSource(dictPath);
        final SAMFileHeader actualHeader = readsDataSource.getHeader();
        final SAMSequenceDictionary actualDictionary = actualHeader.getSequenceDictionary();
        dictionary.assertSameDictionary(actualDictionary);
    }

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
}
