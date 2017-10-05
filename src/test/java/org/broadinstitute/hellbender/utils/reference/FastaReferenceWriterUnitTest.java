package org.broadinstitute.hellbender.utils.reference;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.nio.file.Path;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Unit tests for {@link FastaReferenceWriter}.
 */
public class FastaReferenceWriterUnitTest extends GATKBaseTest {

    @Test(expectedExceptions = IllegalStateException.class)
    public void testEmptySequence() throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        try (final FastaReferenceWriter writer = new FastaReferenceWriter(testOutput.toPath(), false, false)) {
            writer.startSequence("seq1");
            writer.appendBases(new RandomDNA(113).nextBases(100));
            writer.startSequence("seq2");
            writer.startSequence("seq3");
        } finally {
            Assert.assertTrue(testOutput.delete());
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testEmptyReference() throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        try (final FastaReferenceWriter writer = new FastaReferenceWriter(testOutput.toPath(), false, false)) {
        } finally {
            Assert.assertTrue(testOutput.delete());
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testStartSequenceAfterClose() throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        final FastaReferenceWriter writer = new FastaReferenceWriter(testOutput.toPath(), false, false);
        writer.startSequence("seq1").appendBases(new byte[] { 'A', 'C', 'G', 'T'});
        writer.close();
        try {
            writer.startSequence("seq2");
        } finally {
            Assert.assertTrue(testOutput.delete());
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testAddBasesAfterClose() throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        final FastaReferenceWriter writer = new FastaReferenceWriter(testOutput.toPath(), false, false);
        writer.startSequence("seq1").appendBases(new byte[] { 'A', 'C', 'G', 'T'});
        writer.close();
        try {
            writer.appendBases(new byte[]{'A', 'A', 'A'});
        } finally {
            Assert.assertTrue(testOutput.delete());
        }
    }

    @Test(dataProvider = "invalidBplData", expectedExceptions = IllegalArgumentException.class)
    public void testBadDefaultBasesPerLine(final int invalidBpl) throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        try {
            new FastaReferenceWriter(testOutput.toPath(), invalidBpl, false, false);
        } finally {
            @SuppressWarnings("unused")  // it might create it or it might not.
            final boolean dummy = testOutput.delete();
        }
    }

    @Test(dataProvider = "invalidBplData", expectedExceptions = IllegalArgumentException.class)
    public void testBadSequenceBasesPerLine(final int invalidBpl) throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        try (final FastaReferenceWriter writer = new FastaReferenceWriter(testOutput.toPath(), false, false)) {
            writer.startSequence("seq1", invalidBpl);
        } finally {
            Assert.assertTrue(testOutput.delete());
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testEmptySequenceAtTheEnd() throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        try (final FastaReferenceWriter writer = new FastaReferenceWriter(testOutput.toPath(), false, false)) {
            writer.startSequence("seq1");
            writer.appendBases(new RandomDNA(113).nextBases(100));
            writer.startSequence("seq2");
            writer.appendBases(new RandomDNA(13).nextBases(1001));
            writer.startSequence("seq3");
        } finally {
            Assert.assertTrue(testOutput.delete());
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testAppendBasesBeforeStartingSequence() throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        try (final FastaReferenceWriter writer = new FastaReferenceWriter(testOutput.toPath(), false, false)) {
            writer.appendBases(new RandomDNA(113).nextBases(100));
        } finally {
            Assert.assertTrue(testOutput.delete());
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testAddingSameSequenceTwice() throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        try (final FastaReferenceWriter writer = new FastaReferenceWriter(testOutput.toPath(), false, false)) {
            writer.startSequence("seq1");
            writer.appendBases(new RandomDNA(113).nextBases(100));
            writer.startSequence("seq2");
            writer.appendBases(new RandomDNA(114).nextBases(300));
            writer.startSequence("seq1");
        } finally {
            Assert.assertTrue(testOutput.delete());
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testAddingSameSequenceRightAfter() throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        try (final FastaReferenceWriter writer = new FastaReferenceWriter(testOutput.toPath(), false, false)) {
            writer.startSequence("seq1");
            writer.appendBases(new RandomDNA(113).nextBases(100));
            writer.startSequence("seq1");
        } finally {
            Assert.assertTrue(testOutput.delete());
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class, dataProvider = "invalidNameData")
    public void testAddingInvalidSequenceName(final String invalidName) throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        try (final FastaReferenceWriter writer = new FastaReferenceWriter(testOutput.toPath(), false, false)) {
            writer.startSequence(invalidName);
        } finally {
            Assert.assertTrue(testOutput.delete());
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class, dataProvider = "invalidDescriptionData")
    public void testAddingInvalidDescription(final String invalidDescription) throws IOException {
        final File testOutput = createTempFile("fwr-test", ".fasta");
        Assert.assertTrue(testOutput.delete());
        try (final FastaReferenceWriter writer = new FastaReferenceWriter(testOutput.toPath(), false, false)) {
            writer.startSequence("seq1", invalidDescription);
        } finally {
            Assert.assertTrue(testOutput.delete());
        }
    }

    @DataProvider(name = "invalidBplData")
    public Object[][] invalidBplData() {
        return IntStream.of(0, -1, -110)
                .mapToObj(i -> new Object[] { i }).toArray(Object[][]::new);
    }

    @DataProvider(name = "invalidNameData")
    public Object[][] invalidNameData() {
        return Stream.of("seq with spaces", "seq\twith\ttabs", "with blank", " ", "", "nnn\n", "rrr\r", null)
                .map(s -> new Object[] {s}).toArray(Object[][]::new);
    }

    @DataProvider(name = "invalidDescriptionData")
    public Object[][] invalidDescriptionData() {
        return Stream.of("\nwith control chars\nthat are not\0tabs\r", "with the null\0", "with nl\n")
                .map(s -> new Object[] {s}).toArray(Object[][]::new);
    }



    @Test(dataProvider = "testData")
    public void testWriter(final SAMSequenceDictionary dictionary, final boolean withIndex, final boolean withDictionary,
                           final boolean withDescriptions, final int defaultBpl,
                                        final int minBpl, final int maxBpl, final int seed)
            throws IOException, GeneralSecurityException, URISyntaxException {
        final Map<String, byte[]> bases = new LinkedHashMap<>(dictionary.getSequences().size());
        final Map<String, Integer> bpl = new LinkedHashMap<>(dictionary.getSequences().size());
        final Random rdn = new Random(seed);
        generateRandomBasesAndBpls(dictionary, minBpl, maxBpl, bases, bpl, rdn);
        final File fastaFile = createTempFile("fwr-test", ".fa");
        Assert.assertTrue(fastaFile.delete());
        final File fastaIndexFile = new File(fastaFile.getParentFile(), fastaFile.getName() + ".fai");
        final File dictFile = new File(fastaFile.getParentFile(), fastaFile.getName().replaceAll("\\.fa", ".dict"));
        fastaIndexFile.deleteOnExit();
        dictFile.deleteOnExit();

        try (final FastaReferenceWriter writer = defaultBpl < 0
                ? new FastaReferenceWriter(fastaFile.toPath(), withIndex, withDictionary)
                : new FastaReferenceWriter(fastaFile.toPath(), defaultBpl, withIndex, withDictionary)) {
            writeReference(writer, withDescriptions, rdn, dictionary, bases, bpl);
        }
        assertOutput(fastaFile.toPath(), withIndex, withDictionary, withDescriptions, dictionary, defaultBpl, bases, bpl);
        Assert.assertTrue(fastaFile.delete());
        Assert.assertEquals(fastaIndexFile.delete(), withIndex);
        Assert.assertEquals(dictFile.delete(), withDictionary);
    }

    private void generateRandomBasesAndBpls(SAMSequenceDictionary dictionary, int minBpl, int maxBpl, Map<String, byte[]> bases, Map<String, Integer> bpl, Random rdn) {
        final RandomDNA rdnDNA = new RandomDNA(rdn.nextLong());
        // We avoid to use the obvious first choice {@link RandomDNA#nextFasta} as these may actually use
        // this writer to do its job eventually.
        for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
            bases.put(sequence.getSequenceName(), rdnDNA.nextBases(sequence.getSequenceLength()));
            if (rdn.nextDouble() < 0.333333) { // 1/3 of the time we will use the default.
                bpl.put(sequence.getSequenceName(), -1);
            } else {
                bpl.put(sequence.getSequenceName(), rdn.nextInt(maxBpl - minBpl + 1) + minBpl);
            }
        }
    }

    private void writeReference(final FastaReferenceWriter writer, final boolean withDescriptions,
                               final Random rdn, final SAMSequenceDictionary dictionary,
                               final Map<String, byte[]> seqs,
                               final Map<String, Integer> basesPerLine)
        throws IOException
    {
        for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
            final int bpl = basesPerLine.get(sequence.getSequenceName());
            if (withDescriptions) {
               final String description = String.format("index=%d\tlength=%d",
                            dictionary.getSequenceIndex(sequence.getSequenceName()),
                            sequence.getSequenceLength());
               if (bpl < 0) {
                   writer.startSequence(sequence.getSequenceName(), description);
               } else {
                   writer.startSequence(sequence.getSequenceName(), description, bpl);
               }
            } else {
                if (bpl < 0) {
                    writer.startSequence(sequence.getSequenceName());
                } else {
                    writer.startSequence(sequence.getSequenceName(), bpl);
                }
            }
            final boolean onOneGo = rdn.nextDouble() < 0.25; // 25% of times we just write the whole sequence of one go.
            if (onOneGo) {
                writer.appendBases(seqs.get(sequence.getSequenceName()));
            } else {
                int done = 0;
                while (done < seqs.get(sequence.getSequenceName()).length) {
                    final boolean useBpl = bpl > 0 && rdn.nextDouble() < 0.10; // 10% of times we exactly add the same number of bases as bases-per-line, remaining bases permitting.
                    int left = sequence.getSequenceLength() - done;
                    final int length = useBpl ? Math.min(bpl, left) : rdn.nextInt(left) + 1;
                    Assert.assertSame(writer.appendBases(seqs.get(sequence.getSequenceName()), done, length), writer);
                    done += length;
                    left -= length;
                    if (useBpl && rdn.nextDouble() < 0.10) { // 10% of the time we align with bpl so that it will start a new line on the next write.
                        final int lengthToEndOfLine = Math.min(left, bpl - (done % bpl));
                        Assert.assertSame(writer.appendBases(seqs.get(sequence.getSequenceName()), done, lengthToEndOfLine), writer);
                        done += lengthToEndOfLine;
                    }
                    if (rdn.nextDouble() < 0.10) { // 10% of the time we do a stupid zero length append.
                        Assert.assertSame(writer.appendBases(seqs.get(sequence.getSequenceName()), done, 0), writer);
                    }
                }
            }
        }
    }

    private void assertOutput(final Path path, final boolean mustHaveIndex, final boolean mustHaveDictionary,
                              final boolean withDescriptions, final SAMSequenceDictionary dictionary, final int defaultBpl,
                              final Map<String, byte[]> bases, final Map<String, Integer> basesPerLine)
            throws GeneralSecurityException, IOException, URISyntaxException {
        assertFastaContent(path, withDescriptions, dictionary, defaultBpl, bases, basesPerLine);
        if (mustHaveDictionary) {
            assertFastaDictionaryContent(path, dictionary);
        }
        if (mustHaveIndex) {
            assertFastaIndexContent(path, dictionary, bases, basesPerLine);
        }
    }


    private void assertFastaContent(final Path path, final boolean withDescriptions, final SAMSequenceDictionary dictionary, final int defaultBpl,
                                    final Map<String, byte[]> bases, final Map<String, Integer> basesPerLine)
        throws IOException
    {
        try (final BufferedReader reader = new BufferedReader(new InputStreamReader(path.getFileSystem().provider().newInputStream(path)))) {
            for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
                final String description = String.format("index=%d\tlength=%d",
                        dictionary.getSequenceIndex(sequence.getSequenceName()), sequence.getSequenceLength());
                final String expectedHeader =
                        FastaReferenceWriter.HEADER_START_CHAR + sequence.getSequenceName()
                        + ((withDescriptions) ? FastaReferenceWriter.HEADER_NAME_AND_DESCRIPTION_SEPARATOR + description : "");
                Assert.assertEquals(reader.readLine(), expectedHeader);
                final byte[] expectedBases = bases.get(sequence.getSequenceName());
                final int bpl_ = basesPerLine.get(sequence.getSequenceName());
                final int bpl = bpl_ < 0 ? (defaultBpl < 0 ? FastaReferenceWriter.DEFAULT_BASES_PER_LINE : defaultBpl) : bpl_;
                int offset = 0;
                while (offset < expectedBases.length) {
                    final int expectedLength = Math.min(expectedBases.length - offset, bpl);
                    final byte[] expectedBaseLine = SequenceUtil.upperCase(Arrays.copyOfRange(expectedBases, offset , offset + expectedLength));
                    final byte[] actualBaseLine = SequenceUtil.upperCase(reader.readLine().getBytes());
                    Assert.assertEquals(actualBaseLine, expectedBaseLine);
                    offset += expectedLength;
                }
            }
        }
    }

    private void assertFastaIndexContent(final Path path, final SAMSequenceDictionary dictionary,
                                    final Map<String, byte[]> bases, final Map<String, Integer> basesPerLine)
        throws IOException
    {
        final IndexedFastaSequenceFile indexedFasta = new IndexedFastaSequenceFile(path);
        for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
            final String name = sequence.getSequenceName();
            final int length = sequence.getSequenceLength();
            final ReferenceSequence start = indexedFasta.getSubsequenceAt(name, 1, Math.min(length, 30));
            final ReferenceSequence end = indexedFasta.getSubsequenceAt(name, Math.max(1, length - 29), length);
            final int middlePos = Math.max(1, Math.min(length, length / 2));
            final ReferenceSequence middle = indexedFasta.getSubsequenceAt(name, middlePos, Math.min(middlePos + 29, length));
            Assert.assertEquals(start.getBases(), Arrays.copyOfRange(bases.get(name), 0, start.length()));
            Assert.assertEquals(end.getBases(), Arrays.copyOfRange(bases.get(name), Math.max(0, length - 30), length));
            Assert.assertEquals(middle.getBases(), Arrays.copyOfRange(bases.get(name), middlePos - 1, middlePos -1 + middle.length()));
        }
    }

    private void assertFastaDictionaryContent(final Path path, final SAMSequenceDictionary dictionary)
            throws IOException, GeneralSecurityException, URISyntaxException {
        final Path dictPath = path.resolveSibling(path.getFileName().toString().replaceAll("\\.(fa|fasta)",".dict"));
        final ReadsDataSource readsDataSource = new ReadsDataSource(dictPath);
        final SAMFileHeader actualHeader = readsDataSource.getHeader();
        final SAMSequenceDictionary actualDictionary = actualHeader.getSequenceDictionary();
        dictionary.assertSameDictionary(actualDictionary);
    }

    @DataProvider(name="testData")
    public Object[][] testData() {
        // data-signature: (SAMSequenceDictionary dictionary, boolean withDescriptions, int defaultBpl, int minBpl, int maxBpl, int seed
        // defaultBpl == -1 means to use the default {@link FastaReferenceWriter#DEFAULT_BASE_PER_LINE}.
        // [minBpl , manBpl] range for possible bpl when the default for the file is not to be used.
        final Random rdn = new Random(113);
        final SAMSequenceDictionary typicalDictionary = new SAMSequenceDictionary(
                Arrays.asList(new SAMSequenceRecord("chr1", 10_000),
                              new SAMSequenceRecord("chr2", 20_000),
                              new SAMSequenceRecord("chr3", 20_000),
                              new SAMSequenceRecord("chr4", 2_000),
                              new SAMSequenceRecord("chr5", 200),
                              new SAMSequenceRecord("chr6", 3_010),
                              new SAMSequenceRecord("X", 1_000)
        ));

        final SAMSequenceDictionary manyBPLMatchingSequences = new SAMSequenceDictionary(
                IntStream.range(0, 100)
                        .mapToObj(i -> new SAMSequenceRecord("" + (i + 1), FastaReferenceWriter.DEFAULT_BASES_PER_LINE * (rdn.nextInt(10) + 1) ))
                        .collect(Collectors.toList())
        );

        final SAMSequenceDictionary singleSequence = new SAMSequenceDictionary(Collections.singletonList(new SAMSequenceRecord("seq", 2_000)));

        final SAMSequenceDictionary oneBaseSequencesContaining = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("chr1", 10_000),
                new SAMSequenceRecord("chr2", 20_000),
                new SAMSequenceRecord("chr2.small", 1),
                new SAMSequenceRecord("X", 1_000),
                new SAMSequenceRecord("MT", 1))
        );

        final SAMSequenceDictionary[] testDictionaries = new SAMSequenceDictionary[] {typicalDictionary, manyBPLMatchingSequences, singleSequence, oneBaseSequencesContaining};
        final int[] testBpls = new int[] { -1 , FastaReferenceWriter.DEFAULT_BASES_PER_LINE, 1, 100, 51, 63};
        final boolean[] testWithDescriptions = new boolean[] { true, false};
        final boolean[] testWithIndex = new boolean[] { true, false};
        final boolean[] testWithDictionary = new boolean[] { true, false};
        final int[] testSeeds = new int [] { 31, 113, 73 };
        final List<Object[]> result = new ArrayList<>();
        for (final SAMSequenceDictionary dictionary : testDictionaries) {
            for (final boolean withIndex : testWithIndex) {
                for (final boolean withDictionary : testWithDictionary) {
                    for (final boolean withDescriptions : testWithDescriptions) {
                        for (final int bpl : testBpls) {
                            for (final int seed : testSeeds) {
                                result.add(new Object[]{dictionary, withIndex, withDictionary, withDescriptions, bpl, 1, (bpl < 0 ? FastaReferenceWriter.DEFAULT_BASES_PER_LINE : bpl) * 2, seed});
                            }
                        }
                    }
                }
            }
        }
        return result.stream().toArray(Object[][]::new);
    }

}
