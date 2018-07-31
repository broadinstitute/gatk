package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.SamFiles;

import java.nio.file.Path;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


public class HaplotypeBAMWriterUnitTest extends GATKBaseTest {
    private final SAMFileHeader samHeader = ArtificialReadUtils.createArtificialSamHeader(20, 1, 1000);
    private final String expectedFilePath = getToolTestDataDir() + "/expected/";

    @Test
    public void testSimpleCreate() throws Exception {
        final MockValidatingDestination writer = new MockValidatingDestination(null);
        Assert.assertTrue(HaplotypeBAMWriter.create(HaplotypeBAMWriter.WriterType.CALLED_HAPLOTYPES, writer) instanceof CalledHaplotypeBAMWriter);
        Assert.assertTrue(HaplotypeBAMWriter.create(HaplotypeBAMWriter.WriterType.ALL_POSSIBLE_HAPLOTYPES, writer) instanceof AllHaplotypeBAMWriter);
    }

    @DataProvider(name="ReadsLikelikhoodData")
    public Object[][] makeReadsLikelikhoodData() {
        final String haplotypeBaseSignature = "ACTGAAGGTTCC";
        final Haplotype haplotype = makeHaplotype(haplotypeBaseSignature);
        final int sampleCounts[] = {2, 2};
        final Locatable loc = new SimpleInterval("20", 10, 20);
        final ReadLikelihoods<Haplotype> readLikelihoods = generateReadLikelihoods(sampleCounts);

        return new Object[][]{
                { haplotypeBaseSignature, Collections.singletonList(haplotype), loc, readLikelihoods }
        };
    }

    @Test(dataProvider = "ReadsLikelikhoodData")
    public void testWriteToSAMFile
        (
            @SuppressWarnings("unused") final String haplotypeBaseSignature,
            final List<Haplotype> haplotypes,
            final Locatable genomeLoc,
            final ReadLikelihoods <Haplotype> readLikelihoods
        ) throws IOException
    {
        // create OUTPUT SAM file
        final Path outPath = testWriteToFile(".sam", haplotypes, genomeLoc, readLikelihoods, false, false);
        final File expectedFile = new File(expectedFilePath, "testSAM.sam");
        IntegrationTestSpec.assertEqualTextFiles(outPath.toFile(), expectedFile);
    }

    @Test(dataProvider = "ReadsLikelikhoodData")
    public void testWriteToBAMFile
        (
            @SuppressWarnings("unused") final String haplotypeBaseSignature,
            final List<Haplotype> haplotypes,
            final Locatable genomeLoc,
            final ReadLikelihoods <Haplotype> readLikelihoods
        ) throws IOException
    {
        // create OUTPUT BAM file
        final Path outPath = testWriteToFile(".bam", haplotypes, genomeLoc, readLikelihoods, false, false);
        final File expectedFile = new File(expectedFilePath, "testBAM.bam");
        SamAssertionUtils.assertEqualBamFiles(outPath.toFile(), expectedFile, false, ValidationStringency.DEFAULT_STRINGENCY);
    }

    @Test(dataProvider = "ReadsLikelikhoodData")
    public void testWriteToBAMFileWithIndex
            (
                    @SuppressWarnings("unused") final String haplotypeBaseSignature,
                    final List<Haplotype> haplotypes,
                    final Locatable genomeLoc,
                    final ReadLikelihoods <Haplotype> readLikelihoods
            ) throws IOException
    {
        // create OUTPUT BAM file
        final Path outPath = testWriteToFile(".bam", haplotypes, genomeLoc, readLikelihoods, true, false);
        final File expectedFile = new File(expectedFilePath, "testBAM.bam");
        SamAssertionUtils.assertEqualBamFiles(outPath.toFile(), expectedFile, false, ValidationStringency.DEFAULT_STRINGENCY);
    }

    @Test(dataProvider = "ReadsLikelikhoodData")
    public void testWriteToBAMFileWithMD5
            (
                    @SuppressWarnings("unused") final String haplotypeBaseSignature,
                    final List<Haplotype> haplotypes,
                    final Locatable genomeLoc,
                    final ReadLikelihoods <Haplotype> readLikelihoods
            ) throws IOException
    {
        // create OUTPUT BAM file
        final Path outPath = testWriteToFile(".bam", haplotypes, genomeLoc, readLikelihoods, false, true);
        final File expectedFile = new File(expectedFilePath, "testBAM.bam");
        SamAssertionUtils.assertEqualBamFiles(outPath.toFile(), expectedFile, false, ValidationStringency.DEFAULT_STRINGENCY);
    }

    private Path testWriteToFile
        (
            final String outputFileExtension,
            final List<Haplotype> haplotypes,
            final Locatable genomeLoc,
            final ReadLikelihoods <Haplotype> readLikelihoods,
            final boolean createIndex,
            final boolean createMD5
        ) throws IOException
    {
        final Path outPath = GATKBaseTest.createTempFile("haplotypeBamWriterTest", outputFileExtension).toPath();
        final HaplotypeBAMDestination fileDest = new SAMFileDestination(outPath, createIndex, createMD5, samHeader, "TestHaplotypeRG");

        try (final HaplotypeBAMWriter haplotypeBAMWriter = HaplotypeBAMWriter.create(HaplotypeBAMWriter.WriterType.ALL_POSSIBLE_HAPLOTYPES, fileDest)) {
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                    haplotypes,
                    genomeLoc,
                    haplotypes,
                    null, // called haplotypes
                    readLikelihoods);
        }

        Assert.assertEquals(getReadCounts(outPath), 5);

        final File expectedMD5File = new File(outPath.toFile().getAbsolutePath() + ".md5");
        Assert.assertEquals(expectedMD5File.exists(), createMD5);
        if (createIndex) {
            Assert.assertNotNull(SamFiles.findIndex(outPath));
        } else {
            Assert.assertNull(SamFiles.findIndex(outPath));
        }
        return outPath;
    }

    @Test(dataProvider = "ReadsLikelikhoodData")
    public void testCreateSAMFromHeader(
            @SuppressWarnings("unused")  final String haplotypeBaseSignature,
            final List<Haplotype> haplotypes,
            final Locatable genomeLoc,
            final ReadLikelihoods <Haplotype> readLikelihoods
    ) throws IOException {
        final Path outputPath = GATKBaseTest.createTempFile("fromHeaderSAM", ".sam").toPath();

        try (HaplotypeBAMWriter haplotypeBAMWriter = HaplotypeBAMWriter.create(HaplotypeBAMWriter.WriterType.ALL_POSSIBLE_HAPLOTYPES, outputPath, false, false, samHeader)) {
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                    haplotypes,
                    genomeLoc,
                    haplotypes,
                    null, // called haplotypes
                    readLikelihoods);
        }

        final File expectedFile = new File(expectedFilePath, "fromHeaderSAM.sam");
        IntegrationTestSpec.assertEqualTextFiles(outputPath.toFile(), expectedFile);
    }

    @Test(dataProvider = "ReadsLikelikhoodData")
    public void testAllHaplotypeWriter
        (
            final String haplotypeBaseSignature,
            final List<Haplotype> haplotypes,
            final Locatable genomeLoc,
            final ReadLikelihoods <Haplotype> readLikelihoods
        )
    {
        final MockValidatingDestination mockDest = new MockValidatingDestination(haplotypeBaseSignature);

        try (final HaplotypeBAMWriter haplotypeBAMWriter = HaplotypeBAMWriter.create(HaplotypeBAMWriter.WriterType.ALL_POSSIBLE_HAPLOTYPES, mockDest)) {
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                    haplotypes,
                    genomeLoc,
                    haplotypes,
                    null, // called haplotypes
                    readLikelihoods);
        }

        Assert.assertTrue(mockDest.foundBases);
        Assert.assertTrue(mockDest.readCount == 5); // 4 samples + 1 haplotype
    }

    @Test(dataProvider = "ReadsLikelikhoodData")
    public void testCalledHaplotypeWriter
        (
            final String haplotypeBaseSignature,
            final List<Haplotype> haplotypes,
            final Locatable genomeLoc,
            final ReadLikelihoods <Haplotype> readLikelihoods
        )
    {
        final MockValidatingDestination mockDest = new MockValidatingDestination(haplotypeBaseSignature);

        Set<Haplotype> calledHaplotypes = new LinkedHashSet<>(1);
        calledHaplotypes.addAll(haplotypes);

        try (final HaplotypeBAMWriter haplotypeBAMWriter = HaplotypeBAMWriter.create(HaplotypeBAMWriter.WriterType.CALLED_HAPLOTYPES, mockDest)) {
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                    haplotypes,
                    genomeLoc,
                    haplotypes,
                    calledHaplotypes,
                    readLikelihoods);
        }

        Assert.assertTrue(mockDest.foundBases);
        Assert.assertTrue(mockDest.readCount==5); // 4 samples + 1 haplotype
    }

    @Test(dataProvider = "ReadsLikelikhoodData")
    public void testNoCalledHaplotypes
        (
            final String haplotypeBaseSignature,
            final List<Haplotype> haplotypes,
            final Locatable genomeLoc,
            final ReadLikelihoods <Haplotype> readLikelihoods
        )
    {
        final MockValidatingDestination mockDest = new MockValidatingDestination(haplotypeBaseSignature);

        try (final HaplotypeBAMWriter haplotypeBAMWriter = HaplotypeBAMWriter.create(HaplotypeBAMWriter.WriterType.CALLED_HAPLOTYPES, mockDest)) {
            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                    haplotypes,
                    genomeLoc,
                    haplotypes,
                    new LinkedHashSet<>(),
                    readLikelihoods);
        }

        Assert.assertTrue(mockDest.readCount == 0); // no called haplotypes, no reads
    }

    @Test(dataProvider = "ReadsLikelikhoodData")
    public void testDontWriteHaplotypes
        (
            final String haplotypeBaseSignature,
            final List<Haplotype> haplotypes,
            final Locatable genomeLoc,
            final ReadLikelihoods <Haplotype> readLikelihoods
        )
    {
        final MockValidatingDestination mockDest = new MockValidatingDestination(haplotypeBaseSignature);

        try (final HaplotypeBAMWriter haplotypeBAMWriter = HaplotypeBAMWriter.create(HaplotypeBAMWriter.WriterType.ALL_POSSIBLE_HAPLOTYPES, mockDest)) {
            haplotypeBAMWriter.setWriteHaplotypes(false);

            haplotypeBAMWriter.writeReadsAlignedToHaplotypes(
                    haplotypes,
                    genomeLoc,
                    haplotypes,
                    null, // called haplotypes
                    readLikelihoods);
        }

        Assert.assertFalse(mockDest.foundBases);
        Assert.assertTrue(mockDest.readCount == 4); // 4 samples + 0 haplotypes
    }

    private Haplotype makeHaplotype(final String bases) {
        final String cigar = bases.length() + "M";
        final Haplotype hap = new Haplotype(bases.getBytes());
        hap.setCigar(TextCigarCodec.decode(cigar));
        return hap;
    }

    private class MockValidatingDestination extends HaplotypeBAMDestination {
        private final String expectedBaseSignature;  // bases expected for the synthesized haplotype read

        public int readCount = 0;           // number of reads written to this destination
        public boolean foundBases = false;  // true we've seen a read that contains the expectedBaseSignature

        private MockValidatingDestination(String baseSignature) {
            super(samHeader, "testGroupID");
            expectedBaseSignature = baseSignature;
        }

        @Override
        void close() {}

        @Override
        public void add(GATKRead read) {
            readCount++;

            if (read.getBasesString().equals(expectedBaseSignature)) { // found haplotype base signature
                foundBases = true;
            }
        }
    }

    private ReadLikelihoods<Haplotype> generateReadLikelihoods(final int[] readCount) {
        final AlleleList<Haplotype> haplotypeList = generateHaplotypeList();
        final SampleList sampleList = generateSampleList(readCount.length);
        final Map<String,List<GATKRead>> readSamples = new LinkedHashMap<>(readCount.length);

        for (int i = 0; i < readCount.length; i++) {
            readSamples.put(sampleList.getSample(i), generateReadsList(i, readCount[i]));
        }

        return new ReadLikelihoods<>(sampleList, haplotypeList, readSamples);
    }

    private AlleleList<Haplotype> generateHaplotypeList() {
        final Haplotype[] haps = {
                makeHaplotype("AAAA"),
                makeHaplotype("AAAG")
        };
        return new IndexedAlleleList<>(haps);
    }

    private SampleList generateSampleList(final int sampleCount) {
        if (sampleCount < 0) {
            throw new IllegalArgumentException("the number of sample cannot be negative");
        }
        final List<String> result = new ArrayList<>(sampleCount);
        for (int i = 0; i < sampleCount; i++)
            result.add("SAMPLE_" + i);
        return new IndexedSampleList(result);
    }

    private List<GATKRead> generateReadsList(final int sampleIndex, final int readCount) {
        final List<GATKRead> reads = new ArrayList<>(readCount);
        int readIndex = 0;
        for (int j = 0; j < readCount; j++)
            reads.add(ArtificialReadUtils.createArtificialRead(samHeader, "READ_" + sampleIndex + "_" + (readIndex++), 1, 1, 100));
        return reads;
    }

    private int getReadCounts(final Path result) throws IOException {
        IOUtil.assertFileIsReadable(result);

        int count = 0;
        try (final SamReader in = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(result)) {
            for (@SuppressWarnings("unused") final SAMRecord rec : in) {
                count++;
            }
        }
        return count;
    }

}
