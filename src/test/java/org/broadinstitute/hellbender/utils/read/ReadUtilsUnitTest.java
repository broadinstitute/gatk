package org.broadinstitute.hellbender.utils.read;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.RandomDNA;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;


public final class ReadUtilsUnitTest extends GATKBaseTest {
    private interface GetAdaptorFunc {
        public int getAdaptor(final GATKRead record);
    }

    @DataProvider(name = "AdaptorGetter")
    public Object[][] makeActiveRegionCutTests() {
        final List<Object[]> tests = new LinkedList<>();

        tests.add( new Object[]{ new GetAdaptorFunc() {
            @Override public int getAdaptor(final GATKRead record) { return ReadUtils.getAdaptorBoundary(record); }
        }});

        tests.add( new Object[]{ new GetAdaptorFunc() {
            @Override public int getAdaptor(final GATKRead record) { return ReadUtils.getAdaptorBoundary(record); }
        }});

        return tests.toArray(new Object[][]{});
    }

    private GATKRead makeRead(final int fragmentSize, final int mateStart) {
        final byte[] bases = {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'};
        final byte[] quals = {30, 30, 30, 30, 30, 30, 30, 30};
        final String cigar = "8M";
        GATKRead read = ArtificialReadUtils.createArtificialRead(bases, quals, cigar);
        read.setIsProperlyPaired(true);
        read.setIsPaired(true);
        read.setMatePosition(read.getContig(), mateStart);
        read.setFragmentLength(fragmentSize);
        return read;
    }

    @Test(dataProvider = "AdaptorGetter")
    public void testGetAdaptorBoundary(final GetAdaptorFunc get) {
        final int fragmentSize = 10;
        final int mateStart = 1000;
        final int BEFORE = mateStart - 2;
        final int AFTER = mateStart + 2;
        int myStart, boundary;
        GATKRead read;

        // Test case 1: positive strand, first read
        read = makeRead(fragmentSize, mateStart);
        myStart = BEFORE;
        read.setPosition(read.getContig(), myStart);
        read.setIsReverseStrand(false);
        read.setMateIsReverseStrand(true);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, myStart + fragmentSize);

        // Test case 2: positive strand, second read
        read = makeRead(fragmentSize, mateStart);
        myStart = AFTER;
        read.setPosition(read.getContig(), myStart);
        read.setIsReverseStrand(false);
        read.setMateIsReverseStrand(true);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, myStart + fragmentSize);

        // Test case 3: negative strand, second read
        read = makeRead(fragmentSize, mateStart);
        myStart = AFTER;
        read.setPosition(read.getContig(), myStart);
        read.setIsReverseStrand(true);
        read.setMateIsReverseStrand(false);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, mateStart - 1);

        // Test case 4: negative strand, first read
        read = makeRead(fragmentSize, mateStart);
        myStart = BEFORE;
        read.setPosition(read.getContig(), myStart);
        read.setIsReverseStrand(true);
        read.setMateIsReverseStrand(false);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, mateStart - 1);

        // Test case 5: mate is mapped to another chromosome (test both strands)
        read = makeRead(fragmentSize, mateStart);
        read.setFragmentLength(0);
        read.setIsReverseStrand(true);
        read.setMateIsReverseStrand(false);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);
        read.setIsReverseStrand(false);
        read.setMateIsReverseStrand(true);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);
        read.setFragmentLength(10);

        // Test case 6: read is unmapped
        read = makeRead(fragmentSize, mateStart);
        read.setIsUnmapped();
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);

        // Test case 7:  reads don't overlap and look like this:
        //    <--------|
        //                 |------>
        // first read:
        read = makeRead(fragmentSize, mateStart);
        myStart = 980;
        read.setPosition(read.getContig(), myStart);
        read.setFragmentLength(20);
        read.setIsReverseStrand(true);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);

        // second read:
        read = makeRead(fragmentSize, mateStart);
        myStart = 1000;
        read.setPosition(read.getContig(), myStart);
        read.setFragmentLength(20);
        read.setMatePosition(read.getContig(), 980);
        read.setIsReverseStrand(false);
        boundary = get.getAdaptor(read);
        Assert.assertEquals(boundary, ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);

        // Test case 8: read doesn't have proper pair flag set
        read = makeRead(fragmentSize, mateStart);
        read.setIsPaired(true);
        read.setIsProperlyPaired(false);
        Assert.assertEquals(get.getAdaptor(read), ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY);

        // Test case 9: read and mate have same negative flag setting
        for ( final boolean negFlag: Arrays.asList(true, false) ) {
            read = makeRead(fragmentSize, mateStart);
            read.setPosition(read.getContig(), BEFORE);
            read.setIsPaired(true);
            read.setIsProperlyPaired(true);
            read.setIsReverseStrand(negFlag);
            read.setMateIsReverseStrand(!negFlag);
            Assert.assertTrue(get.getAdaptor(read) != ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY, "Get adaptor should have succeeded");

            read = makeRead(fragmentSize, mateStart);
            read.setPosition(read.getContig(), BEFORE);
            read.setIsPaired(true);
            read.setIsProperlyPaired(true);
            read.setIsReverseStrand(negFlag);
            read.setMateIsReverseStrand(negFlag);
            Assert.assertEquals(get.getAdaptor(read), ReadUtils.CANNOT_COMPUTE_ADAPTOR_BOUNDARY, "Get adaptor should have failed for reads with bad alignment orientation");
        }
    }

    @Test
    public void testGetBasesReverseComplement() {
        int iterations = 1000;
        Random random = Utils.getRandomGenerator();
        while(iterations-- > 0) {
            final int l = random.nextInt(1000);
            GATKRead read = ArtificialReadUtils.createRandomRead(l);
            byte [] original = read.getBases();
            byte [] reconverted = new byte[l];
            String revComp = ReadUtils.getBasesReverseComplement(read);
            for (int i=0; i<l; i++) {
                reconverted[l-1-i] = BaseUtils.getComplement((byte) revComp.charAt(i));
            }
            Assert.assertEquals(reconverted, original);
        }
    }

    @Test
    public void testGetMaxReadLength() {
        for( final int minLength : Arrays.asList( 5, 30, 50 ) ) {
            for( final int maxLength : Arrays.asList( 50, 75, 100 ) ) {
                final List<GATKRead> reads = new ArrayList<>();
                for( int readLength = minLength; readLength <= maxLength; readLength++ ) {
                    reads.add( ArtificialReadUtils.createRandomRead( readLength ) );
                }
                Assert.assertEquals(ReadUtils.getMaxReadLength(reads), maxLength, "max length does not match");
            }
        }

        final List<GATKRead> reads = new LinkedList<>();
        Assert.assertEquals(ReadUtils.getMaxReadLength(reads), 0, "Empty list should have max length of zero");
    }

    @Test
    public void testReadWithNsRefIndexInDeletion() throws FileNotFoundException {

        final IndexedFastaSequenceFile seq = new CachingIndexedFastaSequenceFile(IOUtils.getPath(exampleReference));
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final int readLength = 76;

        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 8975, readLength);
        read.setBases(Utils.dupBytes((byte) 'A', readLength));
        read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
        read.setCigar("3M414N1D73M");

        final int result = ReadUtils.getReadCoordinateForReferenceCoordinateUpToEndOfRead(read, 9392, ReadUtils.ClippingTail.LEFT_TAIL);
        Assert.assertEquals(result, 2);
    }

    @Test
    public void testReadWithNsRefAfterDeletion() throws FileNotFoundException {

        final IndexedFastaSequenceFile seq = new CachingIndexedFastaSequenceFile(IOUtils.getPath(exampleReference));
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final int readLength = 76;

        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "myRead", 0, 8975, readLength);
        read.setBases(Utils.dupBytes((byte) 'A', readLength));
        read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
        read.setCigar("3M414N1D73M");

        final int result = ReadUtils.getReadCoordinateForReferenceCoordinateUpToEndOfRead(read, 9393, ReadUtils.ClippingTail.LEFT_TAIL);
        Assert.assertEquals(result, 3);
    }

    @DataProvider(name = "HasWellDefinedFragmentSizeData")
    public Object[][] makeHasWellDefinedFragmentSizeData() throws Exception {
        final List<Object[]> tests = new LinkedList<>();

        // setup a basic read that will work
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "read1", 0, 10, 10);
        read.setIsPaired(true);
        read.setIsProperlyPaired(true);
        read.setPosition(read.getContig(), 100);
        read.setCigar("50M");
        read.setMatePosition(read.getContig(), 130);
        read.setFragmentLength(80);
        read.setIsFirstOfPair();
        read.setIsReverseStrand(false);
        read.setMateIsReverseStrand(true);

        tests.add( new Object[]{ "basic case", read.copy(), true });

        {
            final GATKRead bad1 = read.copy();
            bad1.setIsPaired(false);
            tests.add( new Object[]{ "not paired", bad1, false });
        }

        {
            final GATKRead bad = read.copy();
            bad.setIsProperlyPaired(false);
            // we currently don't require the proper pair flag to be set
            tests.add( new Object[]{ "not proper pair", bad, true });
//            tests.add( new Object[]{ "not proper pair", bad, false });
        }

        {
            final GATKRead bad = read.copy();
            bad.setIsUnmapped();
            tests.add( new Object[]{ "read is unmapped", bad, false });
        }

        {
            final GATKRead bad = read.copy();
            bad.setMateIsUnmapped();
            tests.add( new Object[]{ "mate is unmapped", bad, false });
        }

        {
            final GATKRead bad = read.copy();
            bad.setMateIsReverseStrand(false);
            tests.add( new Object[]{ "read and mate both on positive strand", bad, false });
        }

        {
            final GATKRead bad = read.copy();
            bad.setIsReverseStrand(true);
            tests.add( new Object[]{ "read and mate both on negative strand", bad, false });
        }

        {
            final GATKRead bad = read.copy();
            bad.setFragmentLength(0);
            tests.add( new Object[]{ "insert size is 0", bad, false });
        }

        {
            final GATKRead bad = read.copy();
            bad.setPosition(bad.getContig(), 1000);
            tests.add( new Object[]{ "positve read starts after mate end", bad, false });
        }

        {
            final GATKRead bad = read.copy();
            bad.setIsReverseStrand(true);
            bad.setMateIsReverseStrand(false);
            bad.setMatePosition(bad.getMateContig(), 1000);
            tests.add( new Object[]{ "negative strand read ends before mate starts", bad, false });
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "HasWellDefinedFragmentSizeData")
    private void testHasWellDefinedFragmentSize(final String name, final GATKRead read, final boolean expected) {
        Assert.assertEquals(ReadUtils.hasWellDefinedFragmentSize(read), expected);
    }

    @DataProvider(name = "ReadsWithReadGroupData")
    public Object[][] readsWithReadGroupData() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(2, 1, 1000000);
        final SAMReadGroupRecord readGroup = new SAMReadGroupRecord("FOO");
        readGroup.setPlatform("FOOPLATFORM");
        readGroup.setPlatformUnit("FOOPLATFORMUNIT");
        readGroup.setLibrary("FOOLIBRARY");
        readGroup.setSample("FOOSAMPLE");
        header.addReadGroup(readGroup);

        final GATKRead samBackedRead = new SAMRecordToGATKReadAdapter(ArtificialReadUtils.createArtificialSAMRecord(header, "sam", header.getSequenceIndex("1"), 5, new byte[]{'A', 'C', 'G', 'T'}, new byte[]{1, 2, 3, 4}, "4M"));
        samBackedRead.setReadGroup("FOO");

        return new Object[][] {
                { samBackedRead, header, "FOO" }
        };
    }

    @Test(dataProvider = "ReadsWithReadGroupData")
    public void testReadGroupOperations( final GATKRead read, final SAMFileHeader header, final String expectedReadGroupID ) {
        final SAMReadGroupRecord readGroup = ReadUtils.getSAMReadGroupRecord(read, header);
        Assert.assertEquals(readGroup.getId(), expectedReadGroupID, "Wrong read group returned from ReadUtils.getSAMReadGroupRecord()");

        Assert.assertEquals(ReadUtils.getPlatform(read, header), readGroup.getPlatform(), "Wrong platform returned from ReadUtils.getPlatform()");
        Assert.assertEquals(ReadUtils.getPlatformUnit(read, header), readGroup.getPlatformUnit(), "Wrong platform unit returned from ReadUtils.getPlatformUnit()");
        Assert.assertEquals(ReadUtils.getLibrary(read, header), readGroup.getLibrary(), "Wrong library returned from ReadUtils.getLibrary()");
        Assert.assertEquals(ReadUtils.getSampleName(read, header), readGroup.getSample(), "Wrong sample name returned from ReadUtils.getSampleName()");
    }

    @Test(dataProvider = "ReadsWithReadGroupData")
    public void testReadGroupOperationsOnReadWithNoReadGroup( final GATKRead read, final SAMFileHeader header, final String expectedReadGroupID ) {
        read.setReadGroup(null);

        Assert.assertNull(ReadUtils.getSAMReadGroupRecord(read, header), "Read group should be null");
        Assert.assertNull(ReadUtils.getPlatform(read, header), "Platform should be null");
        Assert.assertNull(ReadUtils.getPlatformUnit(read, header), "Platform unit should be null");
        Assert.assertNull(ReadUtils.getLibrary(read, header), "Library should be null");
        Assert.assertNull(ReadUtils.getSampleName(read, header), "Sample name should be null");
    }

    @DataProvider(name = "ReferenceIndexTestData")
    public Object[][] referenceIndexTestData() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(2, 1, 1000000);

        final GATKRead samBackedRead = new SAMRecordToGATKReadAdapter(ArtificialReadUtils.createArtificialSAMRecord(header, "sam", header.getSequenceIndex("2"), 5, new byte[]{'A', 'C', 'G', 'T'}, new byte[]{1, 2, 3, 4}, "4M"));
        samBackedRead.setMatePosition("1", 5);

        final GATKRead unmappedSamBackedRead = new SAMRecordToGATKReadAdapter(ArtificialReadUtils.createArtificialSAMRecord(header, "sam", header.getSequenceIndex("1"), 5, new byte[]{'A', 'C', 'G', 'T'}, new byte[]{1, 2, 3, 4}, "4M"));
        unmappedSamBackedRead.setIsUnmapped();
        unmappedSamBackedRead.setMateIsUnmapped();

        return new Object[][] {
                { samBackedRead, header, 1, 0 },
                { unmappedSamBackedRead, header, -1, -1 }
        };
    }

    @Test(dataProvider = "ReferenceIndexTestData")
    public void testGetReferenceIndex( final GATKRead read, final SAMFileHeader header, final int expectedReferenceIndex, final int expectedMateReferenceIndex ) {
        Assert.assertEquals(ReadUtils.getReferenceIndex(read, header), expectedReferenceIndex, "Wrong reference index for read");
        Assert.assertEquals(ReadUtils.getMateReferenceIndex(read, header), expectedMateReferenceIndex, "Wrong reference index for read's mate");
    }

    @Test
    public void testGetAssignedReferenceIndex() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final GATKRead mappedRead = ArtificialReadUtils.createArtificialRead(header, "foo", 0, 5, 10);
        final GATKRead unmappedRead = ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{'A'}, new byte[]{30});
        final GATKRead unmappedReadWithAssignedPosition = ArtificialReadUtils.createArtificialUnmappedReadWithAssignedPosition(header, "2", 10, new byte[]{'A'}, new byte[]{30});

        Assert.assertEquals(ReadUtils.getAssignedReferenceIndex(mappedRead, header), 0);
        Assert.assertEquals(ReadUtils.getAssignedReferenceIndex(unmappedRead, header), -1);
        Assert.assertEquals(ReadUtils.getAssignedReferenceIndex(unmappedReadWithAssignedPosition, header), 1);
    }

    @DataProvider(name = "ReadHasNoAssignedPositionTestData")
    public Object[][] readHasNoAssignedPositionTestData() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();

        // We'll also test the improbable case of a SAMRecord marked as mapped, but with an unset contig/start
        final SAMRecord mappedSAMWithUnsetContigSetStart = new SAMRecord(header);
        mappedSAMWithUnsetContigSetStart.setReferenceName(ReadConstants.UNSET_CONTIG);
        mappedSAMWithUnsetContigSetStart.setAlignmentStart(10);
        mappedSAMWithUnsetContigSetStart.setReadUnmappedFlag(false);
        final GATKRead mappedReadWithUnsetContigSetStart = new SAMRecordToGATKReadAdapter(mappedSAMWithUnsetContigSetStart);

        final SAMRecord mappedSAMWithSetContigUnsetStart = new SAMRecord(header);
        mappedSAMWithSetContigUnsetStart.setReferenceName("1");
        mappedSAMWithSetContigUnsetStart.setAlignmentStart(ReadConstants.UNSET_POSITION);
        mappedSAMWithSetContigUnsetStart.setReadUnmappedFlag(false);
        final GATKRead mappedReadWithSetContigUnsetStart = new SAMRecordToGATKReadAdapter(mappedSAMWithSetContigUnsetStart);

        return new Object[][] {
                // Mapped read with position
                { ArtificialReadUtils.createArtificialRead(header, "foo", 0, 5, 10), false },
                // Basic unmapped read with no position
                { ArtificialReadUtils.createArtificialUnmappedRead(header, new byte[]{'A'}, new byte[]{30}), true },
                // Unmapped read with set position (contig and start)
                { ArtificialReadUtils.createArtificialUnmappedReadWithAssignedPosition(header, "1", 10, new byte[]{'A'}, new byte[]{30}), false },
                // Unmapped read with "*" contig, set start
                { ArtificialReadUtils.createArtificialUnmappedReadWithAssignedPosition(header, ReadConstants.UNSET_CONTIG, 10, new byte[]{'A'}, new byte[]{30}), true },
                // Unmapped read with set contig, unset start
                { ArtificialReadUtils.createArtificialUnmappedReadWithAssignedPosition(header, "1", ReadConstants.UNSET_POSITION, new byte[]{'A'}, new byte[]{30}), true },
                // Unmapped read with "*" contig, unset start
                { ArtificialReadUtils.createArtificialUnmappedReadWithAssignedPosition(header, ReadConstants.UNSET_CONTIG, ReadConstants.UNSET_POSITION, new byte[]{'A'}, new byte[]{30}), true },
                // "Mapped" read with unset contig, set start
                { mappedReadWithUnsetContigSetStart, true },
                // "Mapped" read with set contig, unset start
                { mappedReadWithSetContigUnsetStart, true }
        };
    }

    @Test(dataProvider = "ReadHasNoAssignedPositionTestData")
    public void testReadHasNoAssignedPosition( final GATKRead read, final boolean expectedResult ) {
        Assert.assertEquals(ReadUtils.readHasNoAssignedPosition(read), expectedResult);
    }

    @DataProvider(name="createSAMWriter")
    public Object[][] createSAMWriterData() {
        return new Object[][] {
                // Note: We expect to silently fail to create an index if createIndex is true but sort order is not coord.
            {getTestFile("print_reads.sam"),     getTestFile("print_reads.fasta"), ".cram", false,  true, true, false},
            {getTestFile("query_sorted.bam"),     null, ".bam", false,  true, true, false},
            {getTestFile("coordinate_sorted.bam"),null, ".bam", false,  true, true, true},
            {getTestFile("query_sorted.bam"),     null, ".bam", true,   true, true, false},
            {getTestFile("coordinate_sorted.bam"),null, ".bam", true,   true, true, true},
            {getTestFile("query_sorted.bam"),     null, ".bam", true,   true, false, false},
            {getTestFile("coordinate_sorted.bam"),null, ".bam", true,   true, false, true},
            {getTestFile("coordinate_sorted.bam"),null, ".bam", true,   false, false, false}        };
    }

    @Test(dataProvider="createSAMWriter")
    public void testCreateLegacySAMWriter(
            final File bamFile,
            final File referenceFile,
            final String outputExtension,
            final boolean preSorted,
            final boolean createIndex,
            final boolean createMD5,
            final boolean expectIndex) throws Exception {
        final File outputFile = createTempFile("samWriterTest",  outputExtension);

        try (final SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(referenceFile).open(bamFile)) {
             final SAMFileHeader header = samReader.getFileHeader();
            if (expectIndex) { // ensure test condition
                Assert.assertEquals(expectIndex, header.getSortOrder() == SAMFileHeader.SortOrder.coordinate);
            }

            try (final SAMFileWriter samWriter = ReadUtils.createCommonSAMWriter
                            (outputFile, referenceFile, samReader.getFileHeader(), preSorted, createIndex, createMD5)) {
                final Iterator<SAMRecord> samRecIt = samReader.iterator();
                while (samRecIt.hasNext()) {
                    samWriter.addAlignment(samRecIt.next());
                }
            }
        }

        final File md5File = new File(outputFile.getAbsolutePath() + ".md5");
        if (md5File.exists()) {
            md5File.deleteOnExit();
        }
        Assert.assertEquals(expectIndex, null != SamFiles.findIndex(outputFile));
        Assert.assertEquals(createMD5, md5File.exists());

        // now check the contents are the same
        try (final SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(referenceFile).open(bamFile);
            final SamReader outputReader = SamReaderFactory.makeDefault().referenceSequence(referenceFile).open(outputFile)) {
            final Iterator<SAMRecord> samRecIt = samReader.iterator();
            final Iterator<SAMRecord> outRecIt = outputReader.iterator();
            Assert.assertEquals(samRecIt, outRecIt);
        }
    }

    @Test(dataProvider="createSAMWriter")
    public void testCreatePathSAMWriter(
        final File bamFile,
        final File referenceFile,
        final String outputExtension,
        final boolean preSorted,
        final boolean createIndex,
        final boolean createMD5,
        final boolean expectIndex) throws Exception {

        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {

            final Path outputPath = jimfs.getPath("samWriterTest" + outputExtension);
            final Path md5Path = jimfs.getPath("samWriterTest" + outputExtension + ".md5");
            final Path referencePath;
            if( referenceFile == null) {
                referencePath = null;
            } else {
                referencePath = Files.copy(referenceFile.toPath(), jimfs.getPath(referenceFile.getName()));
                final Path fastaIndexLocal = ReferenceSequenceFileFactory.getFastaIndexFileName(referenceFile.toPath());
                Files.copy(fastaIndexLocal, jimfs.getPath(fastaIndexLocal.getFileName().toString()));
            }

            try (final SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(referencePath).open(bamFile)) {

                final SAMFileHeader header = samReader.getFileHeader();
                if (expectIndex) { // ensure test condition
                    Assert.assertEquals(expectIndex, header.getSortOrder() == SAMFileHeader.SortOrder.coordinate);
                }

                try (final SAMFileWriter samWriter = ReadUtils.createCommonSAMWriter
                    (outputPath, referencePath, samReader.getFileHeader(), preSorted, createIndex, createMD5)) {
                    final Iterator<SAMRecord> samRecIt = samReader.iterator();
                    while (samRecIt.hasNext()) {
                        samWriter.addAlignment(samRecIt.next());
                    }
                }

                Assert.assertEquals(expectIndex, null != SamFiles.findIndex(outputPath));
                Assert.assertEquals(createMD5, Files.exists(md5Path));
            }

            // now check the contents are the same
            try (final SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(referenceFile).open(bamFile);
                final SamReader outputReader = SamReaderFactory.makeDefault().referenceSequence(referenceFile).open(outputPath)) {
                final Iterator<SAMRecord> samRecIt = samReader.iterator();
                final Iterator<SAMRecord> outRecIt = outputReader.iterator();
                Assert.assertEquals(samRecIt, outRecIt);
            }
        }
    }

    @DataProvider(name="hasCRAMFileContents")
    public Object[][] createHasCRAMFileContentsData() {
        return new Object[][] {
                {getTestFile("valid.sam"), false},
                {getTestFile("coordinate_sorted.bam"), false},
                {getTestFile("valid.cram"), true},
                {getTestFile("fake_cram_with_bam_contents.cram"), false}
        };
    }

    @Test(dataProvider = "hasCRAMFileContents")
    public void testHasCRAMFileContents(final File testFile, final boolean expected) {
        Assert.assertEquals(ReadUtils.hasCRAMFileContents(testFile), expected);
    }

    @Test
    public void testGetSamplesFromHeader() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(5, 1, 100);
        final List<SAMReadGroupRecord> readGroups = new ArrayList<>();
        for ( int i = 1; i <= 5; ++i ) {
            SAMReadGroupRecord readGroup = new SAMReadGroupRecord("ReadGroup" + i);
            readGroup.setSample("Sample" + i);
            readGroups.add(readGroup);
        }
        header.setReadGroups(readGroups);

        final Set<String> samples = ReadUtils.getSamplesFromHeader(header);
        Assert.assertEquals(samples.size(), 5, "Wrong number of samples returned from ReadUtils.getSamplesFromHeader()");
        for ( int i = 1; i <= 5; ++i ) {
            Assert.assertTrue(samples.contains("Sample" + i), "Missing Sample" + i + " in samples returned from ReadUtils.getSamplesFromHeader()");
        }
    }

    @Test
    public void testGetSamplesFromHeaderNoReadGroups() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(5, 1, 100);
        Assert.assertTrue(header.getReadGroups().isEmpty());

        Assert.assertTrue(ReadUtils.getSamplesFromHeader(header).isEmpty(), "Non-empty Set returned from ReadUtils.getSamplesFromHeader() for a header with no read groups");
    }

    @Test
    public void testGetSamplesFromHeaderNoSamples() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(5, 1, 100);
        header.setReadGroups(Arrays.asList(new SAMReadGroupRecord("ReadGroup1")));
        Assert.assertEquals(header.getReadGroups().size(), 1);
        Assert.assertNull(header.getReadGroups().get(0).getSample());

        Assert.assertTrue(ReadUtils.getSamplesFromHeader(header).isEmpty(), "Non-empty Set returned from ReadUtils.getSamplesFromHeader() for a header with no samples");
    }

    @DataProvider(name="testValidateSortOrder")
    public Object[][] createValidateSortOrderData() {
        return new Object[][] {
                { SAMFileHeader.SortOrder.coordinate, SAMFileHeader.SortOrder.coordinate, true, true},
                { SAMFileHeader.SortOrder.coordinate, SAMFileHeader.SortOrder.coordinate, false, true},
                { SAMFileHeader.SortOrder.queryname, SAMFileHeader.SortOrder.queryname, true, true},
                { SAMFileHeader.SortOrder.queryname, SAMFileHeader.SortOrder.queryname, false, true},
                { SAMFileHeader.SortOrder.coordinate, SAMFileHeader.SortOrder.unsorted, true, true},
                { SAMFileHeader.SortOrder.coordinate, SAMFileHeader.SortOrder.unsorted, false, true},
                { SAMFileHeader.SortOrder.queryname, SAMFileHeader.SortOrder.coordinate, true, false},
        };
    }

    @Test(dataProvider = "testValidateSortOrder")
    public void testValidateExpectedReadOrder(
            final SAMFileHeader.SortOrder actualSortOrder,
            final SAMFileHeader.SortOrder expectedSortOrder,
            final boolean assumeSorted,
            final boolean expectedValid) {
        boolean isValid = ReadUtils.validateExpectedSortOrder(
                actualSortOrder,
                expectedSortOrder,
                assumeSorted,
                "test"
        );
        Assert.assertEquals(isValid, expectedValid);
    }

    @DataProvider(name="testInvalidSortOrder")
    public Object[][] createInvalidValidateSortOrderData() {
        return new Object[][] {
                { SAMFileHeader.SortOrder.coordinate, SAMFileHeader.SortOrder.queryname, false, false},
                { SAMFileHeader.SortOrder.queryname, SAMFileHeader.SortOrder.coordinate, false, false},
        };
    }

    @Test(dataProvider = "testInvalidSortOrder", expectedExceptions = UserException.class)
    public void testNotValidExpectedReadOrder(
        final SAMFileHeader.SortOrder actualSortOrder,
        final SAMFileHeader.SortOrder expectedSortOrder,
        final boolean assumeSorted,
        final boolean expectedValid) {
            ReadUtils.validateExpectedSortOrder(
                    actualSortOrder,
                    expectedSortOrder,
                    assumeSorted,
                    "test"
            );
    }

    @Test
    public void testOptionalIntAttributeOnSAMRecord() {
        final SAMRecord record = ArtificialReadUtils.createArtificialSAMRecord(
                new Cigar(Collections.singletonList(new CigarElement(100, CigarOperator.M))));
        Assert.assertFalse(ReadUtils.getOptionalIntAttribute(record, SAMTag.AS.name()).isPresent());
        record.setAttribute(SAMTag.AS.name(), -10);
        Assert.assertEquals(ReadUtils.getOptionalIntAttribute(record, SAMTag.AS.name()).orElse(-1), -10);
        record.setAttribute(SAMTag.AS.name(), "-20");
        Assert.assertEquals(ReadUtils.getOptionalIntAttribute(record, SAMTag.AS.name()).orElse(-1), -20);

        record.setAttribute(SAMTag.AS.name(), 10.213f);
        try {
            ReadUtils.getOptionalIntAttribute(record, "AS").isPresent();
            Assert.fail("expected and exception");
        } catch (final Throwable ex) {
            Assert.assertTrue(ex instanceof GATKException.ReadAttributeTypeMismatch, "wrong ex class: " + ex.getClass());
        }

        record.setAttribute(SAMTag.AS.name(), 10.0f);
        try {
            ReadUtils.getOptionalIntAttribute(record, "AS").isPresent();
            Assert.fail("expected and exception");
        } catch (final Throwable ex) {
            Assert.assertTrue(ex instanceof GATKException.ReadAttributeTypeMismatch, "wrong ex class: " + ex.getClass());
        }

        record.setAttribute(SAMTag.AS.name(), "10.213");
        try {
            ReadUtils.getOptionalIntAttribute(record, "AS").isPresent();
            Assert.fail("expected and exception");
        } catch (final Throwable ex) {
            Assert.assertTrue(ex instanceof GATKException.ReadAttributeTypeMismatch, "wrong ex class: " + ex.getClass());
        }
        record.setAttribute(SAMTag.AS.name(), null);
        Assert.assertFalse(ReadUtils.getOptionalIntAttribute(record, "AS").isPresent());
    }

    @Test
    public void testGetOptionalIntAttributeOnGATKRead() {
        final SAMRecord record = ArtificialReadUtils.createArtificialSAMRecord(
                new Cigar(Collections.singletonList(new CigarElement(100, CigarOperator.M))));
        final GATKRead read = new SAMRecordToGATKReadAdapter(record);
        Assert.assertFalse(ReadUtils.getOptionalIntAttribute(read, SAMTag.AS.name()).isPresent());
        read.setAttribute(SAMTag.AS.name(), -10);
        Assert.assertEquals(ReadUtils.getOptionalIntAttribute(read, SAMTag.AS.name()).orElse(-1), -10);
        read.setAttribute(SAMTag.AS.name(), "-20");
        Assert.assertEquals(ReadUtils.getOptionalIntAttribute(record, SAMTag.AS.name()).orElse(-1), -20);
        read.setAttribute(SAMTag.AS.name(), "10.213");
        try {
            ReadUtils.getOptionalIntAttribute(read, SAMTag.AS.name()).isPresent();
            Assert.fail("expected and exception");
        } catch (final Throwable ex) {
            Assert.assertTrue(ex instanceof GATKException.ReadAttributeTypeMismatch, "wrong ex class: " + ex.getClass());
        }

        read.setAttribute(SAMTag.AS.name(), "10.0");
        try {
            ReadUtils.getOptionalIntAttribute(read, SAMTag.AS.name()).isPresent();
            Assert.fail("expected and exception");
        } catch (final Throwable ex) {
            Assert.assertTrue(ex instanceof GATKException.ReadAttributeTypeMismatch, "wrong ex class: " + ex.getClass());
        }

        read.clearAttribute(SAMTag.AS.name());
        Assert.assertFalse(ReadUtils.getOptionalIntAttribute(read, SAMTag.AS.name()).isPresent());
    }

    @Test(dataProvider = "mappedGatkReadsData")
    public void testHasBasesAligned(final GATKRead read) {
        Assert.assertTrue(ReadUtils.hasBasesAlignedAgainstTheReference(read));
    }

    @Test(dataProvider = "mappedGatkReadsData")
    public void testGetFirstAlignedBaseOffset(final GATKRead read) {
        final int actual = ReadUtils.getFirstAlignedBaseOffset(read);
        final int expected = CigarUtils.countLeftClippedBases(read.getCigar()) - CigarUtils.countLeftHardClippedBases(read.getCigar());
        Assert.assertEquals(actual, expected);
    }

    @Test(dataProvider = "mappedGatkReadsData")
    public void testGetAfterLastAlignedBaseOffset(final GATKRead read) {
        final int actual = ReadUtils.getAfterLastAlignedBaseOffset(read);
        final int expected = read.getLength() - (CigarUtils.countRightClippedBases(read.getCigar()) - CigarUtils.countRightHardClippedBases(read.getCigar()));
        Assert.assertEquals(actual, expected);
    }

    @Test(dataProvider = "mappedGatkReadsData")
    public void testGetFirstPositionAligned(final GATKRead read) {
        final int actual = ReadUtils.getFirstAlignedReadPosition(read);
        final int expected = (!read.isReverseStrand()
                    ? CigarUtils.countLeftClippedBases(read.getCigar()) + 1
                    : CigarUtils.countRightClippedBases(read.getCigar()) + 1);
        Assert.assertEquals(actual, expected);
    }

    @Test(dataProvider = "mappedGatkReadsData")
    public void testGetLastPositionAligned(final GATKRead read) {
        final int actual = ReadUtils.getLastAlignedReadPosition(read);
        final int expected = (!read.isReverseStrand()
                ? CigarUtils.countUnclippedReadBases(read.getCigar()) - CigarUtils.countRightClippedBases(read.getCigar())
                : CigarUtils.countUnclippedReadBases(read.getCigar()) - CigarUtils.countLeftClippedBases(read.getCigar()));
        Assert.assertEquals(actual, expected);
    }

    @Test(dataProvider = "unmappedGatkReadData")
    public void testHashBasesNotAligned(final GATKRead read) {
        Assert.assertFalse(ReadUtils.hasBasesAlignedAgainstTheReference(read));
    }

    @Test(dataProvider = "unmappedGatkReadData", expectedExceptions = IllegalArgumentException.class)
    public void testGetFirstAlignedBaseOffsetOnUnmapped(final GATKRead read) {
        ReadUtils.getFirstAlignedBaseOffset(read);
    }

    @Test(dataProvider = "unmappedGatkReadData", expectedExceptions = IllegalArgumentException.class)
    public void testGetAfterLastAlignedBaseOffsetOnUnmapped(final GATKRead read) {
        ReadUtils.getAfterLastAlignedBaseOffset(read);
    }

    @Test(dataProvider = "unmappedGatkReadData", expectedExceptions = IllegalArgumentException.class)
    public void testGetFirtAlignedReadPositionOnUnmapped(final GATKRead read) {
        ReadUtils.getFirstAlignedReadPosition(read);
    }

    @Test(dataProvider = "unmappedGatkReadData", expectedExceptions = IllegalArgumentException.class)
    public void testGetLastAlignedReadPositionOnUnmapped(final GATKRead read) {
        ReadUtils.getLastAlignedReadPosition(read);
    }

    @DataProvider(name="unmappedGatkReadData")
    public Object[][] unmappedGatkReadData() {
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final Random rdn = new Random(17);
        final RandomDNA rdnDna = new RandomDNA(11);
        return Stream.of("*", "30I", "10S1I10S", "1D")
                .map(TextCigarCodec::decode)
                .map(c -> {
                    final SAMRecord record = new SAMRecord(header);
                    record.setReadUnmappedFlag(c.isEmpty());
                    record.setCigar(c);
                    record.setReadName("test-read");
                    record.setReadBases(rdnDna.nextBases(c.getCigarElements().stream()
                            .filter(ce -> ce.getOperator().consumesReadBases())
                            .mapToInt(CigarElement::getLength).sum()));
                    if (!c.isEmpty()) {
                        record.setReferenceIndex(rdn.nextInt(header.getSequenceDictionary().getSequences().size()));
                        record.setAlignmentStart(10_000 + header.getSequence(record.getReferenceIndex()).getSequenceLength() - 20_000);
                        record.setReadNegativeStrandFlag(rdn.nextBoolean());
                    }
                    return record;
                })
                .map(SAMRecordToGATKReadAdapter::new)
                .map(r -> new Object[] { r })
                .toArray(Object[][]::new);
    }

    @DataProvider(name="mappedGatkReadsData")
    public Object[][] mappedGatkReadsData() {
        final Random rdn = new Random(17);
        final RandomDNA rdnDna = new RandomDNA(11);
        final List<Cigar> randomCigars = CigarTestUtils.randomValidCigars(rdn, 1_000, 10, 100);
        final List<Object[]> result = new ArrayList<>(randomCigars.size());
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        for (final Cigar cigar : randomCigars) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadUnmappedFlag(false);
            record.setCigar(cigar);
            record.setReadName("test-read");
            record.setReadBases(rdnDna.nextBases(cigar.getCigarElements().stream()
                    .filter(ce -> ce.getOperator().consumesReadBases())
                    .mapToInt(CigarElement::getLength).sum()));
            record.setReferenceIndex(rdn.nextInt(header.getSequenceDictionary().getSequences().size()));
            record.setAlignmentStart(10_000 + header.getSequence(record.getReferenceIndex()).getSequenceLength() - 20_000);
            record.setReadNegativeStrandFlag(rdn.nextBoolean());
            result.add(new Object[] {new SAMRecordToGATKReadAdapter(record)});
        }

        return result.toArray(new Object[result.size()][]);
    }

    @Test
    public void testGetReadToMateMap() {
        final String bam = largeFileTestDir + "mutect/dream_synthetic_bams/tumor_1.bam";
        final ReadsDataSource readsDataSource = new ReadsDataSource(IOUtils.getPath(bam));

        // to keep this example small, we choose a place that only four reads overlap
        final SimpleInterval interval = new SimpleInterval("20", 576_940, 577_000);
        final ReadsContext readsContext = new ReadsContext(readsDataSource, interval, ReadFilterLibrary.ALLOW_ALL_READS);

        final Map<GATKRead, GATKRead> smallFragmentResult = ReadUtils.getReadToMateMap(readsContext, 10);
        Assert.assertEquals(smallFragmentResult.size(), 1);
        final Map.Entry<GATKRead, GATKRead> readAndMate = smallFragmentResult.entrySet().iterator().next();
        Assert.assertEquals(readAndMate.getKey().getName(), "0685fae6-eddc-4257-81f5-c7b9eab84f2f");
        Assert.assertTrue(readAndMate.getKey().isReverseStrand());

        final Map<GATKRead, GATKRead> mediumFragmentResult = ReadUtils.getReadToMateMap(readsContext, 200);
        Assert.assertEquals(mediumFragmentResult.size(), 2);
        final List<Map.Entry<GATKRead, GATKRead>> readsAndMates = Utils.stream(mediumFragmentResult.entrySet().iterator()).collect(Collectors.toList());
        Assert.assertTrue(readsAndMates.stream().anyMatch(pair -> pair.getKey().getName().equals("0685fae6-eddc-4257-81f5-c7b9eab84f2f")));
        Assert.assertTrue(readsAndMates.stream().anyMatch(pair -> pair.getKey().getName().equals("9eaf071f-6c59-4f10-8b77-467d359ac702")));

        final Map<GATKRead, GATKRead> bigFragmentResult = ReadUtils.getReadToMateMap(readsContext, 1000);
        Assert.assertEquals(bigFragmentResult.size(), 4);
    }
}
