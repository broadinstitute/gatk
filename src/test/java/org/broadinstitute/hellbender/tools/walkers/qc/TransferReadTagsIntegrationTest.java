package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2TestingUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class TransferReadTagsIntegrationTest extends CommandLineProgramTest {
    @DataProvider(name = "testSAMPairs")
    public Object[][] getTestSAMPairs(){
        final File alignedSAMFile = createTempFile("sam1", ".bam");
        final File alignedSAMWithExistingTagFile = createTempFile("sam2", ".bam");
        final File unmappedSAMFile = createTempFile("sam3", ".bam");

        final SAMFileHeader samHeader = M2TestingUtils.createSamHeader("sample");
        samHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        // Often an unaligned bam does not have the SortOrder field populated.
        final SAMFileHeader unalignedHeader = M2TestingUtils.createSamHeader("sample");

        final List<String> umis = Arrays.asList("TCC-ATG", "CTA-GGA", "CCT-AGC", "GCA-AAA");
        final int readLength = 30;
        final byte[] bases = new byte[readLength];
        final byte[] quals = new byte[readLength];
        Arrays.fill(bases, (byte)'A');
        Arrays.fill(quals, (byte)30);

        try (final SAMFileGATKReadWriter alignedSAMWriter = new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(alignedSAMFile, null, samHeader, true, false, false));
             final SAMFileGATKReadWriter alignedSAMWithExisingTagWriter = new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(alignedSAMWithExistingTagFile, null, samHeader, true, false, false));
             final SAMFileGATKReadWriter unmappedWriter = new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(unmappedSAMFile, null, unalignedHeader, true, false, false))) {

            for (int i = 0; i < umis.size(); i++) {
                final GATKRead unmappedRead = ArtificialReadUtils.createArtificialUnmappedRead(samHeader, bases, quals);
                unmappedRead.setName("read" + i);
                unmappedRead.setAttribute("RX", umis.get(i));
                unmappedWriter.addRead(unmappedRead);
            }

            for (int i = 0; i < umis.size(); i++) {
                final GATKRead alignedRead = ArtificialReadUtils.createArtificialRead(samHeader, bases, quals, "30M");
                alignedRead.setName("read" + i);

                // Simulate the case where an aligner drops a read
                if (i == 1) {
                    continue;
                }

                alignedSAMWriter.addRead(alignedRead);

                // When there the same tag exists in the aligned bam, it should be overwritten.
                final GATKRead alignedReadWithExistingTag = alignedRead.copy();
                alignedReadWithExistingTag.setAttribute("RX", "ZZZ-ZZZ");
                alignedSAMWithExisingTagWriter.addRead(alignedReadWithExistingTag);
            }

        } catch (RuntimeIOException e){
            throw new UserException("Failed to open SAM writers", e);
        }
        return new Object[][]{{alignedSAMFile, unmappedSAMFile, umis },
                {alignedSAMWithExistingTagFile, unmappedSAMFile, umis }
        };
    }

    @Test(dataProvider = "testSAMPairs")
    public void test(final File alignedSAMFile, final File unmappedSAMFile, List<String> umis ) {
        final File outputFile = createTempFile("output", ".bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("I", alignedSAMFile.getAbsolutePath())
                .add("unmapped-sam", unmappedSAMFile.getAbsolutePath())
                .add("read-tags", "RX")
                .add("O", outputFile.getAbsolutePath());
        runCommandLine(args, TransferReadTags.class.getSimpleName());

        final ReadsPathDataSource outputReadSource = new ReadsPathDataSource(outputFile.toPath());
        final Iterator<GATKRead> outputSamIterator = outputReadSource.iterator();

        for (int i = 0; i < umis.size(); i++){
            // Skip the same read we removed above
            if (i == 1){
                continue;
            }

            final GATKRead outputRead = outputSamIterator.next();
            Assert.assertEquals(outputRead.getName(), "read" + i);
            Assert.assertEquals(outputRead.getAttributeAsString("RX"), umis.get(i));
        }
    }
}