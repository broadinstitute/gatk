package org.broadinstitute.hellbender.tools.walkers.qc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2TestingUtils;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMFileGATKReadWriter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class RecoverReadTagsIntegrationTest extends CommandLineProgramTest {

    @Test
    public void test2() throws IOException {
        final File alignedSAMFile = File.createTempFile("sam1", ".bam");
        final File unmappedSAMFile = File.createTempFile("sam2", ".bam");

        final SAMFileHeader samHeader = M2TestingUtils.createSamHeader("sample");
        final SAMFileGATKReadWriter writer = new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(alignedSAMFile, null, samHeader, true, true, false));
        final SAMFileGATKReadWriter unmappedWriter = new SAMFileGATKReadWriter(ReadUtils.createCommonSAMWriter(unmappedSAMFile, null, samHeader, true, true, false));


        final List<String> umis = Arrays.asList("TCC-ATG", "CTA-GGA", "CCT-AGC", "GCA-AAA");
        final int readLength = 30;
        final byte[] bases = new byte[readLength];
        final byte[] quals = new byte[readLength];
        Arrays.fill(bases, (byte)'A');
        Arrays.fill(quals, (byte)30);

        for (int i = 0; i < umis.size(); i++){
            final GATKRead unmappedRead = ArtificialReadUtils.createArtificialUnmappedRead(samHeader, bases, quals);
            unmappedRead.setName("read" + i);
            unmappedRead.setAttribute("RX", umis.get(i));
            unmappedWriter.addRead(unmappedRead);

//            read.setReadGroup(DEFAULT_READ_GROUP_NAME);
//            read.setMappingQuality(DEFAULT_MAPQ);
//            read.setIsFirstOfPair();
//            read.setIsReverseStrand(i % 2 == 0);
        }

        for (int i = 0; i < umis.size(); i++){
            final GATKRead alignedRead = ArtificialReadUtils.createArtificialRead(samHeader, bases, quals, "30M");
            alignedRead.setName("read" + i);

            // Simulated the case where an aligner drops a read
            if (i == 1){
                continue;
            }

            writer.addRead(alignedRead);

//            read.setReadGroup(DEFAULT_READ_GROUP_NAME);
//            read.setMappingQuality(DEFAULT_MAPQ);
//            read.setIsFirstOfPair();
//            read.setIsReverseStrand(i % 2 == 0);
        }

        // sato: just do tr with resources...
        writer.close();
        unmappedWriter.close();

        final File outputFile = File.createTempFile("output", ".bam");
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("I", alignedSAMFile.getAbsolutePath())
                .add("unmapped-sam", unmappedSAMFile.getAbsolutePath())
                .add("O", outputFile.getAbsolutePath());
        runCommandLine(args, RecoverReadTags.class.getSimpleName());

        int d = 3;
        final ReadsPathDataSource outputReadSource = new ReadsPathDataSource(outputFile.toPath());
        final Iterator<GATKRead> outputSamIterator = outputReadSource.iterator();

        for (int i = 0; i < umis.size(); i++){
            // Skip the same read we removed above
            if (i == 1){
                continue;
            }
            final GATKRead outputRead = outputSamIterator.next();
            Assert.assertEquals(outputRead.getName(), "read" + i);
            outputRead.setName("read" + i);
            Assert.assertEquals(outputRead.getAttributeAsString("RX"), umis.get(i));
        }
    }

    @Test
    public void test(){
        final String home = "/Volumes/dsde_working/tsato/hydro.gen/Analysis/874_twist_RNA/add_umi/";

        final File alignedBam = new File(home + "SM-KYN26_SSIV_SM-LQZZ2_transcriptome.grouped.queryname_sorted.bam");
        final File unmappedSam = new File(home + "SM-KYN26_SSIV_SM-LQZZ2_UMI_extracted.bam");
        final File out = new File(home + "SM-KYN26_SSIV_SM-LQZZ2_transcriptome_with_UMI.bam");

//        final File alignedBam = new File(home + "SM-LVFDV_15_Min_Low_High_transcriptome.grouped.queryname_sorted.bam");
//        final File unmappedSam = new File(home + "SM-LVFDV_15_Min_Low_High_UMI_extracted.bam");
//        final File out = new File(home + "SM-LVFDV_15_Min_Low_High_transcriptome_with_UMI.bam");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .add("I", alignedBam)
                .add("unmapped-sam", unmappedSam)
                .add("O", out.getAbsolutePath());
        runCommandLine(args, RecoverReadTags.class.getSimpleName());
        int d = 3;
    }


}