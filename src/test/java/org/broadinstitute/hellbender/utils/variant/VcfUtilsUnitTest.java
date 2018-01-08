package org.broadinstitute.hellbender.utils.variant;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;
import java.nio.file.Path;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class VcfUtilsUnitTest extends GATKBaseTest {

    VCFHeaderVersion vcfHeaderVersion = VCFHeaderVersion.VCF4_2;

    @DataProvider(name = "testData")
    public Object[][] getInitialData() {
    return new Object[][] {
            { createHeaderLines(), createSequenceDictonary(), IOUtils.getPath(b37_reference_20_21), true,  "human_g1k_v37.20.21"},
            { createHeaderLines(), createSequenceDictonary(), null, true,  "human_g1k_v37.20.21"},
            { createHeaderLines(), createSequenceDictonary(), IOUtils.getPath(b37_reference_20_21), false,
                    IOUtils.getPath(b37_reference_20_21).toUri().toString() }
        };
    }

    @Test(dataProvider = "testData")
    public void testUpdateContigsReferenceNameOnly(
            final Set<VCFHeaderLine> inHeaderLines,
            final SAMSequenceDictionary seqDict,
            final Path referencePath,
            final boolean refNameOnly,
            final String expectedRefName) {

        Set<VCFHeaderLine> resultLines = VcfUtils.updateHeaderContigLines(
                inHeaderLines, referencePath, seqDict, refNameOnly
        );

        Assert.assertEquals(resultLines.size(), referencePath == null ? 2 : 3);
        for (VCFHeaderLine resultLine : resultLines) {
            if (resultLine instanceof VCFContigHeaderLine) {
                VCFContigHeaderLine headerLine = (VCFContigHeaderLine) resultLine;
                SAMSequenceRecord samSeqRec = headerLine.getSAMSequenceRecord();
                if (samSeqRec.getSequenceName().equals("contig1")) {
                    Assert.assertEquals(samSeqRec.getSequenceLength(), 100);
                } else if (samSeqRec.getSequenceName().equals("contig2")) {
                    Assert.assertEquals(samSeqRec.getSequenceLength(), 200);
                } else {
                    Assert.fail("Bad sequence name in header lines");
                }
            } else {
                if (referencePath != null) {
                    Assert.assertEquals(resultLine.getValue(), expectedRefName);
                }
                else {
                    Assert.fail("Unexpected reference name in header lines");
                }
            }
        }
    }

    private Set<VCFHeaderLine> createHeaderLines() {
        Set<VCFHeaderLine> headerLines = new HashSet<>(2);
        headerLines.add(new VCFContigHeaderLine(
                "##contig=<ID=1,length=249250621,assembly=b37>",
                vcfHeaderVersion,
                "",
                0));
        headerLines.add(new VCFContigHeaderLine(
                "##contig=<ID=2,length=249250622,assembly=b37>",
                vcfHeaderVersion,
                "",
                0));
        return headerLines;
    }

    private SAMSequenceDictionary createSequenceDictonary() {
        List<SAMSequenceRecord> seqRecList = new ArrayList<>(2);
        seqRecList.add(new SAMSequenceRecord("contig1", 100));
        seqRecList.add(new SAMSequenceRecord("contig2", 200));
        return new SAMSequenceDictionary(seqRecList);
    }

}