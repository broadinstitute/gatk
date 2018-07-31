package org.broadinstitute.hellbender.tools.spark.sv.utils;

import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

public class GATKSVVCFUtilsUnitTest extends GATKBaseTest {


    @Test(groups = "sv")
    public void testVCFConstants() {
        Assert.assertEquals(GATKSVVCFConstants.expectedHeaderLinesInVCF,
                Stream.of(SVTYPE, SVLEN, IMPRECISE, CIPOS, CIEND, BND_MATEID_STR, SYMB_ALT_ALLELE_INV, READ_PAIR_SUPPORT,
                        SPLIT_READ_SUPPORT, SYMB_ALT_ALLELE_DEL, SYMB_ALT_ALLELE_INS, SYMB_ALT_ALLELE_DUP, SYMB_ALT_ALLELE_INVDUP,
                        CONTIG_NAMES, TOTAL_MAPPINGS, MAPPING_QUALITIES, HQ_MAPPINGS, ALIGN_LENGTHS, MAX_ALIGN_LENGTH,
                        SEQ_ALT_HAPLOTYPE, INSERTED_SEQUENCE, INSERTED_SEQUENCE_LENGTH, INSERTED_SEQUENCE_MAPPINGS, HOMOLOGY, HOMOLOGY_LENGTH,
                        DUP_REPEAT_UNIT_REF_SPAN, DUP_SEQ_CIGARS, DUPLICATION_NUMBERS, DUP_ANNOTATIONS_IMPRECISE, DUP_IMPRECISE_AFFECTED_RANGE,
                        DUP_TAN_CONTRACTION_STRING, DUP_TAN_EXPANSION_STRING, DUP_ORIENTATIONS, INV33, INV55, EXTERNAL_CNV_CALLS,
                        CTG_GOOD_NONCANONICAL_MAPPING, LINK)
                        .sorted().collect(Collectors.toList()));
    }

    @Test(groups = "sv")
    public void testHeaderLines() {
        Assert.assertTrue(GATKSVVCFHeaderLines.getFilterLines().isEmpty());
        Assert.assertTrue(GATKSVVCFHeaderLines.getFormatLines().isEmpty());
        final Stream<String> infoHeaders = GATKSVVCFHeaderLines.getInfoLines().stream().map(VCFInfoHeaderLine::getID);
        final Stream<String> altAlleleHeaders = GATKSVVCFHeaderLines.getSymbAltAlleleLines().stream().map(VCFSimpleHeaderLine::getID);
        Assert.assertEquals(Stream.concat(infoHeaders, altAlleleHeaders).sorted().collect(Collectors.toList()),
                GATKSVVCFConstants.expectedHeaderLinesInVCF);
    }

    @DataProvider(name = "svVcfFiles")
    private Object[][] testDataForSVVCFFiles() throws IOException {
        final List<Object[]> data = new ArrayList<>(20);
        final List<Path> vcfFiles = Files.walk(Paths.get(toolsTestDir + "spark/sv/integration/outputs"))
                .filter(filePath -> filePath.toString().endsWith(".vcf") || filePath.toString().endsWith(".vcf.gz"))
                .collect(Collectors.toList());
        vcfFiles.addAll( Files.walk(Paths.get(toolsTestDir + "spark/sv/utils"))
                .filter(filePath -> filePath.toString().endsWith(".vcf") || filePath.toString().endsWith(".vcf.gz"))
                .collect(Collectors.toList()));
        vcfFiles.forEach(p -> data.add( new Object[]{p}));
        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "svVcfFiles", groups = "sv")
    public void checkTestVcfFiles(final Path svVCFFilePath) {
        try (final VCFFileReader reader = new VCFFileReader(svVCFFilePath.toFile(), false)) {
            final VCFHeader fileHeader = reader.getFileHeader();
            Assert.assertNotNull(fileHeader.getSequenceDictionary());
            Assert.assertTrue(fileHeader.getFilterLines().isEmpty());
            final List<String> refContigs = fileHeader.getContigLines().stream().map(VCFContigHeaderLine::getID).collect(Collectors.toList());
            Assert.assertFalse(refContigs.isEmpty());
            Assert.assertTrue(fileHeader.getFormatHeaderLines().isEmpty());
            final List<String> headerKeys = fileHeader.getIDHeaderLines().stream().map(VCFIDHeaderLine::getID).sorted().collect(Collectors.toList());
            Assert.assertTrue(headerKeys.remove(VCFConstants.END_KEY));
            Assert.assertTrue(headerKeys.removeAll(refContigs));
            Assert.assertEquals(headerKeys, GATKSVVCFConstants.expectedHeaderLinesInVCF);
        }
    }
}
