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


    private static final Stream<String> expectedAltAlleleHeaderKeysInVCF
            = Stream.of("INV", "DEL", "INS", "DUP", "DUP:INV", "CPX");
    private static final Stream<String> expectedInfoHeaderKeysInVCF
            = Stream.of("SVTYPE", "SVLEN", "MATEID",
            "CIPOS", "CIEND", "IMPRECISE", "READ_PAIR_SUPPORT", "SPLIT_READ_SUPPORT", "LINK",
            "CTG_NAMES", "TOTAL_MAPPINGS", "MAPPING_QUALITIES", "HQ_MAPPINGS", "ALIGN_LENGTHS", "MAX_ALIGN_LENGTH",
            "SEQ_ALT_HAPLOTYPE", "INSSEQ", "INSLEN", "INSSEQ_MAP", "HOMSEQ", "HOMLEN", "DUP_REPEAT_UNIT_REF_SPAN",
            "DUP_SEQ_CIGARS", "DUP_NUM", "DUP_ANNOTATIONS_IMPRECISE", "CONTRACTION", "EXPANSION", "DUP_ORIENTATIONS",
            "INV33", "INV55", "EXTERNAL_CNV_CALLS", "DUP_IMPRECISE_AFFECTED_RANGE", "CTG_GOOD_NONCANONICAL_MAPPING",
            "ALT_ARRANGEMENT", "SEGMENTS", "CPX_EVENT");
    private static final Stream<String> expectedFormatHeaderKeysInVCF
            = Stream.of("CN", "CNQ");
    private static final Stream<String> expectedFilterHeaderKeysInVCF
            = Stream.of("LOW_MQ", "SHORT_ALN");
    static final List<String> expectedHeaderKeysInVCF
            = Stream.of(expectedAltAlleleHeaderKeysInVCF, expectedInfoHeaderKeysInVCF, expectedFormatHeaderKeysInVCF,
                        expectedFilterHeaderKeysInVCF)
            .flatMap(i->i).sorted().collect(Collectors.toList());

    @Test(groups = "sv")
    public void testVCFConstants() {
        Assert.assertEquals(expectedHeaderKeysInVCF,
                Stream.of(SVTYPE, SVLEN, IMPRECISE, CIPOS, CIEND, BND_MATEID_STR, SYMB_ALT_ALLELE_INV, READ_PAIR_SUPPORT,
                        SPLIT_READ_SUPPORT, SYMB_ALT_ALLELE_DEL, SYMB_ALT_ALLELE_INS, SYMB_ALT_ALLELE_DUP, SYMB_ALT_ALLELE_INVDUP,
                        CONTIG_NAMES, TOTAL_MAPPINGS, MAPPING_QUALITIES, HQ_MAPPINGS, ALIGN_LENGTHS, MAX_ALIGN_LENGTH,
                        SEQ_ALT_HAPLOTYPE, INSERTED_SEQUENCE, INSERTED_SEQUENCE_LENGTH, INSERTED_SEQUENCE_MAPPINGS, HOMOLOGY, HOMOLOGY_LENGTH,
                        DUP_REPEAT_UNIT_REF_SPAN, DUP_SEQ_CIGARS, DUPLICATION_NUMBERS, DUP_ANNOTATIONS_IMPRECISE, DUP_IMPRECISE_AFFECTED_RANGE,
                        DUP_TAN_CONTRACTION_STRING, DUP_TAN_EXPANSION_STRING, DUP_ORIENTATIONS, INV33, INV55, EXTERNAL_CNV_CALLS,
                        CTG_GOOD_NONCANONICAL_MAPPING, LINK, COPY_NUMBER_FORMAT, COPY_NUMBER_QUALITY_FORMAT,
                        CPX_SV_SYB_ALT_ALLELE_STR, CPX_EVENT_ALT_ARRANGEMENTS, CPX_SV_REF_SEGMENTS, CPX_EVENT_KEY,
                        ASSEMBLY_BASED_VARIANT_MQ_FILTER_KEY, ASSEMBLY_BASED_VARIANT_ALN_LENGTH_FILTER_KEY)
                        .sorted().collect(Collectors.toList()));
    }

    @Test(groups = "sv")
    public void testHeaderLines() {

        final Stream<String> infoHeaders = GATKSVVCFHeaderLines.getInfoLines().stream().map(VCFInfoHeaderLine::getID);
        final Stream<String> altAlleleHeaders = GATKSVVCFHeaderLines.getSymbAltAlleleLines().stream().map(VCFSimpleHeaderLine::getID);
        final Stream<String> formatHeaders = GATKSVVCFHeaderLines.getFormatLines().stream().map((VCFCompoundHeaderLine::getID));
        final Stream<String> filterHeaders = GATKSVVCFHeaderLines.getFilterLines().stream().map((VCFSimpleHeaderLine::getID));

        Assert.assertEquals(
                Stream.of(infoHeaders, altAlleleHeaders, formatHeaders, filterHeaders).flatMap(i -> i).sorted().collect(Collectors.toList()),
                expectedHeaderKeysInVCF);
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

            final List<String> refContigs = fileHeader.getContigLines().stream().map(VCFContigHeaderLine::getID).collect(Collectors.toList());
            Assert.assertFalse(refContigs.isEmpty());

            final List<String> headerKeys = fileHeader.getIDHeaderLines().stream().map(VCFIDHeaderLine::getID).sorted().collect(Collectors.toList());
            Assert.assertTrue(headerKeys.remove(VCFConstants.END_KEY));
            Assert.assertTrue(headerKeys.removeAll(refContigs));
            Assert.assertEquals(headerKeys, expectedHeaderKeysInVCF);
        }
    }
}
