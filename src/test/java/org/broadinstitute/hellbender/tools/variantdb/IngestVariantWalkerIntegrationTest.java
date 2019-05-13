package org.broadinstitute.hellbender.tools.variantdb;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.variantdb.IngestPetCreation;
import org.broadinstitute.hellbender.tools.variantdb.IngestVetExomeCreation;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

@Test(groups = {"variantcalling"})
public class IngestVariantWalkerIntegrationTest extends CommandLineProgramTest {
    private static String TEST_RESOURCES_PATH = "src/test/resources/org/broadinstitute/hellbender/tools/walkers/variantutils/ReblockGVCF/NA12878/";
    private static final File truthVET = new File(TEST_RESOURCES_PATH + "vet/vet.tsv");
    private static final File truthPET = new File(TEST_RESOURCES_PATH + "pet/pet.tsv");
    private static final File truthMETADATA = new File(TEST_RESOURCES_PATH + "metadata/metadata.tsv");
    private static final File inputVCF = new File(TEST_RESOURCES_PATH + "expected.NA12878.AS.chr20snippet.reblocked.g.vcf");
    private static final File inputINTERVAL_LIST = new File(TEST_RESOURCES_PATH + "expected.NA12878.interval_list");

    // TODO should there be a Missing test one too? // createMissingTSV
    private static final int PET_COL_COUNT = 3;

    private static final String SAMPLE_NAME = "rose";
    private static final String SAMPLE_ID = "1";
    private static final Allele ARef = Allele.create("A", true);
    private static final Allele C = Allele.create("C");
    private static final Allele NON_REF = Allele.NON_REF_ALLELE;
    private static List<Allele> variantBases = Lists.newArrayList(ARef, C);
    private static List<Allele> alleles = Lists.newArrayList(ARef, NON_REF);
    private static final int GQ30 = 30;
    private static final Genotype GENOTYPE = new GenotypeBuilder(SAMPLE_NAME).alleles(Arrays.asList(C, C)).DP(10).AD(new int[]{10,0}).GQ(GQ30).make();
    private static final Genotype REF_GENOTYPE = new GenotypeBuilder(SAMPLE_NAME).GQ(GQ30).make();


    final VariantContext VARIANT = new VariantContextBuilder(
            SAMPLE_NAME,
            "1",
            1,
            1,
            variantBases
    ).attribute(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, "AS_RAW_MQ")
            .attribute("AS_QUALapprox", "AS_QUALapprox")
            .attribute(GATKVCFConstants.AS_SB_TABLE_KEY, "AS_SB_TABLE")
            .attribute("AS_VarDP", "AS_VarDP")
            .genotypes(GENOTYPE)
            .make();

    final VariantContext VARIANT_MISSING_REQUIRED = new VariantContextBuilder(
            SAMPLE_NAME,
            "1",
            1,
            1,
            variantBases
    ).attribute(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, "AS_RAW_MQ")
            .genotypes(GENOTYPE)
            .make();

    final VariantContext REF_VARIANT = new VariantContextBuilder(
            SAMPLE_NAME,
            "1",
            1,
            1,
            alleles
    ).genotypes(REF_GENOTYPE)
            .attribute(VCFConstants.END_KEY, 1)
            .make();

    //the variant walker -- tsv creation --are the tsvs valid -- 3 cols
     // √ check the tsvs are valid: PET is 3 columns & expected length
    // do I get the expected truth values?
    // data provider with different intervals (missing and non-missing)
    // test for removing GQ60 or not remove anything
    // truth of nothing removed and something removed and negative tests -- 2 truth files, 4 tests
    // add var for ""

    @Test(dataProvider = "inputVCF")
    public void testBlahVariantWalker(
            final File truthPET,
            final File truthVET,
            final File inputVCF,
            final File inputIntervalList,
            final Optional<String> missingArg,
            final Boolean successBool
    ) {
        final File vetDirectory = createTempDir("vet");
        final File petDirectory = createTempDir("pet");
        final File metadataWritten = createTempFile("metadata", ".tsv");
        final List<String> args = Arrays.asList(
                "-V", inputVCF.getAbsolutePath(), // the original gVCF for ingest
                "-VO", vetDirectory.getAbsolutePath(), // Path to the directory where the variants table should be written
                "-PO", petDirectory.getAbsolutePath(), // Path to the directory where the positions expanded table should be written
                "-SMO", metadataWritten.getAbsolutePath(), // Path to where the positions expanded table should be written
                missingArg.isPresent()? "-IG" : "",
                missingArg.isPresent()? missingArg.get() : "",
                "-L", inputIntervalList.getAbsolutePath()); // TODO do this with missing and non-missing
        // TODO what is the successBool expecting?

        runCommandLine(args);

        // verify the output pet tsv
        File[] petFiles = petDirectory.listFiles();
        File petFile = petFiles[0];
        Assert.assertTrue(petFile.isFile());
        try {
            long lineCount = Files.lines(Paths.get(petFile.getAbsolutePath())).count();
            long truthLineCount = Files.lines(Paths.get(truthPET.getAbsolutePath())).count();

            if (missingArg == null) {
                Assert.assertEquals(lineCount, truthLineCount);
            } else {
                Assert.assertEquals(lineCount, truthLineCount); // truthline count should be way less if missing gq60?
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        try (Stream<String> petStream = Files.lines(Paths.get(petFile.getAbsolutePath()))) {
            final String[] petRowValues = petStream.findAny().get().split("\t");
            Assert.assertEquals(petRowValues.length, PET_COL_COUNT);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // verify the output vet tsv
        File[] vetFiles = vetDirectory.listFiles();
        File vetFile = vetFiles[0];
        Assert.assertTrue(vetFile.isFile());
        try {
            long lineCount = Files.lines(Paths.get(vetFile.getAbsolutePath())).count();
            long truthLineCount = Files.lines(Paths.get(truthVET.getAbsolutePath())).count();
            Assert.assertEquals(lineCount, truthLineCount);
        } catch (IOException e) {
            e.printStackTrace();
        }
        try (Stream<String> stream = Files.lines(Paths.get(vetFile.getAbsolutePath()))) {
            // TODO do I need to bother to assert any actual values? stream.forEach(System.out::println);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @DataProvider(name = "inputVCF")
    public Object[][] inputVCF() {
        return new Object[][]{
                // truthPET, truthVET, inputVCF, input interval list, missingArg, successBool
                //{truthPET, truthVET, inputVCF, inputINTERVAL_LIST, null, true},
                //{truthPET, truthVET, inputVCF, inputINTERVAL_LIST, Optional.ofNullable("SIXTY"), true},
                {truthPET, truthVET, inputVCF, inputINTERVAL_LIST, Optional.ofNullable("THIRTY"), true},
               // {truthPET, truthVET, inputVCF, inputINTERVAL_LIST, null, false},
                //{truthPET, truthVET, inputVCF, inputINTERVAL_LIST, IngestPetCreation.GQStateEnum.SIXTY, false},
        };
    }

    @Test
    public void testPetCreation() throws Exception {
        // start with a variant
        final List<List<String>> positionExpandedVariantRows = IngestPetCreation.createPositionRows(1, 1, VARIANT, SAMPLE_NAME);

        final List<String> positionExpandedVariantRowTruth = new ArrayList<>();
        positionExpandedVariantRowTruth.add(String.valueOf(1));
        positionExpandedVariantRowTruth.add(SAMPLE_NAME);
        positionExpandedVariantRowTruth.add(IngestPetCreation.GQStateEnum.VARIANT.value);

        final List<List<String>> positionExpandedRowsTruth = new ArrayList<>();
        positionExpandedRowsTruth.add(positionExpandedVariantRowTruth);

        Assert.assertEquals(positionExpandedVariantRows, positionExpandedRowsTruth);

        // TODO then add some GQ30 or whatever block
        final List<List<String>> positionExpandedRefRows = IngestPetCreation.createPositionRows(1, 1, REF_VARIANT, SAMPLE_NAME);
        final List<String> positionExpandedRefRowTruth = new ArrayList<>();
        positionExpandedRefRowTruth.add(String.valueOf(1));
        positionExpandedRefRowTruth.add(SAMPLE_NAME);
        positionExpandedRefRowTruth.add(IngestPetCreation.GQStateEnum.THIRTY.value);

        final List<List<String>> positionExpandedRefRowsTruth = new ArrayList<>();
        positionExpandedRefRowsTruth.add(positionExpandedRefRowTruth);
        System.out.print(REF_VARIANT.getAlternateAlleles());
        System.out.print(REF_VARIANT.getAlternateAllele(0));
        System.out.print(REF_VARIANT.getAlternateAlleles().size());

        Assert.assertTrue(REF_VARIANT.getAttribute(VCFConstants.END_KEY) != null, "threee");

        Assert.assertEquals(positionExpandedRefRows, positionExpandedRefRowsTruth);


        // IngestPetCreation.createSpanDelRows();

    }

    @Test
    public void testVetCreation() throws Exception {
        // final List<String> variantRow = IngestVetExomeCreation.createVariantRow(VARIANT);
        final List<String> variantRow = IngestVetExomeCreation.createVariantRow(1, VARIANT, SAMPLE_ID);
        final List<String> variantRowTruth = new ArrayList<>();
        variantRowTruth.add(String.valueOf(VARIANT.getStart()));
        variantRowTruth.add(VARIANT.getReference().getBaseString());
        variantRowTruth.add(VARIANT.getAlternateAlleles().get(0).toString());
        variantRowTruth.add(String.valueOf(VARIANT.getAttribute(GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY)));
        variantRowTruth.add("");
        variantRowTruth.add(String.valueOf(VARIANT.getAttribute("AS_QUALapprox")));
        variantRowTruth.add("");
        //variantRowTruth.add(String.valueOf(VARIANT.getAttribute(GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY))); // TODO should not be null
        variantRowTruth.add(String.valueOf(VARIANT.getAttribute(GATKVCFConstants.AS_SB_TABLE_KEY)));
        variantRowTruth.add(String.valueOf(VARIANT.getAttribute("AS_VarDP")));
        variantRowTruth.add(VARIANT.getSource());
        ArrayList<Integer> allele_indices = new ArrayList<Integer>();
        for (Allele allele : VARIANT.getGenotype(0).getAlleles()){
            allele_indices.add(GATKVariantContextUtils.indexOfAllele(VARIANT, allele, true, true, true  ));
        }
        variantRowTruth.add(String.valueOf(VARIANT.getGenotype(0).isPhased() ? StringUtils.join(allele_indices, VCFConstants.PHASED) : StringUtils.join(allele_indices, VCFConstants.UNPHASED)));
        variantRowTruth.add(Arrays.toString(VARIANT.getGenotype(0).getAD()));  // TODO is this a bug?
        variantRowTruth.add(String.valueOf(VARIANT.getGenotype(0).getDP()));
        variantRowTruth.add(String.valueOf(VARIANT.getGenotype(0).getGQ()));
        variantRowTruth.add("");
        variantRowTruth.add("");
        variantRowTruth.add("");
        System.out.print(variantRow);
        System.out.print(variantRowTruth);
        Assert.assertEquals(variantRow, variantRowTruth);

        // make sure if there is missing req field that we throw an error
        Assert.assertThrows(() -> IngestVetExomeCreation.createVariantRow(1, VARIANT_MISSING_REQUIRED, SAMPLE_ID));
    }
}
