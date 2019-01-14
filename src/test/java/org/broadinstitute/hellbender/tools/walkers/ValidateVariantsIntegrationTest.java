package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

import static org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants.ValidationType.*;

public final class ValidateVariantsIntegrationTest extends CommandLineProgramTest {

    private static final String MITO_REF = toolsTestDir + "mutect/mito/Homo_sapiens_assembly38.mt_only.fasta";

    public String baseTestString(final boolean sharedFile, final String file, final boolean exclude, final ValidateVariants.ValidationType type) {

        return baseTestString(sharedFile, file, exclude, type, null, hg19_chr1_1M_Reference);
    }

    public String baseTestString(boolean sharedFile, String file, boolean exclude, ValidateVariants.ValidationType type, String region, String reference) {
        final String filePath = sharedFile ? file: getToolTestDataDir() + file;
        final String typeArgString = exclude ? " --validation-type-to-exclude " + type.name() : excludeValidationTypesButString(type);
        final String intervals = region == null ? "" : " -L " + region;
        final String referenceString = reference == null ? "" : " -R " + reference;

        return referenceString + intervals + " --variant " + filePath + typeArgString;
    }

    public String baseTestStringWithoutReference(boolean sharedFile, String file, boolean exclude, ValidateVariants.ValidationType type) {
        final String filePath = sharedFile ? file: getToolTestDataDir() + file;
        final String typeArgString = exclude ? " --validation-type-to-exclude " + type.name() : excludeValidationTypesButString(type);
        return " --variant " + filePath + typeArgString;
    }

    private static String excludeValidationTypesButString(ValidateVariants.ValidationType type) {
        if (type.equals(ALL)) {
            return "";
        }
        final StringBuilder sbuilder = new StringBuilder();
        for (final ValidateVariants.ValidationType t : CONCRETE_TYPES) {
            if (t != type) {
                sbuilder.append(" --validation-type-to-exclude ").append(t.toString());
            }
        }
        return sbuilder.toString();
    }

    @Test
    public void testGoodFile() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationExampleGood.vcf", false, ALL),
                Collections.emptyList()
        );

        spec.executeTest("test good file", this);
    }

    @Test
    public void testGoodFile2() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(true, hg19_chr1_1M_exampleVCF, false, ALL),
                Collections.emptyList()
        );

        spec.executeTest("test good file", this);
    }

    @Test
    public void testBadRefBase1() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationExampleBad.vcf", false, REF),
                0,
                UserException.FailsStrictValidation.class
        );

        spec.executeTest("test bad ref base #1", this);
    }

    @Test
    public void testMissingReference() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestStringWithoutReference(false, "validationExampleBad.vcf", false, REF),
                0,
                UserException.MissingReference.class
        );

        spec.executeTest("test bad ref base #1", this);
    }

    @Test
    public void testBadRefBase2() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationExampleBad2.vcf", false, REF),
                0,
                UserException.FailsStrictValidation.class
        );

        spec.executeTest("test bad ref base #2", this);
    }

    @Test
    public void testBadChrCount1() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationExampleBad.vcf", false, CHR_COUNTS),
                0,
                UserException.FailsStrictValidation.class
        );

        spec.executeTest("test bad chr counts #1", this);
    }

    @Test
    public void testBadChrCount2() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationExampleBad2.vcf", false, CHR_COUNTS),
                0,
                UserException.FailsStrictValidation.class
        );

        spec.executeTest("test bad chr counts #2", this);
    }

    @Test
    public void testBadID() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationExampleBadRSID.vcf", false, IDS) + " --dbsnp " + hg19_chr1_1M_dbSNP_modified,
                0,
                UserException.FailsStrictValidation.class
        );
        spec.executeTest("test bad RS ID", this);
    }

    @Test
    public void testBadID2_OKif_noDBSNPArgument() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationExampleBadRSID.vcf", false, IDS),
                Collections.emptyList()
        );
        spec.executeTest("test bad RS ID is OK if there's no dbsnp argument", this);
    }

    @Test
    public void testBadID2_OKif_notInDBSNP() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationExampleRSIDonPositionNotInDBSNP.vcf", false, IDS) + " --dbsnp " + hg19_chr1_1M_dbSNP_modified,
                Collections.emptyList()
        );
        spec.executeTest("test bad RS ID is OK when not in dbSNP even with dbsnp arg", this);
    }

    @Test
    public void testBadID2_okIfIDnotInDBSNP_withoutDbDNParg() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationExampleRSIDonPositionNotInDBSNP.vcf", false, IDS),
                Collections.emptyList()
        );
        spec.executeTest("test bad RS ID is OK when not in dbSNP when no dbsnp arg", this);
    }

    @Test
    public void testBadAllele() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
            baseTestString(false, "validationExampleBad.vcf", false, ALLELES),
            0,
            UserException.FailsStrictValidation.class
        );

        spec.executeTest("test bad alt allele", this);
    }

    @Test
    public void testBadAllele2() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
            baseTestString(false, "validationExampleBad3.vcf", false, REF),
            0,
            UserException.FailsStrictValidation.class
        );

        spec.executeTest("test bad ref allele in deletion", this);
    }

    @Test
    public void testComplexEventsDictError() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "complexEvents_incompatibleDict.vcf", false, ALL),
                0,
                UserException.IncompatibleSequenceDictionaries.class
        );

        spec.executeTest("test validating complex events", this);
    }

    @Test
    public void testNoValidation() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationExampleBad.vcf", true, ALL),
                Collections.emptyList()
        );

        spec.executeTest("test no validation", this);
    }

    @Test
    public void testComplexEvents() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "complexEvents.vcf", false, ALL),
                Collections.emptyList()
        );

        spec.executeTest("test validating complex events", this);
    }

    //The problem was that a gvcf contained a variant with <NON_REF> alt and no sample had that genotype.
    // This is not allowed in the vcf spec but OK in gvcf. The fix is to not fail when ALLELES is not checked.
    @Test(description = "Fixes '''bug''' reported in story https://www.pivotaltracker.com/story/show/68725164")
    public void testUnusedAlleleFix() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationUnusedAllelesBugFix.vcf", true, ALLELES,"1:1-739000",hg19_chr1_1M_Reference), Collections.emptyList());
        spec.executeTest("test unused allele bug fix", this);
    }

    //The problem was that a gvcf contained a variant with <NON_REF> alt and no sample had that genotype.
    // This is not allowed in the vcf spec but OK in gvcf. The fix is to not fail when ALLELES is not checked.
    @Test(description = "Checks '''bug''' reported in story https://www.pivotaltracker.com/story/show/68725164")
    public void testUnusedAlleleError() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "validationUnusedAllelesBugFix.vcf", false, ALLELES,"1:1-739000",hg19_chr1_1M_Reference),0, UserException.FailsStrictValidation.class);
        spec.executeTest("test unused allele bug fix", this);
    }

    @Test
    public void testGoodGvcf() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "NA12891.AS.chr20snippet.g.vcf", false, ALL, "20:10433000-10437000", b37_reference_20_21) + " -gvcf ",
                Collections.emptyList());
        spec.executeTest("tests correct gvcf", this);
    }

    @Test
    public void testGoodGvcfExcludingAlleles() throws IOException  {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "NA12891.AS.chr20snippet.g.vcf", true, ALLELES, "20:10433000-10437000", b37_reference_20_21) + " --gvcf ",
                Collections.emptyList());
        spec.executeTest("tests correct gvcf", this);
    }

    @Test
    public void testNoIntervalFullGenomeGVCF() throws IOException  {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "NA12891.AS.chr20snippet.g.vcf", true, REF, null, null) + " --gvcf ",
                0,
                UserException.class);
        spec.executeTest("tests gvcf that doesnt cover full genome but doesnt give an interval", this);
    }

    @Test
    public void testNoIntervalShortenedHeaderDict() throws IOException  {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "NA12891.AS.fullchr20_with_chr20_dict_only.g.vcf", true, REF, null, null) + " --gvcf ",
                Collections.emptyList());
        spec.executeTest("tests gvcf that has a header dictionary of just chr20 and doesn't give an interval", this);
    }

    @Test
    public void testBadGvcfMissingNON_REF() throws IOException  {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "NA12891.AS.chr20snippet.BAD_MISSING_NON_REF.g.vcf", true, ALLELES, "20:10433000-10437000", b37_reference_20_21) + " -gvcf ",
                0, UserException.class);
        spec.executeTest("tests capture of missing NON_REF allele", this);
    }

    //need to find file for this test
    @Test()
    public void testBadGvcfRegions() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "NA12891.AS.chr20snippet.missingrefblock.g.vcf", true, ALLELES, "20:10433000-10437000", b37_reference_20_21) + " -gvcf  ",
                0, UserException.class);
        spec.executeTest("tests capture of a gvcf missing a reference block", this);
    }

    @Test
    public void testBadGvcfOutOfOrder() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "badGVCF.outOfOrder.g.vcf", true, ALLELES, "chrM:1-1000", MITO_REF) + " -gvcf  ",
                0, UserException.class);
        spec.executeTest("tests capture of a gvcf missing a reference block", this);
    }

    @Test()
    public void testNonOverlappingRegions() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "NA12891.AS.chr20snippet_BAD_INCOMPLETE_REGION.g.vcf", true, ALLELES, "Y:4966254-4967190", b37_reference_20_21) + " -gvcf ",
                0, UserException.class);
        spec.executeTest("tests capture of non-complete region", this);
    }

    @Test
    public void testNonOverlappingRegionsBP_RESOLUTION() throws IOException {
        IntegrationTestSpec spec = new IntegrationTestSpec(
                baseTestString(false, "gvcf.basepairResolution.vcf", true, ALLELES, "20:10000000-10002158", b37_reference_20_21) + " -gvcf ",
                Collections.emptyList());
        spec.executeTest("tests capture of non-complete region, on BP_RESOLUTION gvcf", this);
    }
}
