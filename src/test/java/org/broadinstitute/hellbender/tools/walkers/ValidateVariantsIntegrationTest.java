package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Collections;

import static org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants.ValidationType.*;

public final class ValidateVariantsIntegrationTest extends CommandLineProgramTest {

    public String baseTestString(final boolean sharedFile, final String file, final boolean exclude, final ValidateVariants.ValidationType type) {
        final String defaultRegion = "1:1-1000000";

        return baseTestString(sharedFile, file, exclude, type, defaultRegion, hg19_chr1_1M_Reference);
    }

    public String baseTestString(boolean sharedFile, String file, boolean exclude, ValidateVariants.ValidationType type, String region, String reference) {
        final String filePath = sharedFile ? file: getToolTestDataDir() + file;
        final String typeArgString = exclude ? " --validationTypeToExclude " + type.name() : excludeValidationTypesButString(type);
        final String intervals = "";//TODO enable this: " -L " + region;

        return "-R " + reference + intervals + " --variant " + filePath + typeArgString;
    }

    private static String excludeValidationTypesButString(ValidateVariants.ValidationType type) {
        if (type.equals(ALL)) {
            return "";
        }
        final StringBuilder sbuilder = new StringBuilder();
        for (final ValidateVariants.ValidationType t : CONCRETE_TYPES) {
            if (t != type) {
                sbuilder.append(" --validationTypeToExclude ").append(t.toString());
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
}
