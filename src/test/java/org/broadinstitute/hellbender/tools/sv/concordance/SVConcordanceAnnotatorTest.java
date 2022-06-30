package org.broadinstitute.hellbender.tools.sv.concordance;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVTestUtils;
import org.testng.Assert;
import org.testng.TestException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.vcf.GenotypeConcordanceStates;

import java.util.Arrays;
import java.util.Collections;

public class SVConcordanceAnnotatorTest {

    @DataProvider(name = "testMatchCnvData")
    public Object[][] testMatchCnvData() {
        return new Object[][]{

                // Note we only match on CN of truth genotypes with respect to the val expected CN

                ////////////////////////////////////////
                // CNV - Ploidy 0
                ////////////////////////////////////////

                // truth equal
                {
                        0,                                              // eval expected copy number (ie ploidy)
                        0,                                              // eval copy number
                        0,                                              // truth copy number
                        true                                            // expected result
                },

                // truth non-equal
                {
                        0,
                        0,
                        1,
                        false
                },
                {
                        0,
                        0,
                        2,
                        false
                },
                {
                        0,
                        0,
                        3,
                        false
                },

                ////////////////////////////////////////
                // CNV - eval ploidy 1
                ////////////////////////////////////////

                // truth ref
                {
                        1,
                        1,
                        1,
                        true
                },

                // truth var
                {
                        1,
                        1,
                        0,
                        false
                },
                {
                        1,
                        1,
                        2,
                        false
                },
                // same var as eval
                {
                        1,
                        0,
                        0,
                        true
                },
                {
                        1,
                        2,
                        2,
                        true
                },
                // different CN as eval and therefore false
                {
                        1,
                        2,
                        3,
                        false
                },
                {
                        1,
                        2,
                        1,
                        false
                },
                {
                        2,
                        1,
                        0,
                        false
                },

                ////////////////////////////////////////
                // CNV - eval ploidy 2
                ////////////////////////////////////////

                // truth hom ref
                {
                        2,
                        2,
                        2,
                        true
                },

                // true - same types
                {
                        2,
                        1,
                        1,
                        true
                },
                {
                        2,
                        3,
                        3,
                        true
                },
                // false - different types or different CN
                {
                        2,
                        1,
                        3,
                        false
                },
                {
                        2,
                        1,
                        0,
                        false
                },
                {
                        2,
                        3,
                        1,
                        false
                },
                {
                        2,
                        3,
                        4,
                        false
                }
        };
    }

    @Test(dataProvider= "testMatchCnvData")
    public void testMatchCnvNonNull(
            final int evalExpectedCopyNumber,
            final int evalCopyNumber,
            final int truthCopyNumber,
            final boolean expected) {

        final SVCallRecord evalRecord = SVTestUtils.newCallRecordWithAlleles(
                Collections.nCopies(evalExpectedCopyNumber, Allele.NO_CALL),
                SVTestUtils.getCNVAlleles(StructuralVariantType.CNV),
                StructuralVariantType.CNV,
                evalExpectedCopyNumber,
                evalCopyNumber
        );

        final String sample = "sample";
        final int truthExpectedCopyNumber = 2;  // this shouldn't matter
        final SVCallRecord truthRecord = SVTestUtils.newCallRecordWithAlleles(
                Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL),
                SVTestUtils.getCNVAlleles(StructuralVariantType.CNV),
                StructuralVariantType.CNV,
                truthExpectedCopyNumber,
                truthCopyNumber
        );

        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator(false);
        final boolean actual = collapser.cnvGenotypesMatch(sample, evalRecord, truthRecord);
        Assert.assertEquals(actual, expected);
    }

    @Test
    public void testMatchCnvNull() {
        final String sample = "sample";
        final SVCallRecord evalRecord = SVTestUtils.newCallRecordWithAlleles(
                Arrays.asList(Allele.NO_CALL, Allele.NO_CALL),
                SVTestUtils.getCNVAlleles(StructuralVariantType.CNV),
                StructuralVariantType.CNV,
                2,
                2
        );

        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator(false);

        // Null truth results in null call
        Assert.assertNull(collapser.cnvGenotypesMatch(sample, evalRecord, null));
    }

    @DataProvider(name = "testGetStateFromGenotypeData")
    public Object[][] testGetStateFromGenotypeData() {
        return new Object[][]{

                ////////////////////////////////////////
                // Minimal edge cases
                ////////////////////////////////////////

                // no genotype alleles but not a null genotype
                { new Allele[]{}, GenotypeConcordanceStates.TruthState.NO_CALL },

                ////////////////////////////////////////
                // Haploid
                ////////////////////////////////////////

                // no-call
                { new Allele[]{ Allele.NO_CALL }, GenotypeConcordanceStates.TruthState.NO_CALL },
                // hom ref
                { new Allele[]{ Allele.REF_N }, GenotypeConcordanceStates.TruthState.HOM_REF },
                // hom var
                { new Allele[]{ Allele.SV_SIMPLE_INS }, GenotypeConcordanceStates.TruthState.HOM_VAR1 },

                ////////////////////////////////////////
                // Diploid
                ////////////////////////////////////////

                // hom ref
                { new Allele[]{ Allele.REF_N, Allele.REF_N }, GenotypeConcordanceStates.TruthState.HOM_REF },
                // hom var
                { new Allele[]{ Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS }, GenotypeConcordanceStates.TruthState.HOM_VAR1 },
                // het
                { new Allele[]{ Allele.REF_N, Allele.SV_SIMPLE_INS }, GenotypeConcordanceStates.TruthState.HET_REF_VAR1 },
                // mixed cases
                { new Allele[]{ Allele.REF_N, Allele.NO_CALL }, GenotypeConcordanceStates.TruthState.IS_MIXED },
                { new Allele[]{ Allele.NO_CALL, Allele.SV_SIMPLE_INS }, GenotypeConcordanceStates.TruthState.IS_MIXED },
        };
    }

    @Test(dataProvider= "testGetStateFromGenotypeData")
    public void testGetStateFromGenotype(final Allele[] genotypeAlleles,
                                         final GenotypeConcordanceStates.TruthState expectedTruthState) {
        final Genotype genotype = alleleArrayToGenotype(genotypeAlleles, null, null);
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator(false);
        final GenotypeConcordanceStates.TruthState actualTruthState = collapser.getTruthState(genotype);
        Assert.assertEquals(actualTruthState, expectedTruthState);

        // Equivalent call states
        final GenotypeConcordanceStates.CallState expectedCallState;
        if (expectedTruthState == GenotypeConcordanceStates.TruthState.HOM_REF) {
            expectedCallState = GenotypeConcordanceStates.CallState.HOM_REF;
        } else if (expectedTruthState == GenotypeConcordanceStates.TruthState.HET_REF_VAR1) {
            expectedCallState = GenotypeConcordanceStates.CallState.HET_REF_VAR1;
        } else if (expectedTruthState == GenotypeConcordanceStates.TruthState.HOM_VAR1) {
            expectedCallState = GenotypeConcordanceStates.CallState.HOM_VAR1;
        } else if (expectedTruthState == GenotypeConcordanceStates.TruthState.NO_CALL) {
            expectedCallState = GenotypeConcordanceStates.CallState.NO_CALL;
        } else if (expectedTruthState == GenotypeConcordanceStates.TruthState.IS_MIXED) {
            expectedCallState = GenotypeConcordanceStates.CallState.IS_MIXED;
        } else {
            throw new TestException("Unexpected truth state: " + expectedTruthState.name());
        }
        final GenotypeConcordanceStates.CallState actualCallState = collapser.getEvalState(genotype);
        Assert.assertEquals(actualCallState, expectedCallState);
    }

    @Test
    public void testGetTruthNullGenotypeState() {
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator(false);
        final GenotypeConcordanceStates.TruthState actualTruthState = collapser.getTruthState(null);
        Assert.assertEquals(actualTruthState, GenotypeConcordanceStates.TruthState.HOM_REF);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testGetNullEvalGenotypeState() {
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator(false);
        collapser.getEvalState(null);
    }

    private Genotype alleleArrayToGenotype(final Allele[] allelesArr,
                                           final Integer expectedCopyNumber,
                                           final Integer copyNumber) {
        final GenotypeBuilder builder = new GenotypeBuilder("");
        if (expectedCopyNumber != null) {
            builder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, expectedCopyNumber);
        }
        if (copyNumber != null) {
            builder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, copyNumber);
        }
        return builder.alleles(Arrays.asList(allelesArr)).make();
    }
}
