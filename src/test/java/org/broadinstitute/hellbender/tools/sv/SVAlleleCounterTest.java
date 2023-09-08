package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SVAlleleCounterTest extends GATKBaseTest {

    @DataProvider(name = "testCounterData")
    public Object[][] testCounterData() {
        return new Object[][]{
                // Empty cases
                {
                        new Allele[]{},
                        new Allele[][]{},
                        new int[]{},
                        new double[]{},
                        0
                },
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{},
                        new int[]{0},
                        new double[]{Double.NaN},
                        0
                },
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{{}},
                        new int[]{0},
                        new double[]{Double.NaN},
                        0
                },

                // Haploid
                {
                        new Allele[]{},
                        new Allele[][]{
                                {Allele.REF_N}
                        },
                        new int[]{},
                        new double[]{},
                        1
                },
                // Non-ref allele but not defined in alts
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{
                                {Allele.SV_SIMPLE_DEL}
                        },
                        new int[]{0},
                        new double[]{0.},
                        1
                },
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{
                                {Allele.SV_SIMPLE_INS}
                        },
                        new int[]{1},
                        new double[]{1.},
                        1
                },
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{
                                {Allele.REF_N}
                        },
                        new int[]{0},
                        new double[]{0.},
                        1
                },
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{
                                {Allele.SV_SIMPLE_INS},
                                {Allele.REF_N}
                        },
                        new int[]{1},
                        new double[]{0.5},
                        2
                },
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{
                                {Allele.SV_SIMPLE_INS},
                                {Allele.SV_SIMPLE_INS}
                        },
                        new int[]{2},
                        new double[]{1.},
                        2
                },
                // No-call not counted toward AN
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{
                                {Allele.SV_SIMPLE_INS},
                                {Allele.NO_CALL}
                        },
                        new int[]{1},
                        new double[]{1.},
                        1
                },
                // Multi-allelic
                {
                        new Allele[]{Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP},
                        new Allele[][]{
                                {Allele.SV_SIMPLE_DEL},
                                {Allele.SV_SIMPLE_DUP}
                        },
                        new int[]{1, 1},
                        new double[]{0.5, 0.5},
                        2
                },
                {
                        new Allele[]{Allele.SV_SIMPLE_DEL, Allele.SV_SIMPLE_DUP},
                        new Allele[][]{
                                {Allele.SV_SIMPLE_DUP},
                                {Allele.REF_N}
                        },
                        new int[]{0, 1},
                        new double[]{0., 0.5},
                        2
                },

                // Diploid
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{
                                {Allele.REF_N, Allele.REF_N}
                        },
                        new int[]{0},
                        new double[]{0.},
                        2
                },
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{
                                {Allele.REF_N, Allele.SV_SIMPLE_INS}
                        },
                        new int[]{1},
                        new double[]{0.5},
                        2
                },
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{
                                {Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS}
                        },
                        new int[]{2},
                        new double[]{1.},
                        2
                },
                {
                        new Allele[]{Allele.SV_SIMPLE_INS},
                        new Allele[][]{
                                {Allele.REF_N, Allele.SV_SIMPLE_INS},
                                {Allele.SV_SIMPLE_INS, Allele.SV_SIMPLE_INS}
                        },
                        new int[]{3},
                        new double[]{0.75},
                        4
                },
        };
    }

    @Test(dataProvider= "testCounterData")
    public void testCounter(final Allele[] altAlleles,
                            final Allele[][] genotypeAllelesArr,
                            final int[] expectedAlleleCounts,
                            final double[] expectedAlleleFrequencies,
                            final int expectedAlleleNumber) {
        final List<Genotype> genotypesList = new ArrayList<>();
        for (int i = 0; i < genotypeAllelesArr.length; i++) {
            final List<Allele> alleles = Arrays.asList(genotypeAllelesArr[i]);
            genotypesList.add(new GenotypeBuilder(String.valueOf(i), alleles)
                    .attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, alleles.size())
                    .make());
        }

        final List<Allele> alleles = Arrays.asList(altAlleles);
        final SVAlleleCounter counter = new SVAlleleCounter(alleles, genotypesList);

        Assert.assertEquals(counter.getAlleles(), alleles);
        Assert.assertEquals(counter.getCounts(), expectedAlleleCounts);
        Assert.assertEquals(counter.getFrequencies(), expectedAlleleFrequencies);
        Assert.assertEquals(counter.getNumber(), expectedAlleleNumber);
    }

}