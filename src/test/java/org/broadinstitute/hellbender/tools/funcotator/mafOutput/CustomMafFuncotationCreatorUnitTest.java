package org.broadinstitute.hellbender.tools.funcotator.mafOutput;

import com.google.common.collect.ImmutableMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.TumorNormalPair;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class CustomMafFuncotationCreatorUnitTest extends GATKBaseTest {

    /**
     * The ground truth field is ordered by {@link CustomMafFuncotationCreator#COUNT_FIELD_NAMES}
     * @return array of objects to be used in tests.  Never {@code null}
     */
    @DataProvider
    public Object[][] provideVariantsAndPairs() {
        return new Object[][] {
                {new VariantContextBuilder("foo", "3", 1000000, 1000000, Arrays.asList(Allele.create("A", true), Allele.create("T", false)))
                        .genotypes(
                                Arrays.asList(
                                        new GenotypeBuilder("T1",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false)))
                                                .AD(new int[]{10,2})
                                                .attribute("AF", "0.166").make(),
                                        new GenotypeBuilder("N1",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false)))
                                                .AD(new int[]{20,0}).make()
                                        )
                        )
                        .make(),
                        Collections.singletonList(new TumorNormalPair("T1", "N1")),
                        ImmutableMap.of(Allele.create("T", false), Arrays.asList("2", "10", "0", "20", "0.166"))
                },{new VariantContextBuilder("foo", "3", 1000000, 1000000, Arrays.asList(Allele.create("A", true), Allele.create("T", false)))
                    .genotypes(
                            Arrays.asList(
                                    new GenotypeBuilder("T1",
                                            Arrays.asList(Allele.create("A", true), Allele.create("T", false)))
                                            .AD(new int[]{10,2})
                                            .attribute("AF", "0.166").make()
                                    )
                            )
                    .make(),
                    Collections.singletonList(new TumorNormalPair("T1", "")),
                    ImmutableMap.of(Allele.create("T", false), Arrays.asList("2", "10", "", "", "0.166"))
                },
                {new VariantContextBuilder("foo", "3", 1000000, 1000000, Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                        .genotypes(
                                Arrays.asList(
                                        new GenotypeBuilder("T1",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                                                .AD(new int[]{10,2,4})
                                                .attribute("AF", "0.166,0.32").make(),
                                        new GenotypeBuilder("N1",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                                                .AD(new int[]{20,0,1}).make()
                                )
                        )
                        .make(),
                        Collections.singletonList(new TumorNormalPair("T1", "N1")),
                        ImmutableMap.of(Allele.create("T", false), Arrays.asList("2", "10", "0", "20", "0.166"),
                                Allele.create("AT", false), Arrays.asList("4", "10", "1", "20", "0.32"))
                },

                // Test when we have no recognizable pairs that the fields are generated but empty.
                {new VariantContextBuilder("foo", "3", 1000000, 1000000, Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                        .genotypes(
                                Arrays.asList(
                                        new GenotypeBuilder("GERMLINE1",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                                                .AD(new int[]{10,2,4})
                                                .attribute("AF", "0.166,0.32").make(),
                                        new GenotypeBuilder("GERMLINE2",
                                                Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("AT", false)))
                                                .AD(new int[]{20,0,1}).make()
                                )
                        )
                        .make(),
                        Collections.emptyList(),
                        ImmutableMap.of(Allele.create("T", false), Collections.nCopies(CustomMafFuncotationCreator.COUNT_FIELD_NAMES.size(), ""),
                                Allele.create("AT", false), Collections.nCopies(CustomMafFuncotationCreator.COUNT_FIELD_NAMES.size(), ""))
                }

        };
    }

    @Test(dataProvider = "provideVariantsAndPairs")
    public void testCreation(final VariantContext variant, final List<TumorNormalPair> tnPair, final Map<Allele, List<String>> gtOrderedCustomFieldsByAllele) {
        final List<Funcotation> guess = CustomMafFuncotationCreator.createCustomMafCountFields(variant, tnPair);
        Assert.assertEquals(guess.size(), variant.getAlternateAlleles().size());
        for (final Allele allele: variant.getAlternateAlleles()) {
            final Funcotation funcotationOfInterest = guess.stream().filter(f -> f.getAltAllele().equals(allele)).findFirst().orElse(null);
            Assert.assertNotNull(funcotationOfInterest);
            final List<String> gtValues = gtOrderedCustomFieldsByAllele.getOrDefault(allele, Collections.emptyList());

            for (int i = 0; i < CustomMafFuncotationCreator.COUNT_FIELD_NAMES.size(); i++) {
                final String countField = CustomMafFuncotationCreator.COUNT_FIELD_NAMES.get(i);
                Assert.assertEquals(funcotationOfInterest.getField(countField), gtValues.get(i));
            }
        }
    }
}
