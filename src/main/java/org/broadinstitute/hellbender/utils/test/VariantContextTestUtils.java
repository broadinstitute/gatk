package org.broadinstitute.hellbender.utils.test;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;

import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public final class VariantContextTestUtils {

    private VariantContextTestUtils() {}

    private static void assertAttributeEquals(final String key, final Object actual, final Object expected) {
        final Object notationCorrectedActual = normalizeScientificNotation(actual);
        final Object notationCorrectedExpected = normalizeScientificNotation(expected);
        if (notationCorrectedExpected instanceof Double && notationCorrectedActual instanceof Double) {
            // must be very tolerant because doubles are being rounded to 2 sig figs
            BaseTest.assertEqualsDoubleSmart((Double) notationCorrectedActual, (Double) notationCorrectedExpected, 1e-2, "Attribute " + key);
        } else {
            Assert.assertEquals(notationCorrectedActual, notationCorrectedExpected, "Attribute " + key);
        }
    }

    /**
     * Normalizes the representation of Strings containing doubles.
     * This is necessary to deal with the fact that variant context attributes are deserialized from vcf as Strings
     * instead of their original type.  Some versions of gatk3 output double attributes in scientific notation while gatk4
     * doesn't do so.
     *
     * @param attribute an attribute to attempt to normalize
     * @return if attribute is a String, try to parse it as a Double and return that value, else return the original attribute
     */
    @VisibleForTesting
    static Object normalizeScientificNotation(final Object attribute){
        if (attribute instanceof String){
            try {
                return Double.parseDouble((String) attribute);
            } catch ( final NumberFormatException e) {
                return attribute;
            }
        }
        return attribute;
    }

    public static void assertGenotypesAreEqual(final Genotype actual, final Genotype expected) {
        Assert.assertEquals(actual.getSampleName(), expected.getSampleName(), "Genotype names");
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "Genotype alleles");
        Assert.assertEquals(actual.getGenotypeString(), expected.getGenotypeString(), "Genotype string");
        Assert.assertEquals(actual.getType(), expected.getType(), "Genotype type");

        // filters are the same
        Assert.assertEquals(actual.getFilters(), expected.getFilters(), "Genotype fields");
        Assert.assertEquals(actual.isFiltered(), expected.isFiltered(), "Genotype isFiltered");

        // inline attributes
        Assert.assertEquals(actual.hasDP(), expected.hasDP(), "Genotype hasDP");
        Assert.assertEquals(actual.getDP(), expected.getDP(), "Genotype dp");
        Assert.assertEquals(actual.hasAD(), expected.hasAD(), "Genotype hasAD");
        Assert.assertTrue(Arrays.equals(actual.getAD(), expected.getAD()));
        Assert.assertEquals(actual.hasGQ(), expected.hasGQ(), "Genotype hasGQ");
        Assert.assertEquals(actual.getGQ(), expected.getGQ(), "Genotype gq");
        Assert.assertEquals(actual.hasPL(), expected.hasPL(), "Genotype hasPL");
        Assert.assertTrue(Arrays.equals(actual.getPL(), expected.getPL()));

        Assert.assertEquals(actual.hasLikelihoods(), expected.hasLikelihoods(), "Genotype haslikelihoods");
        Assert.assertEquals(actual.getLikelihoodsString(), expected.getLikelihoodsString(), "Genotype getlikelihoodsString");
        Assert.assertEquals(actual.getLikelihoods(), expected.getLikelihoods(), "Genotype getLikelihoods");

        Assert.assertEquals(actual.getGQ(), expected.getGQ(), "Genotype phredScaledQual");
        assertAttributesEquals(actual.getExtendedAttributes(), expected.getExtendedAttributes());
        Assert.assertEquals(actual.isPhased(), expected.isPhased(), "Genotype isPhased");
        Assert.assertEquals(actual.getPloidy(), expected.getPloidy(), "Genotype getPloidy");
    }

    @SuppressWarnings("unchecked")
    private static void assertAttributesEquals(final Map<String, Object> actual, final Map<String, Object> expected) {
        final Set<String> expectedKeys = new LinkedHashSet<>(expected.keySet());

        for ( final Map.Entry<String, Object> act : actual.entrySet() ) {
            final Object actualValue = act.getValue();
            if ( expected.containsKey(act.getKey()) && expected.get(act.getKey()) != null ) {
                final Object expectedValue = expected.get(act.getKey());
                if ( expectedValue instanceof List) {
                    final List<Object> expectedList = (List<Object>)expectedValue;
                    Assert.assertTrue(actualValue instanceof List, act.getKey() + " should be a list but isn't");
                    final List<Object> actualList = (List<Object>)actualValue;
                    Assert.assertEquals(actualList.size(), expectedList.size(), act.getKey() + " size");
                    for ( int i = 0; i < expectedList.size(); i++ ) {
                        assertAttributeEquals(act.getKey(), actualList.get(i), expectedList.get(i));
                    }
                } else {
                    assertAttributeEquals(act.getKey(), actualValue, expectedValue);
                }
            } else {
                // it's ok to have a binding in x -> null that's absent in y
                Assert.assertNull(actualValue, act.getKey() + " present in one but not in the other");
            }
            expectedKeys.remove(act.getKey());
        }

        // now expectedKeys contains only the keys found in expected but not in actual,
        // and they must all be null
        for ( final String missingExpected : expectedKeys ) {
            final Object value = expected.get(missingExpected);
            Assert.assertTrue(isMissing(value), "Attribute " + missingExpected + " missing in one but not in other" );
        }
    }

    private static boolean isMissing(final Object value) {
        if ( value == null ) return true;
        else if ( value.equals(VCFConstants.MISSING_VALUE_v4) ) return true;
        else if ( value instanceof List ) {
            // handles the case where all elements are null or the list is empty
            for ( final Object elt : (List)value)
                if ( elt != null )
                    return false;
            return true;
        } else
            return false;
    }

    /**
     * Validates that the given lists have variant
     * context that correspond to the same variants in the same order.
     * Compares VariantContext by comparing toStringDecodeGenotypes
     */
    public static void assertEqualVariants(final List<VariantContext> v1, final List<VariantContext> v2) {
        Utils.nonNull(v1, "v1");
        Utils.nonNull(v2, "v2");
        if (v1.size() != v2.size()){
            throw new AssertionError("different sizes " + v1.size()+ " vs " + v2.size());
        }
        for (int i = 0; i < v1.size(); i++) {
            if (! v1.get(i).toStringDecodeGenotypes().equals(v2.get(i).toStringDecodeGenotypes())){
                throw new AssertionError("different element (compared by toStringDecodeGenotypes) " + i + "\n" + v1.get(i) + "\n" + v2.get(i) );
            }
        }
    }

    public static void assertVariantContextsAreEqual(final VariantContext actual, final VariantContext expected, final List<String> attributesToIgnore ) {
        Assert.assertNotNull(actual, "VariantContext expected not null");
        Assert.assertEquals(actual.getContig(), expected.getContig(), "chr");
        Assert.assertEquals(actual.getStart(), expected.getStart(), "start");
        Assert.assertEquals(actual.getEnd(), expected.getEnd(), "end");
        Assert.assertEquals(actual.getID(), expected.getID(), "id");
        Assert.assertEquals(actual.getAlleles(), expected.getAlleles(), "alleles for " + expected + " vs " + actual);

        assertAttributesEquals(filterIgnoredAttributes(actual.getAttributes(), attributesToIgnore),
                               filterIgnoredAttributes(expected.getAttributes(), attributesToIgnore));

        Assert.assertEquals(actual.filtersWereApplied(), expected.filtersWereApplied(), "filtersWereApplied");
        Assert.assertEquals(actual.isFiltered(), expected.isFiltered(), "isFiltered");
        BaseTest.assertEqualsSet(actual.getFilters(), expected.getFilters(), "filters");
        BaseTest.assertEqualsDoubleSmart(actual.getPhredScaledQual(), expected.getPhredScaledQual());

        assertVariantContextsHaveSameGenotypes(actual, expected);
    }

    private static Map<String, Object> filterIgnoredAttributes(final Map<String,Object> attributes, final List<String> attributesToIgnore){
        return attributes.entrySet().stream()
                .filter(p -> !attributesToIgnore.contains(p.getKey()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }

    public static void assertVariantContextsHaveSameGenotypes(final VariantContext actual, final VariantContext expected) {
        Assert.assertEquals(actual.hasGenotypes(), expected.hasGenotypes(), "hasGenotypes");
        if ( expected.hasGenotypes() ) {
            BaseTest.assertEqualsSet(actual.getSampleNames(), expected.getSampleNames(), "sample names set");
            Assert.assertEquals(actual.getSampleNamesOrderedByName(), expected.getSampleNamesOrderedByName(), "sample names");
            final Set<String> samples = expected.getSampleNames();
            for ( final String sample : samples ) {
                assertGenotypesAreEqual(actual.getGenotype(sample), expected.getGenotype(sample));
            }
        }
    }
}
