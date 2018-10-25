package org.broadinstitute.hellbender.testutils;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.io.output.NullOutputStream;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;

// This should be:
//import org.apache.logging.log4j.LogManager;
//import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public final class VariantContextTestUtils {

    private VariantContextTestUtils() {}

    /** Standard Logger.  */
    protected static final Logger logger = LogManager.getLogger(VariantContextTestUtils.class);

    /**
     * Reads an entire VCF into memory, returning both its VCFHeader and all VariantContext records in
     * the vcf. Supports both local files and NIO-supported remote filesystems such as GCS.
     *
     * For unit/integration testing purposes only! Do not call this method from actual tools!
     *
     * @param vcfPath path or URI to a VCF, as a String
     * @return A Pair with the VCFHeader as the first element, and a List of all VariantContexts from the VCF
     *         as the second element
     */
    public static Pair<VCFHeader, List<VariantContext>> readEntireVCFIntoMemory(final String vcfPath) {
        Utils.nonNull(vcfPath);

        try ( final FeatureDataSource<VariantContext> vcfReader = new FeatureDataSource<>(vcfPath) ) {
            final Object header = vcfReader.getHeader();
            if ( ! (header instanceof VCFHeader) ) {
                throw new IllegalArgumentException(vcfPath + " does not have a valid VCF header");
            }

            final List<VariantContext> vcfRecords = new ArrayList<>();
            for ( final VariantContext vcfRecord : vcfReader ) {
                vcfRecords.add(vcfRecord);
            }

            return Pair.of((VCFHeader)header, vcfRecords);
        }
    }

    private static void assertAttributeEquals(final String key, final Object actual, final Object expected) {
        final Object notationCorrectedActual = normalizeScientificNotation(actual);
        final Object notationCorrectedExpected = normalizeScientificNotation(expected);
        if (notationCorrectedExpected instanceof Double && notationCorrectedActual instanceof Double) {
            // must be very tolerant because doubles are being rounded to 2 sig figs
            BaseTest.assertEqualsDoubleSmart((Double) notationCorrectedActual, (Double) notationCorrectedExpected, 1e-2, "Attribute " + key);
        } else if (actual instanceof Integer || expected instanceof Integer) {
            Object actualNormalized = normalizeToInteger(actual);
            Object expectedNormalized = normalizeToInteger(expected);
            Assert.assertEquals(actualNormalized, expectedNormalized, "Attribute " + key);
        } else {
            Assert.assertEquals(notationCorrectedActual, notationCorrectedExpected, "Attribute " + key);
        }
    }

    /**
     * Attempt to convert a String containing a signed integer (no separators) to an integer. If the attribute is not
     * a String, or does not contain an integer, the original object is returned.
     *
     * @param attribute
     * @return An Integer representing the value in the attribute if it contains a parseable integer, otherwise the
     * original attribute.
     */
    @VisibleForTesting
    static Object normalizeToInteger(final Object attribute) {
        if (attribute instanceof String) {
            try {
                return Integer.parseInt((String) attribute);
            } catch ( final NumberFormatException e) {
                return attribute;
            }
        }
        return attribute;
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
                if (((String) attribute).contains("|")) {
                    // If the attribute is an allele specific attribute separated by '|', then we want to remap
                    // each of its contained values (which could be comma separated lists) separately
                    String[] split = ((String) attribute).split("\\|",-1);
                    return Arrays.stream(split).map(
                            s -> {return Arrays.stream(s.split(",",-1))
                                    .map(d -> {if (d.equals("")) return d;
                                    else return Double.toString(Double.parseDouble(d));})
                                    .collect(Collectors.joining(","));})
                            .collect(Collectors.joining("|"));
                } else {
                    return Double.parseDouble((String) attribute);
                }
            } catch ( final NumberFormatException e) {
                return attribute;
            }
        }
        return attribute;
    }

    /**
     * Method which reorders the AltAlleles of a VariantContext so that they are in alphabetical order.
     *
     * It also reorders all annotation fields and the PL,AD, and SAC fields of each sample genotype in the variant context
     * to be consistent with the new ordering
     *
     * NOTE: this method relies on the correctness of AlleleSubsettingUtils.subsetAlleles() for reordering alleles
     *
     * @param vc        Variant context to reorder
     * @param header
     * @return
     */
    public static VariantContext sortAlleles(final VariantContext vc, final VCFHeader header){
        final List<Allele> originalAltAlleles = vc.getAlternateAlleles();
        final List<Allele> sortedAltAlleles = originalAltAlleles.stream().sorted().collect(Collectors.toList());
        final List<Allele> sortedAlleles = new ArrayList<>(vc.getNAlleles());

        sortedAlleles.add(vc.getReference());
        sortedAlleles.addAll(sortedAltAlleles);

        final VariantContextBuilder result = new VariantContextBuilder(vc);
        result.alleles(sortedAlleles);

        GenotypesContext newGT = AlleleSubsettingUtils.subsetAlleles(vc.getGenotypes(),2,vc.getAlleles(),sortedAlleles,
                                                                     GenotypeAssignmentMethod.SET_TO_NO_CALL, vc.getAttributeAsInt(VCFConstants.DEPTH_KEY,0));

        // Asserting that the new genotypes were calculated properly in case AlleleSubsettingUtils behavior changes
        if (newGT.getSampleNames().size() != vc.getGenotypes().size()) throw new IllegalStateException("Sorting this variant context resulted in a different number of genotype alleles, check that AlleleSubsettingUtils still supports reordering:" + vc.toString());
        for (int i =0; i<newGT.size(); i++){
            if (vc.getGenotype(i).hasAD()) {
                if (newGT.get(i).getAD().length != vc.getGenotype(i).getAD().length) throw new IllegalStateException("Sorting this variant context resulted in a different number of genotype alleles, check that AlleleSubsettingUtils still supports reordering:" + vc.toString());
            }
            if (vc.getGenotype(i).hasPL()) {
                if (newGT.get(i).getPL().length != vc.getGenotype(i).getPL().length) throw new IllegalStateException("Sorting this variant context resulted in a different number of genotype alleles, check that AlleleSubsettingUtils still supports reordering:" + vc.toString());
            }
        }

        final HashMap<String, Object> newAttributes = new HashMap<>(vc.getAttributes());
        for (Map.Entry<String, Object> entry : newAttributes.entrySet()) {
            VCFHeaderLineCount type = header.hasInfoLine(entry.getKey())?header.getInfoHeaderLine(entry.getKey()).getCountType():VCFHeaderLineCount.UNBOUNDED;
            int ploidy = vc.getGenotypes().getMaxPloidy(2);

            newAttributes.replace(entry.getKey(), updateAttribute(entry.getKey(), entry.getValue(), vc.getAlleles(), sortedAlleles, type, ploidy));
        }

        // The above will have built new genotype for PL,AD, and SAC fields excluding the GT and GQ field, thus we must re-add them for comparison.
        for (int i = 0; i < newGT.size(); i++) {
            Genotype replacementGenotype = new GenotypeBuilder(newGT.get(i))
                    .GQ(vc.getGenotype(i).getGQ())
                    .phased(vc.getGenotype(i).isPhased())
                    .alleles(vc.getGenotype(i).getAlleles()).make();
            newGT.replace(replacementGenotype);
        }

        result.attributes(newAttributes).genotypes(newGT);

        return result.make();

    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private static Object updateAttribute(final String key, final Object value,
                                          final List<Allele> originalAlleles, final List<Allele> sortedAlleles,
                                          final VCFHeaderLineCount count, int ploidy) {
        if (key.startsWith("AS_")) {
            return remapASValues(value instanceof List? String.join(",", ((List<String>) value)) : (String) value, createAlleleIndexMap(originalAlleles, sortedAlleles));
        }else {
            switch (count) {
                case INTEGER:
                    return value;
                case UNBOUNDED:
                    //doesn't depend on allele ordering
                    return value;
                case A:
                    return remapATypeValues(attributeToList(value), createAlleleIndexMap(originalAlleles, sortedAlleles));
                case R:
                    return remapRTypeValues(attributeToList(value), createAlleleIndexMap(originalAlleles, sortedAlleles));
                case G:
                    return remapGTypeValues(attributeToList(value), originalAlleles, ploidy, sortedAlleles);
                default:
                    throw new GATKException("found unexpected vcf header count type: " + count);
            }
        }
    }

    static List<Integer> createAlleleIndexMap(final List<Allele> originalAlleles, final List<Allele> sortedAlleles){
        final List<Integer> mapping = new ArrayList<>(originalAlleles.size());
        for ( final Allele a: sortedAlleles){
            final int newIndex = originalAlleles.indexOf(a);
            mapping.add(newIndex);
        }
        return mapping;
    }

    static List<Object> remapRTypeValues(List<?> oldValue, List<Integer> mapping){
        return remapListValues(oldValue, mapping, 0);
    }

    private static List<Object> remapListValues(List<?> oldValue, List<Integer> mapping, int offset) {
        for( int i = 0; i < offset; i++){
            Utils.validate(mapping.get(i) == i, "values within the offset must not map outside the offset ");
        }
        final ArrayList<Object> reordered = new ArrayList<>(oldValue.size());
        for(int i = 0; i < oldValue.size(); i++){
            reordered.add(oldValue.get(mapping.get(i+offset) - offset));
        }
        return reordered;
    }

    static List<Object> remapATypeValues(List<?> oldValue, List<Integer> mapping){
        return remapListValues(oldValue, mapping, 1 );
    }

    static List<Object> remapGTypeValues(List<?> oldValue, List<Allele> originalAlleles, int ploidy, List<Allele> remappedAlleles){
        if (oldValue.size() == 1 && oldValue.get(0) instanceof String) {
            oldValue = Arrays.stream(((String) oldValue.get(0)).split(",")).collect(Collectors.toList());
        }

        List<Object> newValues = new ArrayList<>(oldValue.size());
        int[] subsettedGenotypes = AlleleSubsettingUtils.subsettedPLIndices(ploidy, originalAlleles, remappedAlleles);
        List<?> finalOldValue = oldValue;
        newValues.addAll(Arrays.stream(subsettedGenotypes).mapToObj(idx -> finalOldValue.get(idx)).collect(Collectors.toList()));

        return newValues;
    }

    static List<Object> remapASValues(String oldValue, List<Integer> mapping) {
        return remapListValues(Arrays.asList(oldValue.split("\\|")), mapping, 0);
    }

    //copied from htsjdk.variant.variantcontext.CommonInfo.getAttributeAsList for simplicity
    //maybe we should expose this as a static method in htsjdk?
    @SuppressWarnings("unchecked")
    private static List<Object> attributeToList(final Object attribute){
        if ( attribute == null ) return Collections.emptyList();
        if ( attribute instanceof List) return (List<Object>)attribute;
        if ( attribute.getClass().isArray() ) {
            if (attribute instanceof int[]) {
                return Arrays.stream((int[])attribute).boxed().collect(Collectors.toList());
            } else if (attribute instanceof double[]) {
                return Arrays.stream((double[])attribute).boxed().collect(Collectors.toList());
            }
            return Arrays.asList((Object[])attribute);
        }
        return Collections.singletonList(attribute);
    }

    public static void assertGenotypesAreEqual(final Genotype actual, final Genotype expected) {
        Assert.assertEquals(actual.getSampleName(), expected.getSampleName(), "Genotype names");
        Assert.assertTrue(CollectionUtils.isEqualCollection(actual.getAlleles(), expected.getAlleles()), "Genotype alleles");
        Assert.assertEquals(actual.getGenotypeString(false), expected.getGenotypeString(false), "Genotype string");
        Assert.assertEquals(actual.getType(), expected.getType(), "Genotype type");

        // filters are the same
        Assert.assertEquals(actual.getFilters(), expected.getFilters(), "Genotype fields");
        Assert.assertEquals(actual.isFiltered(), expected.isFiltered(), "Genotype isFiltered");

        // inline attributes
        Assert.assertEquals(actual.hasDP(), expected.hasDP(), "Genotype hasDP");
        Assert.assertEquals(actual.getDP(), expected.getDP(), "Genotype dp");
        Assert.assertEquals(actual.hasAD(), expected.hasAD(), "Genotype hasAD");
        Assert.assertEquals(actual.getAD(), expected.getAD(), "Genotype AD");
        Assert.assertEquals(actual.hasGQ(), expected.hasGQ(), "Genotype hasGQ");
        Assert.assertEquals(actual.getGQ(), expected.getGQ(), "Genotype gq");
        Assert.assertEquals(actual.hasPL(), expected.hasPL(), "Genotype hasPL");
        Assert.assertEquals(actual.getPL(), expected.getPL(), "Genotype PL");

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
                if (expectedValue instanceof List && actualValue instanceof List) {
                    // both values are lists, compare element b element
                    List<Object> expectedList = (List<Object>) expectedValue;
                    List<Object> actualList = (List<Object>) actualValue;
                    Assert.assertEquals(actualList.size(), expectedList.size());
                    for (int i = 0; i < expectedList.size(); i++) {
                        assertAttributeEquals(act.getKey(), actualList.get(i), expectedList.get(i));
                    }
                } else if (expectedValue instanceof List) {
                    // expected is a List but actual is not; normalize to String and compare
                    Assert.assertTrue(actualValue instanceof String, "Attempt to compare list to a non-string value");
                    final String expectedString = ((List<Object>) expectedValue).stream().map(v -> v.toString()).collect(Collectors.joining(","));
                    assertAttributeEquals(act.getKey(), actualValue, expectedString);
                }
                else if (actualValue instanceof List) {
                    // actual is a List but expected is not; normalize to String and compare
                    Assert.assertTrue(expectedValue instanceof String, "Attempt to compare list to a non-string value");
                    final String actualString = ((List<Object>) actualValue).stream().map(v -> v.toString()).collect(Collectors.joining(","));
                    assertAttributeEquals(act.getKey(), actualString, expectedValue);
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
        if ( value == null ) { return true; }
        else if ( value.equals(VCFConstants.MISSING_VALUE_v4) ) { return true; }
        else if ( value instanceof List ) {
            // handles the case where all elements are null or the list is empty
            for ( final Object elt : (List)value) {
                if (elt != null) {
                    return false;
                }
            }
            return true;
        } else {
            return false;
        }
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

        boolean passed = true;
        int numFailed = 0;

        for (int i = 0; i < v1.size(); i++) {
            if (! v1.get(i).toStringDecodeGenotypes().equals(v2.get(i).toStringDecodeGenotypes())){
                logger.error("Variant Comparison Error: " + "different element (compared by toStringDecodeGenotypes) " + i + ":\n" + v1.get(i) + "\n" + v2.get(i));
                passed = false;
                ++numFailed;
            }
        }
        if (!passed) {
            throw new AssertionError("Variant comparison failed!  Num non-matching variant pairs: " + numFailed);
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
        Assert.assertEquals(actual.getFilters(), expected.getFilters(), "filters");
        BaseTest.assertEqualsDoubleSmart(actual.getPhredScaledQual(), expected.getPhredScaledQual());

        assertVariantContextsHaveSameGenotypes(actual, expected);
    }

    private static Map<String, Object> filterIgnoredAttributes(final Map<String,Object> attributes, final List<String> attributesToIgnore){
        return attributes.entrySet().stream()
                .filter(p -> !attributesToIgnore.contains(p.getKey()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }

    // Method that determines whether two variant contexts have equivalent allele specific annotations regardless of allele ordering
    public static Boolean alleleSpecificAnnotationEquals(VariantContext actual, VariantContext expected, String annotation) {
        List<Allele> Aalleles = actual.getAlleles();
        String[] actualAnnotation = String.join(",",actual.getAttributeAsStringList(annotation, "")).split("\\|",-1);
        String[] expectedAnnotation = String.join(",",expected.getAttributeAsStringList(annotation, "")).split("\\|",-1);
        if (Arrays.equals(actualAnnotation, expectedAnnotation)) {
            return true;
        }if (actualAnnotation.length!=expectedAnnotation.length) {
            return false;
        }
        for (int i = 0; i < Aalleles.size(); i++) {
            Allele al = Aalleles.get(i);

            int k = expected.getAlleleIndex(al);

            if (!actualAnnotation[i].equals(expectedAnnotation[k])) {
                return false;
            }
        }
        return true;
    }

    public static void assertGenotypePosteriorsAttributeWasRemoved(final VariantContext actual, final VariantContext expected) {
        for (final Genotype g : actual.getGenotypes()) {
            Assert.assertFalse(g.hasExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
        }
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

    /**
     * Method which compares two variant contexts for equality regardless of different allele ordering.
     *
     * It functions by sorting the alleles in each variant context, and using the header to parse its attributes and
     * reorder them based on the new allele ordering.
     *
     * NOTES:
     * - For genotype fields, the order dependant fields PL, AD, and SAC are all recalculated, no guarantee
     *   is made about any other genotype fields which depend on the number of Alleles which might result in false negatives.
     * - This test requires that all attribute keys from the variant context are present, if one is writing a test and needs
     *   a complete header, consider {GATKVCFHeaderLine.getCompleteHeader()}
     *
     * @param actual                Variant context to test for equality
     * @param expected              Expected result
     * @param attributesToIgnore    Attributes we want to exclude from comparision
     * @param header                Header used to map behavior of annotations
     */
    public static void assertVariantContextsAreEqualAlleleOrderIndependent(final VariantContext actual, final VariantContext expected, final List<String> attributesToIgnore, VCFHeader header) {
        if (actual.getAlleles().equals(expected.getAlleles())) {
            assertVariantContextsAreEqual(actual, expected, attributesToIgnore);

        } else {
            VariantContext actualReordered = sortAlleles(actual, header);
            VariantContext expectedReordered = sortAlleles(expected, header);
            assertVariantContextsAreEqual(actualReordered, expectedReordered, attributesToIgnore);
        }
    }

    /**
     * Method which returns a complete header with all the GATK and HTSJDK standard header lines for testing purposes
     *
     * @return
     */
    public static VCFHeader getCompleteHeader() {
        Set<VCFHeaderLine> lines = new HashSet<>();

        // Adding HTSJDK lines
        VCFStandardHeaderLines.addStandardInfoLines(lines,false, Collections.emptyList());
        VCFStandardHeaderLines.addStandardFormatLines(lines,false, Collections.emptyList());

        lines.addAll(GATKVCFHeaderLines.getAllInfoLines());
        lines.addAll(GATKVCFHeaderLines.getAllFormatLines());
        lines.addAll(GATKVCFHeaderLines.getAllFilterLines());

        return new VCFHeader(lines);
    }

    /**
     * Uses the AnnotationsPlugin interface to instantiate and return a list of every annotation currently visible to the gatk
     */
    public static List<Annotation> getAllAnnotations() {
        CommandLineParser clp = new CommandLineArgumentParser(
                new Object(),
                Collections.singletonList(new GATKAnnotationPluginDescriptor(null, null)),
                Collections.emptySet());
        List<String> args = new ArrayList<>();
        args.add("--enable-all-annotations");
        clp.parseArguments(new PrintStream(new NullOutputStream()), args.toArray(new String[args.size()]));
        return instantiateAnnotations(clp);
    }

    private static List<Annotation> instantiateAnnotations(final CommandLineParser clp) {
        GATKAnnotationPluginDescriptor annotationPlugin = clp.getPluginDescriptor(GATKAnnotationPluginDescriptor.class);
        return annotationPlugin.getResolvedInstances();
    }

    public static Stream<VariantContext> streamVcf(final File vcf) {
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(vcf);
        return StreamSupport.stream(featureDataSource.spliterator(), false).onClose(() -> featureDataSource.close());
    }
}
