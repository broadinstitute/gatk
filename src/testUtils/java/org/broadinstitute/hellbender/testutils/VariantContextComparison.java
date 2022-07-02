package org.broadinstitute.hellbender.testutils;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Represents a comparison between two {@link VariantContext} objects.
 * Uses a builder pattern, so must be initialized using the contained {@link Builder} class.
 * <br>
 * {@code VariantContextComparison comparison = new VariantContextComparison.Builder(actual, expected)
 *      .addVariantContextAttributesToIgnore(variantContextAttributesToIgnore)
 *      .addGenotypeExtendedAttributesToIgnore(genotypeExtendedAttributesToIgnore)
 *      .build();}
 * <br>
 * Building makes a call to the {@link #compare()} method.  Then results can then be retrieved as a
 * {@link VariantContextComparisonResults} object.
 * <br>
 * {@code VariantContextComparisonResults results = comparison.getResults();}
 */
public class VariantContextComparison {

    private static final Set<VariantContextAttributeEnum> DEFAULT_ATTRIBUTES_TO_IGNORE = new HashSet<>(Arrays.asList(VariantContextAttributeEnum.SOURCE));

    private VariantContext actual;
    private VariantContext expected;
    private Set<VariantContextAttributeEnum> attributesToIgnore;
    private List<String> extendedAttributesToIgnore;
    private List<String> extendedAttributesWithJitter;
    private Set<GenotypeAttributeEnum> genotypeAttributesToIgnore;
    private List<String> genotypeExtendedAttributesToIgnore;
    private boolean ignoreActualExtraAttributes;

    private VariantContextComparisonResults results;

    /**
     * Class to use to create a {@link VariantContextComparison} object
     */
    public static class Builder{
        private VariantContext actual;
        private VariantContext expected;
        private final Set<VariantContextAttributeEnum> variantContextAttributesToIgnore = new HashSet<>();
        private final List<String> variantContextExtendedAttributesToIgnore = new LinkedList<>();
        private final List<String> variantContextExtendedAttributesWithJitter = new LinkedList<>();
        private final Set<GenotypeAttributeEnum> genotypeAttributesToIgnore = new HashSet<>();
        private final List<String> genotypeExtendedAttributesToIgnore = new LinkedList<>();
        private boolean ignoreActualExtraAttributes = false;

        /**
         * Start a build of a {@link VariantContextComparison} object
         * @param actual The actual {@link VariantContext} value to compare
         * @param expected The expected {@link VariantContext} value
         */
        public Builder(final VariantContext actual, final VariantContext expected){
            this.actual = actual;
            this.expected = expected;
        }

        /**
         * Sets the comparison to use default values for variantContextAttributesToIgnore (VariantContextAttributeEnum.SOURCE)
         * and sets ignoreActualExtraAttributes
         * @return this
         */
        public Builder useDefaultConfiguration(){
            this.variantContextAttributesToIgnore.addAll(DEFAULT_ATTRIBUTES_TO_IGNORE);
            this.ignoreActualExtraAttributes();

            return this;
        }

        /**
         * Adds specified attributes to the list of attributes to ignore when comparing the {@link VariantContext} objects
         * @param variantContextAttributesToIgnore - Attributes to ignore
         * @return this
         */
        public Builder addVariantContextAttributesToIgnore(final Set<VariantContextAttributeEnum> variantContextAttributesToIgnore){
            this.variantContextAttributesToIgnore.addAll(variantContextAttributesToIgnore);

            return this;
        }

        /**
         * Adds specified attribute to the list of attributes to ignore when comparing the {@link VariantContext} objects
         * @param variantContextAttributeToIgnore - Attribute to ignore
         * @return this
         */
        public Builder addVariantContextAttributeToIgnore(final VariantContextAttributeEnum variantContextAttributeToIgnore){
            this.variantContextAttributesToIgnore.add(variantContextAttributeToIgnore);

            return this;
        }

        /**
         * Adds specified attributes to the list of attributes returned by the {@link VariantContext} objects'
         * {@link VariantContext#getAttributes()} method to ignore during comparison
         * @param variantContextExtendedAttributesToIgnore - Attributes to ignore
         * @return this
         */
        public Builder addVariantContextExtendedAttributesToIgnore(final List<String> variantContextExtendedAttributesToIgnore){
            this.variantContextExtendedAttributesToIgnore.addAll(variantContextExtendedAttributesToIgnore);

            return this;
        }

        /**
         * Adds specified attribute to the list of attributes returned by the {@link VariantContext} objects'
         * {@link VariantContext#getAttributes()} method to ignore during comparison
         * @param variantContextExtendedAttributeToIgnore - Attribute to ignore
         * @return this
         */
        public Builder addVariantContextExtendedAttributeToIgnore(final String variantContextExtendedAttributeToIgnore){
            this.variantContextExtendedAttributesToIgnore.add(variantContextExtendedAttributeToIgnore);

            return this;
        }

        /**
         * Adds specified attributes to the list of attributes returned by the {@link VariantContext} objects'
         * {@link VariantContext#getAttributes()} method to only check presence of (without comparing values)
         * during comparison
         * @param variantContextExtendedAttributesWithJitter - Attributes to not compare values
         * @return this
         */
        public Builder addVariantContextExtendedAttributesWithJitter(final List<String> variantContextExtendedAttributesWithJitter){
            this.variantContextExtendedAttributesToIgnore.addAll(variantContextExtendedAttributesWithJitter);
            this.variantContextExtendedAttributesWithJitter.addAll(variantContextExtendedAttributesWithJitter);

            return this;
        }

        /**
         * Adds specified attribute to the list of attributes returned by the {@link VariantContext} objects'
         * {@link VariantContext#getAttributes()} method to only check presence of (without comparing values)
         * during comparison
         * @param variantContextExtendedAttributeWithJitter - Attribute to not compare values
         * @return this
         */
        public Builder addVariantContextExtendedAttributeWithJitter(final String variantContextExtendedAttributeWithJitter){
            this.variantContextExtendedAttributesToIgnore.add(variantContextExtendedAttributeWithJitter);
            this.variantContextExtendedAttributesWithJitter.add(variantContextExtendedAttributeWithJitter);

            return this;
        }

        /**
         * Adds specified attributes to the list of attributes to ignore when comparing the {@link htsjdk.variant.variantcontext.Genotype} objects
         * belonging to these {@link VariantContext} objects
         * @param genotypeAttributesToIgnore - Attributes to ignore
         * @return this
         */
        public Builder addGenotypeAttributesToIgnore(final Set<GenotypeAttributeEnum> genotypeAttributesToIgnore){
            this.genotypeAttributesToIgnore.addAll(genotypeAttributesToIgnore);

            return this;
        }

        /**
         * Adds specified attribute to the list of attributes to ignore when comparing the {@link htsjdk.variant.variantcontext.Genotype} objects
         * belonging to these {@link VariantContext} objects
         * @param genotypeAttributeToIgnore - Attribute to ignore
         * @return this
         */
        public Builder addGenotypeAttributeToIgnore(final GenotypeAttributeEnum genotypeAttributeToIgnore){
            this.genotypeAttributesToIgnore.add(genotypeAttributeToIgnore);

            return this;
        }

        /**
         * Adds specified attributes to the list of attributes to ignore when comparing the attributes returned by the
         * {@link Genotype#getExtendedAttributes()} method of the {@link htsjdk.variant.variantcontext.Genotype} objects
         * belonging to these {@link VariantContext} objects
         * @param genotypeExtendedAttributesToIgnore - Attributes to ignore
         * @return this
         */
        public Builder addGenotypeExtendedAttributesToIgnore(final List<String> genotypeExtendedAttributesToIgnore){
            this.genotypeExtendedAttributesToIgnore.addAll(genotypeExtendedAttributesToIgnore);

            return this;
        }

        /**
         * Adds specified attribute to the list of attributes to ignore when comparing the attributes returned by the
         * {@link Genotype#getExtendedAttributes()} method of the {@link htsjdk.variant.variantcontext.Genotype} objects
         * belonging to these {@link VariantContext} objects
         * @param genotypeExtendedAttributeToIgnore - Attribute to ignore
         * @return this
         */
        public Builder addGenotypeExtendedAttributeToIgnore(final String genotypeExtendedAttributeToIgnore){
            this.genotypeExtendedAttributesToIgnore.add(genotypeExtendedAttributeToIgnore);

            return this;
        }

        /**
         * If called, the comparison will ignore any attributes returned by the actual {@link VariantContext} object's
         * {@link VariantContext#getAttributes()} method, which are not returned by the same method call for expected
         * <br>
         * i.e. If actual has attributes that expected does not have, that will not be considered a mismatch
         *
         * @return this
         */
        public Builder ignoreActualExtraAttributes(){
            this.ignoreActualExtraAttributes = true;

            return this;
        }

        /**
         * Creates a {@link VariantContextComparison} object, filling in attributes based on what was specified
         * during the creation of this {@link Builder} object
         *
         * @return a {@link VariantContextComparison} configured based on this {@link Builder} object
         */
        public VariantContextComparison build(){
            final VariantContextComparison comparison = new VariantContextComparison();
            comparison.actual = this.actual;
            comparison.expected = this.expected;
            comparison.attributesToIgnore = this.variantContextAttributesToIgnore;
            comparison.extendedAttributesToIgnore = this.variantContextExtendedAttributesToIgnore;
            comparison.extendedAttributesWithJitter = this.variantContextExtendedAttributesWithJitter;
            comparison.genotypeAttributesToIgnore = this.genotypeAttributesToIgnore;
            comparison.genotypeExtendedAttributesToIgnore = this.genotypeExtendedAttributesToIgnore;
            comparison.ignoreActualExtraAttributes = this.ignoreActualExtraAttributes;

            comparison.compare();

            return comparison;
        }
    }

    private VariantContextComparison(){};

    public VariantContext getActual() {
        return actual;
    }

    public void setActual(final VariantContext actual) {
        this.actual = actual;
    }

    public VariantContext getExpected() {
        return expected;
    }

    public void setExpected(final VariantContext expected) {
        this.expected = expected;
    }

    public Set<VariantContextAttributeEnum> getAttributesToIgnore() {
        return attributesToIgnore;
    }

    public void setAttributesToIgnore(final Set<VariantContextAttributeEnum> attributesToIgnore) {
        this.attributesToIgnore = attributesToIgnore;
    }

    public List<String> getExtendedAttributesToIgnore() {
        return extendedAttributesToIgnore;
    }

    public void setExtendedAttributesToIgnore(final List<String> extendedAttributesToIgnore) {
        this.extendedAttributesToIgnore = extendedAttributesToIgnore;
    }

    public List<String> getExtendedAttributesWithJitter() {
        return extendedAttributesWithJitter;
    }

    public void setExtendedAttributesWithJitter(final List<String> extendedAttributesWithJitter) {
        this.extendedAttributesWithJitter = extendedAttributesWithJitter;
    }

    public Set<GenotypeAttributeEnum> getGenotypeAttributesToIgnore() {
        return genotypeAttributesToIgnore;
    }

    public void setGenotypeAttributesToIgnore(final Set<GenotypeAttributeEnum> genotypeAttributesToIgnore) {
        this.genotypeAttributesToIgnore = genotypeAttributesToIgnore;
    }

    public List<String> getGenotypeExtendedAttributesToIgnore() {
        return genotypeExtendedAttributesToIgnore;
    }

    public void setGenotypeExtendedAttributesToIgnore(final List<String> genotypeExtendedAttributesToIgnore) {
        this.genotypeExtendedAttributesToIgnore = genotypeExtendedAttributesToIgnore;
    }

    public VariantContextComparisonResults getResults(){
        return results;
    }

    private VariantContextComparison compare(){
        final VariantContextComparisonResults results = new VariantContextComparisonResults(this.actual, this.expected);

        if(actual == null || expected == null){
            throw new IllegalArgumentException("Variant contexts to be compared must not be null");
        }

        for(final VariantContextAttributeEnum attribute : VariantContextAttributeEnum.values()){
            if(!attributesToIgnore.contains(attribute)){
                if(attribute == VariantContextAttributeEnum.ATTRIBUTES){
                    //Compare attributes, ignoring attributes in extendAttributesToIgnore
                    final List<String> errorKeys = VariantContextTestUtils.checkAttributesEquals(
                            VariantContextTestUtils.filterIgnoredAttributes(actual.getAttributes(), extendedAttributesToIgnore),
                            VariantContextTestUtils.filterIgnoredAttributes(expected.getAttributes(), extendedAttributesToIgnore),
                            ignoreActualExtraAttributes
                    );
                    //Verify attributes with jitter exist in both
                    if(!extendedAttributesWithJitter.isEmpty()){
                        errorKeys.addAll(VariantContextTestUtils.listMissingAttributes(expected.getAttributes(), actual.getAttributes(), extendedAttributesWithJitter));
                    }
                    if(!errorKeys.isEmpty()){
                        results.addMismatchedAttribute(VariantContextAttributeEnum.ATTRIBUTES);
                        results.addMismatchedExtendedAttributes(errorKeys);
                    }
                }
                else if(attribute == VariantContextAttributeEnum.GENOTYPES){
                    if((!actual.hasGenotypes() && expected.hasGenotypes()) || (actual.hasGenotypes() && !expected.hasGenotypes())){
                        results.addMismatchedAttribute(VariantContextAttributeEnum.GENOTYPES);
                    }
                    else if(actual.hasGenotypes() && expected.hasGenotypes()){
                        for(String sampleName : expected.getSampleNames()){
                            if(actual.hasGenotype(sampleName)){
                                final GenotypeComparisonResults comparisonResults = new GenotypeComparison.Builder(actual.getGenotype(sampleName), expected.getGenotype(sampleName))
                                        .addAttributesToIgnore(this.genotypeAttributesToIgnore)
                                        .addExtendedAttributesToIgnore(this.genotypeExtendedAttributesToIgnore)
                                        .build()
                                        .getResults();
                                results.addGenotypeComparisonResult(comparisonResults);
                                if(!results.getMismatchedAttributes().contains(VariantContextAttributeEnum.GENOTYPES) && !comparisonResults.isMatch()){
                                    results.addMismatchedAttribute(VariantContextAttributeEnum.GENOTYPES);
                                }
                            }
                            else{
                                if(!results.getMismatchedAttributes().contains(VariantContextAttributeEnum.GENOTYPES)){
                                    results.addMismatchedAttribute(VariantContextAttributeEnum.GENOTYPES);
                                }
                            }
                        }
                    }
                }
                else if(attribute == VariantContextAttributeEnum.PHRED_SCALED_QUALITY){
                    if(!BaseTest.equalsDoubleSmart(actual.getPhredScaledQual(), expected.getPhredScaledQual())){
                        results.addMismatchedAttribute(VariantContextAttributeEnum.PHRED_SCALED_QUALITY);
                    }
                }
                else{
                    if(!VariantContextTestUtils.checkFieldEquals(attribute.getValue(actual), attribute.getValue(expected))){
                        results.addMismatchedAttribute(attribute);
                    }
                }
            }
        }

        this.results = results;
        return this;
    }

}
