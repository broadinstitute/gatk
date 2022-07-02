package org.broadinstitute.hellbender.testutils;

import htsjdk.variant.variantcontext.Genotype;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * Represents a comparison between two {@link Genotype} objects.
 * Uses a builder pattern, so must be initialized using the contained {@link GenotypeComparison.Builder} class.
 * <br>
 * {@code GenotypeComparison comparison = new GenotypeComparison.Builder(actual, expected)
 *      .addAttributesToIgnore(attributesToIgnore)
 *      .addExtendedAttributesToIgnore(extendedAttributesToIgnore)
 *      .build();}
 * <br>
 * Building makes a call to the {@link #compare()} method.  Then results can then be retrieved as a
 * {@link GenotypeComparisonResults} object.
 * <br>
 * {@code GenotypeComparisonResults results = comparison.getResults();}
 */

public class GenotypeComparison {

    private Genotype actual;
    private Genotype expected;
    private Set<GenotypeAttributeEnum> attributesToIgnore;
    private List<String> extendedAttributesToIgnore;
    private boolean ignoreActualExtraAttributes;

    private GenotypeComparisonResults results;

    public static class Builder{

        private Genotype actual;
        private Genotype expected;
        private final Set<GenotypeAttributeEnum> attributesToIgnore = new HashSet<>();
        private final List<String> extendedAttributesToIgnore = new LinkedList<>();
        private boolean ignoreActualExtraAttributes = false;

        /**
         * Start a build of a {@link GenotypeComparison} object
         * @param actual The actual {@link Genotype} value to compare
         * @param expected The expected {@link Genotype} value
         */
        public Builder(final Genotype actual, final Genotype expected){
            if(actual == null || expected == null){
                throw new IllegalArgumentException("Genotypes being compared cannot be null");
            }
            this.actual = actual;
            this.expected = expected;
        }

        /**
         * Adds specified attributes to the list of attributes to ignore when comparing the {@link Genotype} objects
         * @param attributesToIgnore - Attributes to ignore
         * @return this
         */
        public Builder addAttributesToIgnore(final Set<GenotypeAttributeEnum> attributesToIgnore){
            this.attributesToIgnore.addAll(attributesToIgnore);

            return this;
        }

        /**
         * Adds specified attribute to the list of attributes to ignore when comparing the {@link Genotype} objects
         * @param attributeToIgnore - Attributes to ignore
         * @return this
         */
        public Builder addAttributeToIgnore(final GenotypeAttributeEnum attributeToIgnore){
            this.attributesToIgnore.add(attributeToIgnore);

            return this;
        }

        /**
         * Adds specified attributes to the list of attributes to ignore when comparing the attributes returned by the
         * {@link Genotype#getExtendedAttributes()} method of the {@link htsjdk.variant.variantcontext.Genotype} objects
         * @param extendedAttributesToIgnore - Attributes to ignore
         * @return this
         */
        public Builder addExtendedAttributesToIgnore(final List<String> extendedAttributesToIgnore){
            this.extendedAttributesToIgnore.addAll(extendedAttributesToIgnore);

            return this;
        }

        /**
         * Adds specified attribute to the list of attributes to ignore when comparing the attributes returned by the
         * {@link Genotype#getExtendedAttributes()} method of the {@link htsjdk.variant.variantcontext.Genotype} objects
         * @param extendedAttributeToIgnore - Attributes to ignore
         * @return this
         */
        public Builder addExtendedAttributeToIgnore(final String extendedAttributeToIgnore){
            this.extendedAttributesToIgnore.add(extendedAttributeToIgnore);

            return this;
        }

        /**
         * If called, the comparison will ignore any attributes returned by the actual {@link Genotype} object's
         * {@link Genotype#getExtendedAttributes()} method, which are not returned by the same method call for expected
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
         * Creates a {@link GenotypeComparison} object, filling in attributes based on what was specified
         * during the creation of this {@link GenotypeComparison.Builder} object
         *
         * @return a {@link GenotypeComparison} configured based on this {@link GenotypeComparison.Builder} object
         */
        public GenotypeComparison build(){
            final GenotypeComparison comparison = new GenotypeComparison();
            comparison.actual = this.actual;
            comparison.expected = this.expected;
            comparison.attributesToIgnore = this.attributesToIgnore;
            comparison.extendedAttributesToIgnore = this.extendedAttributesToIgnore;
            comparison.ignoreActualExtraAttributes = this.ignoreActualExtraAttributes;

            comparison.compare();

            return comparison;
        }
    }

    public Genotype getActual() {
        return actual;
    }

    public void setActual(final Genotype actual) {
        this.actual = actual;
    }

    public Genotype getExpected() {
        return expected;
    }

    public void setExpected(final Genotype expected) {
        this.expected = expected;
    }

    public Set<GenotypeAttributeEnum> getAttributesToIgnore() {
        return attributesToIgnore;
    }

    public void setAttributesToIgnore(final Set<GenotypeAttributeEnum> attributesToIgnore) {
        this.attributesToIgnore = attributesToIgnore;
    }

    public List<String> getExtendedAttributesToIgnore() {
        return extendedAttributesToIgnore;
    }

    public void setExtendedAttributesToIgnore(final List<String> extendedAttributesToIgnore) {
        this.extendedAttributesToIgnore = extendedAttributesToIgnore;
    }

    public GenotypeComparisonResults getResults(){
        return results;
    }

    private GenotypeComparison compare(){
        final GenotypeComparisonResults results = new GenotypeComparisonResults(actual, expected);

        for(final GenotypeAttributeEnum attribute : GenotypeAttributeEnum.values()){
            if(!attributesToIgnore.contains(attribute)){
                if(attribute == GenotypeAttributeEnum.EXTENDED_ATTRIBUTES){
                    final List<String> errorKeys = VariantContextTestUtils.checkAttributesEquals(
                            VariantContextTestUtils.filterIgnoredAttributes(actual.getExtendedAttributes(), extendedAttributesToIgnore),
                            VariantContextTestUtils.filterIgnoredAttributes(expected.getExtendedAttributes(), extendedAttributesToIgnore),
                            ignoreActualExtraAttributes
                    );
                    if(!errorKeys.isEmpty()){
                        results.addMismatchedAttribute(GenotypeAttributeEnum.EXTENDED_ATTRIBUTES);
                        results.addMismatchedExtendedAttributes(errorKeys);
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
