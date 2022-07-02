package org.broadinstitute.hellbender.testutils;

import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents the results of a comparison between two {@link VariantContext} objects.
 * Can be instantiated via the {@link VariantContextComparison} class.
 * <br>
 * {@code
 *      VariantContextComparisonResults results = new VariantContextComparison.Builder(actual, expected)
 *      .addVariantContextAttributesToIgnore(variantContextAttributesToIgnore)
 *      .addGenotypeExtendedAttributesToIgnore(genotypeExtendedAttributesToIgnore)
 *      .build()
 *      .getResults();
 * }
 */
public class VariantContextComparisonResults {

    private VariantContext actual;
    private VariantContext expected;
    private final Set<VariantContextAttributeEnum> mismatchedAttributes = new HashSet<>();
    private final Set<String> mismatchedExtendedAttributes = new HashSet<>();
    private final List<GenotypeComparisonResults> genotypeComparisonResults = new LinkedList<>();


    public VariantContextComparisonResults(final VariantContext actual, final VariantContext expected){
        this.actual = actual;
        this.expected = expected;
    }

    /**
     * @return true if there were no mismatches in the comparison, false otherwise
     */
    public boolean isMatch(){
        return mismatchedAttributes.isEmpty();
    }

    public void addMismatchedAttribute(final VariantContextAttributeEnum attributeName){
        mismatchedAttributes.add(attributeName);
    }

    public Set<VariantContextAttributeEnum> getMismatchedAttributes() {
        return mismatchedAttributes;
    }

    public void addMismatchedExtendedAttribute(final String attribute) {
        mismatchedExtendedAttributes.add(attribute);
    }

    public void addMismatchedExtendedAttributes(final Collection<String> attributes) {
        mismatchedExtendedAttributes.addAll(attributes);
    }

    public Set<String> getMismatchedExtendedAttributes(){
        return mismatchedExtendedAttributes;
    }

    public void addGenotypeComparisonResult(final GenotypeComparisonResults genotypeComparisonResult){
        this.genotypeComparisonResults.add(genotypeComparisonResult);
    }

    public List<GenotypeComparisonResults> getGenotypeComparisonResults() {
        return genotypeComparisonResults;
    }

    /**
     * Builds a String that gives a human-readable explanation of whether the compared objects match (according to the rules of the comparison)
     * along with a list of mismatched attributes, if there are any
     * @return a message detailing the results of the comparison, in brief
     */
    public String getResultStringConcise(){
        StringBuilder resultString = new StringBuilder();
        if(isMatch()){
            resultString.append("VariantContext actual(" + actual.getSource() + ") matches expected(" + expected.getSource() + ")");
        }
        else{
            resultString.append("Variant actual(" + actual.getSource() + ") does not match expected(" + expected.getSource() + ")");
            resultString.append(" with mismatches in the following attributes: " + mismatchedAttributes.stream().map(attr -> attr.getName()).collect(Collectors.joining(",")));
        }
        return resultString.toString();
    }

    /**
     * Builds a String that gives a human-readable explanation of whether the compared objects match (according to the rules of the comparison)
     * along with a list of mismatched attributes, including the values of those attributes
     * @return a message detailing the results of the comparison, in detail
     */
    public String getResultStringVerbose(){
        StringBuilder resultString = new StringBuilder();
        if(isMatch()){
            resultString.append("VariantContext actual(" + actual.getSource() + ") matches expected(" + expected.getSource() + ")");
        }
        else{
            resultString.append("VariantContext actual(" + actual.getSource() + ") does not match expected(" + expected.getSource() + ") with mismatches in the following attributes:\n");
            Object actualValue = null;
            Object expectedValue = null;
            for(final VariantContextAttributeEnum attribute : mismatchedAttributes){
                actualValue = attribute.getValue(actual);
                expectedValue = attribute.getValue(expected);
                if(attribute.equals(VariantContextAttributeEnum.ATTRIBUTES) && !this.mismatchedExtendedAttributes.isEmpty()){
                    resultString.append("  " + attribute.getName() + ": actual(" + actualValue + ") vs expected(" + expectedValue + ")\n");
                    for(final String extendedAttribute : this.mismatchedExtendedAttributes){
                        resultString.append("    " + extendedAttribute + ": actual(" + actual.getAttribute(extendedAttribute) + ") vs expected(" + expected.getAttribute(extendedAttribute) + ")\n");
                    }
                }
                else if(attribute.equals(VariantContextAttributeEnum.GENOTYPES)){
                    resultString.append("  Genotypes: actual(" + actual.getGenotypes().getSampleNamesOrderedByName().stream().collect(Collectors.joining(",")) + ") vs " +
                            "expected(" + expected.getGenotypes().getSampleNamesOrderedByName().stream().collect(Collectors.joining(",")) + ")\n");
                    for(final GenotypeComparisonResults genotypeResults : genotypeComparisonResults){
                        resultString.append("  " + genotypeResults.getResultStringVerbose());
                    }
                }
                else{
                    resultString.append("  " + attribute.getName() + ": actual(" + actualValue + ") vs expected(" + expectedValue + ")\n");
                }
            }
        }
        return resultString.toString();
    }
}
