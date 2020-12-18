package org.broadinstitute.hellbender.testutils;

import htsjdk.variant.variantcontext.Genotype;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Represents the results of a comparison between two {@link Genotype} objects.
 * Can be instantiated via the {@link GenotypeComparison} class.
 * <br>
 * {@code
 *      GenotypeComparisonResults results = new GenotypeComparison.Builder(actual, expected)
 *      .addAttributesToIgnore(attributesToIgnore)
 *      .addExtendedAttributesToIgnore(extendedAttributesToIgnore)
 *      .build()
 *      .getResults();
 * }
 */
public class GenotypeComparisonResults {

    private Genotype actual;
    private Genotype expected;
    private final Set<GenotypeAttributeEnum> mismatchedAttributes = new HashSet<>();
    private final Set<String> mismatchedExtendedAttributes = new HashSet<>();

    public GenotypeComparisonResults(final Genotype actual, final Genotype expected){
        this.actual = actual;
        this.expected = expected;
    }

    /**
     * @return true if there were no mismatches in the comparison, false otherwise
     */
    public boolean isMatch(){
        return mismatchedAttributes.isEmpty();
    }

    public void addMismatchedAttribute(final GenotypeAttributeEnum attributeName){
        mismatchedAttributes.add(attributeName);
    }

    public Set<GenotypeAttributeEnum> getMismatchedAttributes(){
        return mismatchedAttributes;
    }

    public void addMismatchedExtendedAttribute(final String attribute){
        mismatchedExtendedAttributes.add(attribute);
    }

    public void addMismatchedExtendedAttributes(final Collection<String> attributes){
        mismatchedExtendedAttributes.addAll(attributes);
    }

    public Set<String> getMismatchedExtendedAttributes(){
        return mismatchedExtendedAttributes;
    }

    /**
     * Builds a String that gives a human-readable explanation of whether the compared objects match (according to the rules of the comparison)
     * along with a list of mismatched attributes, if there are any
     * @return a message detailing the results of the comparison, in brief
     */
    public String getResultStringConcise(){
        StringBuilder resultString = new StringBuilder();
        if(isMatch()){
            resultString.append("Genotype actual(" + actual.getSampleName() + ") matches expected(" + expected.getSampleName() + ")");
        }
        else{
            resultString.append("Genotype actual(" + actual.getSampleName() + ") does not match expected(" + expected.getSampleName() + ")");
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
            resultString.append("Genotype actual(" + actual.getSampleName() + ") matches expected(" + expected.getSampleName() + ")");
        }
        else{
            resultString.append("Genotype actual(" + actual.getSampleName() + ") does not match expected(" + expected.getSampleName() + ") with mismatches in the following attributes:\n");
            Object actualValue = null;
            Object expectedValue = null;
            for(final GenotypeAttributeEnum attribute : mismatchedAttributes){
                actualValue = attribute.getValue(actual);
                expectedValue = attribute.getValue(expected);
                resultString.append("  " + attribute.getName() + ": actual(" + actualValue + ") vs expected(" + expectedValue + ")\n");
                if(attribute.equals(GenotypeAttributeEnum.EXTENDED_ATTRIBUTES) && !this.mismatchedExtendedAttributes.isEmpty()){
                    for(final String extendedAttribute : this.mismatchedExtendedAttributes){
                        resultString.append("    " + extendedAttribute + ": actual(" + actual.getExtendedAttribute(extendedAttribute) + ") vs expected(" + expected.getExtendedAttribute(extendedAttribute) + ")\n");
                    }
                }
            }
        }
        return resultString.toString();
    }
}
