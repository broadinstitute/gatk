package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

/**
 * Allele-specific implementation of rank sum test annotations
 */
public abstract class AS_RankSumTest extends RankSumTest implements ReducibleAnnotation {
    private static final Logger logger = Logger.getLogger(AS_RankSumTest.class);
    public static final String SPLIT_DELIM = "\\|"; //String.split takes a regex, so we need to escape the pipe
    public static final String PRINT_DELIM = "|";
    public static final String REDUCED_DELIM = ",";

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        //TODO only raw for now
//        if (AnnotationUtils.walkerRequiresRawData(callingWalker))
            return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getRawKeyName()));
//        else
//            return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                         final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        return annotateRawData(ref, vc, likelihoods);
    }

    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final ReadLikelihoods<Allele> likelihoods ) {
        if ( likelihoods == null) {
            return Collections.emptyMap();
        }

        final Map<String, Object> annotations = new HashMap<>();
        final AlleleSpecificAnnotationData<CompressedDataList<Integer>> myData = initializeNewAnnotationData(vc.getAlleles());
        calculateRawData(vc, likelihoods, myData);
        final String annotationString = makeRawAnnotationString(vc.getAlleles(), myData.getAttributeMap());
        if (annotationString == null){
            return Collections.emptyMap();
        }
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    private AlleleSpecificAnnotationData<CompressedDataList<Integer>> initializeNewAnnotationData(final List<Allele> vcAlleles) {
        final Map<Allele, CompressedDataList<Integer>> perAlleleValues = new HashMap<>();
        for (final Allele a : vcAlleles) {
            perAlleleValues.put(a, new CompressedDataList<>());
        }
        final AlleleSpecificAnnotationData<CompressedDataList<Integer>> ret = new AlleleSpecificAnnotationData<>(vcAlleles, perAlleleValues.toString());
        ret.setAttributeMap(perAlleleValues);
        return ret;
    }

    private String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, CompressedDataList<Integer>> perAlleleValues) {
        if (perAlleleValues.values().stream().allMatch(d -> d.isEmpty())){
            return null;
        }
        final StringBuilder annotationString = new StringBuilder();
        for (int i =0; i< vcAlleles.size(); i++) {
            if (i!=0) {
                annotationString.append(PRINT_DELIM);
            }
            final CompressedDataList<Integer> alleleValues = perAlleleValues.get(vcAlleles.get(i));
            annotationString.append(alleleValues);
        }
        return annotationString.toString();
    }

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    @Override
    public void calculateRawData(VariantContext vc, final ReadLikelihoods<Allele> likelihoods, ReducibleAnnotationData myData) {
        if(likelihoods == null) {
            return;
        }

        final int refLoc = vc.getStart();

        final Map<Allele, CompressedDataList<Integer>> perAlleleValues = myData.getAttributeMap();
        for ( final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles() ) {
            if (bestAllele.isInformative() && isUsableRead(bestAllele.read, refLoc)) {
                final OptionalDouble value = getElementForRead(bestAllele.read, refLoc, bestAllele);
                if (value.isPresent() && value.getAsDouble() != INVALID_ELEMENT_FROM_READ && perAlleleValues.containsKey(bestAllele.allele)) {
                    perAlleleValues.get(bestAllele.allele).add((int) value.getAsDouble());
                }
            }
        }
    }
}
