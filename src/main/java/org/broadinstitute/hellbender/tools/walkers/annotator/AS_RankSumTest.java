package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

/**
 * Allele-specific implementation of rank sum test annotations
 */
public abstract class AS_RankSumTest extends RankSumTest implements ReducibleAnnotation {
    private static final Logger logger = Logger.getLogger(AS_RankSumTest.class);
    protected static final String splitDelim = "\\|"; //String.split takes a regex, so we need to escape the pipe
    protected static final String printDelim = "|";
    protected static final String reducedDelim = ",";

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
                                         final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        return annotateRawData(ref, vc, stratifiedPerReadAlleleLikelihoodMap);
    }

    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap ) {
        if ( perReadAlleleLikelihoodMap == null) {
            return null;
        }

        final Map<String, Object> annotations = new HashMap<>();
        final AlleleSpecificAnnotationData<CompressedDataList<Integer>> myData = initializeNewAnnotationData(vc.getAlleles());
        calculateRawData(vc, perReadAlleleLikelihoodMap, myData);
        final String annotationString = makeRawAnnotationString(vc.getAlleles(), myData.getAttributeMap());
        if (annotationString == null){
            return null;
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
                annotationString.append(printDelim);
            }
            final CompressedDataList<Integer> alleleValues = perAlleleValues.get(vcAlleles.get(i));
            annotationString.append(alleleValues);
        }
        return annotationString.toString();
    }

    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    @Override
    public void calculateRawData(VariantContext vc, Map<String, PerReadAlleleLikelihoodMap> pralm, ReducibleAnnotationData myData) {
        if(pralm == null) {
            return;
        }

        final Map<Allele, CompressedDataList<Integer>> perAlleleValues = myData.getAttributeMap();
        for ( final PerReadAlleleLikelihoodMap likelihoodMap : pralm.values() ) {
            if ( likelihoodMap != null && !likelihoodMap.isEmpty() ) {
                fillQualsFromLikelihoodMap(vc.getAlleles(), vc.getStart(), likelihoodMap, perAlleleValues);
            }
        }
    }

    private void fillQualsFromLikelihoodMap(final List<Allele> alleles,
                                            final int refLoc,
                                            final PerReadAlleleLikelihoodMap likelihoodMap,
                                            final Map<Allele, CompressedDataList<Integer>> perAlleleValues) {
        for ( final Map.Entry<GATKRead, Map<Allele,Double>> el : likelihoodMap.getLikelihoodReadMap().entrySet() ) {
            final MostLikelyAllele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
            if ( ! a.isInformative() ) {
                continue; // read is non-informative
            }

            final GATKRead read = el.getKey();
            if ( isUsableRead(read, refLoc) ) {
                final OptionalDouble value = getElementForRead(read, refLoc, a);
                if (! value.isPresent() ) {
                    continue;
                }

                if(perAlleleValues.containsKey(a.getMostLikelyAllele())) {
                    perAlleleValues.get(a.getMostLikelyAllele()).add((int)value.getAsDouble());
                }
            }
        }
    }
}
