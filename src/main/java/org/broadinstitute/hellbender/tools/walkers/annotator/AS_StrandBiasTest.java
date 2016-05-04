package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.hellbender.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

/**
 * Allele-specific implementation of strand bias annotations
 */
public abstract class AS_StrandBiasTest extends StrandBiasTest implements ReducibleAnnotation {
    private final static Logger logger = Logger.getLogger(AS_StrandBiasTest.class);
    protected static final String splitDelim = "\\|"; //String.split takes a regex, so we need to escape the pipe
    protected static final String printDelim = "|";
    protected static final String reducedDelim = ",";
    protected static final int MIN_COUNT = 2;
    protected static final double MIN_PVALUE = 1E-320;
    protected static final int FORWARD = 0;
    protected static final int REVERSE = 1;
    protected final List<Integer> ZERO_LIST = new ArrayList<>();

    public AS_StrandBiasTest(){
        ZERO_LIST.add(0,0);
        ZERO_LIST.add(1,0);
    }

    @Override
    public List<VCFInfoHeaderLine> getDescriptions() {
        //TODO only raw for now
//        if (AnnotationUtils.walkerRequiresRawData(callingWalker))
            return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getRawKeyName()));
//        else
//            return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0)));
    }

    @Override
    public String getRawKeyName() { return GATKVCFConstants.AS_SB_TABLE_KEY; }

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

        //for allele-specific annotations we only call from HC and we only use perReadAlleleLikelihoodMap
        if ( perReadAlleleLikelihoodMap == null) {
            return new HashMap<>();
            // calculate the annotation from the stratified per read likelihood map
            // stratifiedPerReadAllelelikelihoodMap can come from HaplotypeCaller call to VariantAnnotatorEngine
        } else if (perReadAlleleLikelihoodMap != null) {
            final Map<String, Object> annotations = new HashMap<>();
            final ReducibleAnnotationData<List<Integer>> myData = new AlleleSpecificAnnotationData<>(vc.getAlleles(),null);
            calculateRawData(vc, perReadAlleleLikelihoodMap, myData);
            final Map<Allele, List<Integer>> perAlleleValues = myData.getAttributeMap();
            final String annotationString = makeRawAnnotationString(vc.getAlleles(), perAlleleValues);
            annotations.put(getRawKeyName(), annotationString);
            return annotations;
        } else {
            // for non-snp variants, we need per-read likelihoods.
            // for snps, we can get same result from simple pileup
            // for indels that do not have a computed strand bias (SB) or strand bias by sample (SBBS)
            return null;
        }
    }

    protected String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, List<Integer>> perAlleleValues) {
        String annotationString = "";
        for (final Allele a : vcAlleles) {
            if (!annotationString.isEmpty())
                annotationString += printDelim;
            List<Integer> alleleValues = perAlleleValues.get(a);
            if (alleleValues == null)
                alleleValues = ZERO_LIST;
            annotationString += encode(alleleValues);
        }
        return annotationString;
    }

    protected String encode(List<Integer> alleleValues) {
        String annotationString = "";
        for (int j =0; j < alleleValues.size(); j++) {
            annotationString += alleleValues.get(j);
            if (j < alleleValues.size()-1)
                annotationString += ",";
        }
        return annotationString;
    }

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME
    public void calculateRawData(final VariantContext vc,
                                 final Map<String, PerReadAlleleLikelihoodMap> pralm,
                                 final ReducibleAnnotationData rawAnnotations) {
        if(pralm == null) {
            return;
        }

        getStrandCountsFromLikelihoodMap(vc, pralm, rawAnnotations, MIN_COUNT);
    }

    /**
     Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     *             fw      rc
     *   allele1   #       #
     *   allele2   #       #
     * @return a 2x2 contingency table
     */
    public void getStrandCountsFromLikelihoodMap( final VariantContext vc,
                                                            final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                            final ReducibleAnnotationData<List<Integer>> perAlleleValues,
                                                            final int minCount) {
        if( stratifiedPerReadAlleleLikelihoodMap == null || vc == null ) {
            return;
        }

        final Allele ref = vc.getReference();
        final List<Allele> allAlts = vc.getAlternateAlleles();

        for (final PerReadAlleleLikelihoodMap maps : stratifiedPerReadAlleleLikelihoodMap.values() ) {
            final ReducibleAnnotationData<List<Integer>> sampleTable = new AlleleSpecificAnnotationData<>(vc.getAlleles(),null);
            for (final Map.Entry<GATKRead,Map<Allele,Double>> el : maps.getLikelihoodReadMap().entrySet()) {
                final MostLikelyAllele mostLikelyAllele = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue());
                final GATKRead read = el.getKey();
                updateTable(mostLikelyAllele.getAlleleIfInformative(), read, ref, allAlts, sampleTable);
            }
            //for each sample (value in stratified PRALM), only include it if there are >minCount informative reads
            if ( passesMinimumThreshold(sampleTable, minCount) ) {
                combineAttributeMap(sampleTable, perAlleleValues);
            }
        }
    }


    protected void combineAttributeMap(final ReducibleAnnotationData<List<Integer>> toAdd, final ReducibleAnnotationData<List<Integer>> combined) {
        for (final Allele a : combined.getAlleles()) {
            if (toAdd.hasAttribute(a) && toAdd.getAttribute(a) != null) {
                if (combined.getAttribute(a) != null) {
                    combined.getAttribute(a).set(0, (int) combined.getAttribute(a).get(0) + (int) toAdd.getAttribute(a).get(0));
                    combined.getAttribute(a).set(1, (int) combined.getAttribute(a).get(1) + (int) toAdd.getAttribute(a).get(1));
                }
                else {
                    List<Integer> alleleData = new ArrayList<>();
                    alleleData.add(0, toAdd.getAttribute(a).get(0));
                    alleleData.add(1, toAdd.getAttribute(a).get(1));
                    combined.putAttribute(a,alleleData);
                }
            }
        }
    }

    private void updateTable(final Allele bestAllele, final GATKRead read, final Allele ref, final List<Allele> allAlts, final ReducibleAnnotationData<List<Integer>> perAlleleValues) {

        final boolean matchesRef = bestAllele.equals(ref, true);
        final boolean matchesAnyAlt = allAlts.contains(bestAllele);

        //for uninformative reads
        if(bestAllele.isNoCall())
            return;

        //can happen if a read's most likely allele has been removed when --max_alternate_alleles is exceeded
        if (!( matchesRef || matchesAnyAlt ))
            return;

        final List<Integer> alleleStrandCounts;
        if (perAlleleValues.hasAttribute(bestAllele) && perAlleleValues.getAttribute(bestAllele) != null)
            alleleStrandCounts = perAlleleValues.getAttribute(bestAllele);
        else {
            alleleStrandCounts = new ArrayList<>();
            alleleStrandCounts.add(0,0);
            alleleStrandCounts.add(1,0);
        }
        final boolean isFW = !read.isReverseStrand();
        if (isFW) {
            alleleStrandCounts.set(FORWARD, alleleStrandCounts.get(FORWARD) + 1);
        } else {
            alleleStrandCounts.set(REVERSE, alleleStrandCounts.get(REVERSE) + 1);
        }
        perAlleleValues.putAttribute(bestAllele, alleleStrandCounts);
    }

    /**
     * Does this strand data array pass the minimum threshold for inclusion?
     *
     * @param sampleTable  the per-allele fwd/rev read counts for a single sample
     * @param minCount The minimum threshold of counts in the array
     * @return true if it passes the minimum threshold, false otherwise
     */
    protected boolean passesMinimumThreshold(final ReducibleAnnotationData<List<Integer>> sampleTable, final int minCount) {
        // the read total must be greater than MIN_COUNT
        int readTotal = 0;
        for (final List<Integer> alleleValues : sampleTable.getAttributeMap().values()) {
            if (alleleValues != null) {
                readTotal += alleleValues.get(FORWARD);
                readTotal += alleleValues.get(REVERSE);
            }
        }
        return readTotal > minCount;
    }

    @Override
    //This just calls the non-allele-specific code in StrandBiasTest.java
    protected abstract Map<String, Object> calculateAnnotationFromLikelihoodMap(final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                                                final VariantContext vc);

    @Override
    //Allele-specific annotations cannot be called from walkers other than HaplotypeCaller
    protected Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes){
        return new HashMap<>();
    }

    public static String rawValueAsString(int[][] table) {
        return table[0][1]+","+table[0][1]+printDelim+table[1][0]+","+table[1][1];
    }
}
