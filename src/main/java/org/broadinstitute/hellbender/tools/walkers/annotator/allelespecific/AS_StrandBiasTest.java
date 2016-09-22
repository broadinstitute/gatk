package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasTest;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.*;

/**
 * Allele-specific implementation of strand bias annotations
 */
public abstract class AS_StrandBiasTest extends StrandBiasTest implements ReducibleAnnotation {
    private final static Logger logger = Logger.getLogger(AS_StrandBiasTest.class);
    public static final String SPLIT_DELIM = "\\|"; //String.split takes a regex, so we need to escape the pipe
    public static final String PRINT_DELIM = "|";
    public static final String REDUCED_DELIM = ",";
    public static final int MIN_COUNT = 2;
    public static final double MIN_PVALUE = 1.0E-320;
    public static final int FORWARD = 0;
    public static final int REVERSE = 1;
    private final List<Integer> ZERO_LIST = new ArrayList<>();

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
                                        final ReadLikelihoods<Allele> likelihoods) {
        return annotateRawData(ref, vc, likelihoods);
    }

    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final ReadLikelihoods<Allele> likelihoods ) {

        //for allele-specific annotations we only call from HC and we only use likelihoods
        if ( likelihoods == null) {
            return Collections.emptyMap();
        }
        // calculate the annotation from the likelihoods
        // likelihoods can come from HaplotypeCaller call to VariantAnnotatorEngine
        final Map<String, Object> annotations = new HashMap<>();
        final ReducibleAnnotationData<List<Integer>> myData = new AlleleSpecificAnnotationData<>(vc.getAlleles(),null);
        calculateRawData(vc, likelihoods, myData);
        final Map<Allele, List<Integer>> perAlleleValues = myData.getAttributeMap();
        final String annotationString = makeRawAnnotationString(vc.getAlleles(), perAlleleValues);
        annotations.put(getRawKeyName(), annotationString);
        return annotations;
    }

    protected String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, List<Integer>> perAlleleValues) {
        String annotationString = "";
        for (final Allele a : vcAlleles) {
            if (!annotationString.isEmpty()) {
                annotationString += PRINT_DELIM;
            }
            List<Integer> alleleValues = perAlleleValues.get(a);
            if (alleleValues == null) {
                alleleValues = ZERO_LIST;
            }
            annotationString += encode(alleleValues);
        }
        return annotationString;
    }

    protected String encode(List<Integer> alleleValues) {
        String annotationString = "";
        for (int j =0; j < alleleValues.size(); j++) {
            annotationString += alleleValues.get(j);
            if (j < alleleValues.size()-1) {
                annotationString += ",";
            }
        }
        return annotationString;
    }

    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME
    public void calculateRawData(final VariantContext vc,
                                 final ReadLikelihoods<Allele> likelihoods,
                                 final ReducibleAnnotationData rawAnnotations) {
        if(likelihoods == null) {
            return;
        }

        getStrandCountsFromLikelihoodMap(vc, likelihoods, rawAnnotations, MIN_COUNT);
    }

    /**
     Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     *             fw      rc
     *   allele1   #       #
     *   allele2   #       #
     * @return a 2x2 contingency table
     */
    public void getStrandCountsFromLikelihoodMap( final VariantContext vc,
                                                  final ReadLikelihoods<Allele> likelihoods,
                                                  final ReducibleAnnotationData<List<Integer>> perAlleleValues,
                                                  final int minCount) {
        if( likelihoods == null || vc == null ) {
            return;
        }

        final Allele ref = vc.getReference();
        final List<Allele> allAlts = vc.getAlternateAlleles();

        for (final String sample : likelihoods.samples()) {
            final ReducibleAnnotationData<List<Integer>> sampleTable = new AlleleSpecificAnnotationData<>(vc.getAlleles(),null);
            likelihoods.bestAlleles(sample).stream()
                    .filter(ba -> ba.isInformative())
                    .forEach(ba -> updateTable(ba.allele, ba.read, ref, allAlts, sampleTable));
            if (passesMinimumThreshold(sampleTable, minCount)) {
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

        //can happen if a read's most likely allele has been removed when --max_alternate_alleles is exceeded
        if (!( matchesRef || matchesAnyAlt )) {
            return;
        }

        final List<Integer> alleleStrandCounts;
        if (perAlleleValues.hasAttribute(bestAllele) && perAlleleValues.getAttribute(bestAllele) != null) {
            alleleStrandCounts = perAlleleValues.getAttribute(bestAllele);
        } else {
            alleleStrandCounts = new ArrayList<>();
            alleleStrandCounts.add(0,0);
            alleleStrandCounts.add(1,0);
        }
        final boolean isForward = !read.isReverseStrand();
        if (isForward) {
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
        final int readCount = sampleTable.getAttributeMap().values().stream()
                .filter(alleleValues -> alleleValues != null)
                .mapToInt(alleleValues -> alleleValues.get(FORWARD) + alleleValues.get(REVERSE))
                .sum();
        return readCount > minCount;
    }

    @Override
    //Allele-specific annotations cannot be called from walkers other than HaplotypeCaller
    protected Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes){
        return Collections.emptyMap();
    }

    public static String rawValueAsString(int[][] table) {
        return table[0][0]+","+table[0][1]+ PRINT_DELIM +table[1][0]+","+table[1][1];
    }
}
