package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.samples.MendelianViolation;
import org.broadinstitute.hellbender.utils.samples.Trio;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Existence of a de novo mutation in at least one of the given families
 *
 * <p>This annotation uses the genotype information from individuals in family trios to identify possible de novo mutations and the sample(s) in which they occur. This works best if the genotypes have been processed according to the <a href="https://www.broadinstitute.org/gatk/guide/article?id=4723">Genotype Refinement workflow</a>.</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>The calculation assumes that the organism is diploid.</li>
 *     <li>This annotation requires a valid pedigree file.</li>
 *     <li>Only reports possible de novos for children whose genotypes have not been tagged as filtered (which is most appropriate if parent likelihoods
 * have already been factored in using PhaseByTransmission).</li>
 *     <li>When multiple trios are present, the annotation is simply the maximum of the likelihood ratios, rather than the strict 1-Prod(1-p_i) calculation, as this can scale poorly for uncertain sites and many trios.</li>
 *     <li>This annotation can only be used from the Variant Annotator. If you attempt to use it from the UnifiedGenotyper, the run will fail with an error message to that effect. If you attempt to use it from the HaplotypeCaller, the run will complete successfully but the annotation will not be added to any variants.</li>
 * </ul>
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MVLikelihoodRatio.php">MVLikelihoodRatio</a></b> evaluates whether a site is transmitted from parents to offspring according to Mendelian rules or not.</li>
 * </ul>
 *
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Existence of a de novo mutation in at least one of the given families (hiConfDeNovo, loConfDeNovo)")
public final class PossibleDeNovo extends InfoFieldAnnotation {

    public PossibleDeNovo(final Set<Trio> trios, final double minGenotypeQualityP) {
        this.trios = Collections.unmodifiableSet(new LinkedHashSet<>(trios));
        mendelianViolation = new MendelianViolation(minGenotypeQualityP);
    }

    /**
     * This is dummy constructor that will do nothing until https://github.com/broadinstitute/gatk/issues/1880 is addressed
     */
    public PossibleDeNovo(){
        this(Collections.emptySet(), 0);
    }

    protected final OneShotLogger warning = new OneShotLogger(this.getClass());


    private static final int hi_GQ_threshold = 20; //WARNING - If you change this value, update the description in GATKVCFHeaderLines
    private static final int lo_GQ_threshold = 10; //WARNING - If you change this value, update the description in GATKVCFHeaderLines
    private static final double percentOfSamplesCutoff = 0.001; //for many, many samples use 0.1% of samples as allele frequency threshold for de novos
    private static final int flatNumberOfSamplesCutoff = 4;

    private final MendelianViolation mendelianViolation;
    private final Set<Trio> trios;

    @Override
    public List<String> getKeyNames() { return Arrays.asList(
            GATKVCFConstants.HI_CONF_DENOVO_KEY,
            GATKVCFConstants.LO_CONF_DENOVO_KEY); }

    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        if (trios.isEmpty()){
            warning.warn("Annotation will not be calculated, must provide a valid PED file (-ped) from the command line.");
            return Collections.emptyMap();
        }
        final List<String> highConfDeNovoChildren = new ArrayList<>();
        final List<String> lowConfDeNovoChildren = new ArrayList<>();
        for (final Trio trio : trios) {
            if (vc.isBiallelic() &&
                PossibleDeNovo.contextHasTrioLikelihoods(vc, trio) &&
                mendelianViolation.isViolation(trio.getMother(), trio.getFather(), trio.getChild(), vc) &&
                mendelianViolation.getParentsRefRefChildHet() > 0) {

                final int childGQ = vc.getGenotype(trio.getChildID()).getGQ();
                final int momGQ   = vc.getGenotype(trio.getMaternalID()).getGQ();
                final int dadGQ   = vc.getGenotype(trio.getPaternalID()).getGQ();

                if (childGQ >= hi_GQ_threshold && momGQ >= hi_GQ_threshold && dadGQ >= hi_GQ_threshold) {
                    highConfDeNovoChildren.add(trio.getChildID());
                } else if (childGQ >= lo_GQ_threshold && momGQ > 0 && dadGQ > 0) {
                    lowConfDeNovoChildren.add(trio.getChildID());
                }
            }
        }

        final double percentNumberOfSamplesCutoff = vc.getNSamples()*percentOfSamplesCutoff;
        final double AFcutoff = Math.max(flatNumberOfSamplesCutoff, percentNumberOfSamplesCutoff);
        final int deNovoAlleleCount = vc.getCalledChrCount(vc.getAlternateAllele(0)); //we assume we're biallelic above so use the first alt

        final Map<String,Object> attributeMap = new LinkedHashMap<>(2);
        if ( !highConfDeNovoChildren.isEmpty()  && deNovoAlleleCount < AFcutoff ) {
            attributeMap.put(GATKVCFConstants.HI_CONF_DENOVO_KEY, highConfDeNovoChildren);
        }
        if ( !lowConfDeNovoChildren.isEmpty()  && deNovoAlleleCount < AFcutoff ) {
            attributeMap.put(GATKVCFConstants.LO_CONF_DENOVO_KEY, lowConfDeNovoChildren);
        }
        return attributeMap;
    }

    private static boolean contextHasTrioLikelihoods(final VariantContext vc, final Trio trio) {
        final String mom = trio.getMaternalID();
        final String dad = trio.getPaternalID();
        final String kid = trio.getChildID();

        return   (!mom.isEmpty() && vc.hasGenotype(mom) && vc.getGenotype(mom).hasLikelihoods())
              && (!dad.isEmpty() && vc.hasGenotype(dad) && vc.getGenotype(dad).hasLikelihoods())
              && (!kid.isEmpty() && vc.hasGenotype(kid) && vc.getGenotype(kid).hasLikelihoods());
    }

}
