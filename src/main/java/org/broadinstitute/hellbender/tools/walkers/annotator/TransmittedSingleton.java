package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.samples.Trio;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;


/**
 * Existence of a transmitted or non-transmitted singleton in at least one of the given families
 *
 * <p>This annotation uses the genotype information from individuals in family trios to identify transmitted and non-transmitted singletons and the sample(s) in which they occur.
 * Transmitted singletons occur at sites in a cohort where the allele count is two and these two alleles occur in one parent and the child of a trio. A non-transmitted singleton
 * are sites with an allele count of one and this one allele occurs in a parent, but not the child of a trio. In both cases the other parent must have a high quality hom ref call.
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>The calculation assumes that the organism is diploid.</li>
 *     <li>This annotation requires a valid pedigree file.</li>
 *     <li>Non transmitted singletons are only valid for trios, not quads or other more complicated family structures.</li>
 *     <li>Only reports possible singletons for families where each of the three samples has high GQ (>20) and high depth (>20)</li>
 *     <li>Only reports possible singletons at sites with a Call Rate greater than 90% (meaning less than 10% of the samples at the given site were no-calls)</li>
 * </ul>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Existence of a transmitted (or non-transmitted) singleton in at least one of the given families (transmittedSingleton, nonTransmittedSingleton)")
public final class TransmittedSingleton extends PedigreeAnnotation implements InfoFieldAnnotation {
    protected final Logger warning = LogManager.getLogger(this.getClass());
    private Set<Trio> trios;
    private final int HI_GQ_THRESHOLD = 20;
    private final int HI_DP_THRESHOLD = 20;
    private final double CALL_RATE_THRESHOLD = 0.90;

    @VisibleForTesting
    public TransmittedSingleton(final Set<Trio> trios) {
        super((Set<String>) null);
        this.trios = Collections.unmodifiableSet(new LinkedHashSet<>(trios));
    }

    public TransmittedSingleton(final GATKPath pedigreeFile){
        super(pedigreeFile);
    }

    public TransmittedSingleton(){
        super((Set<String>) null);
    }

    private Set<Trio> initializeAndGetTrios() {
        if (trios == null) {
            trios = getTrios();
        }
        return trios;
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.TRANSMITTED_SINGLETON, GATKVCFConstants.NON_TRANSMITTED_SINGLETON);
    }
    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc);
        Set<Trio> trioSet = initializeAndGetTrios();
        if (!vc.isBiallelic() || trioSet.isEmpty()) {
            return Collections.emptyMap();
        }
        long highQualCalls = vc.getGenotypes().stream().filter(gt -> gt.getGQ() > HI_GQ_THRESHOLD).count();
        if ((double) highQualCalls / vc.getNSamples() < CALL_RATE_THRESHOLD) {
            return Collections.emptyMap();
        }
        final List<String> transmittedSingletonParent = new ArrayList<>();
        final List<String> nonTransmittedSingletonParent = new ArrayList<>();
        for (final Trio trio : trioSet) {
            if (vc.isBiallelic() &&
                    contextHasTrioGQs(vc, trio)) {

                final boolean oneParentHasAllele = (vc.getGenotype(trio.getMaternalID()).isHet() && vc.getGenotype(trio.getPaternalID()).isHomRef()) || (vc.getGenotype(trio.getMaternalID()).isHomRef() && vc.getGenotype(trio.getPaternalID()).isHet());
                final String matchingParentId = vc.getGenotype(trio.getMaternalID()).isHet() ? trio.getMaternalID() : trio.getPaternalID();

                final boolean momIsHighGQ = vc.getGenotype(trio.getMaternalID()).getGQ() >= HI_GQ_THRESHOLD;
                final boolean dadIsHighGQ = vc.getGenotype(trio.getPaternalID()).getGQ() >= HI_GQ_THRESHOLD;

                final boolean childIsHighGQHet = vc.getGenotype(trio.getChildID()).isHet() && vc.getGenotype(trio.getChildID()).getGQ() >= HI_GQ_THRESHOLD;
                final boolean childIsHighGQHomRef = vc.getGenotype(trio.getChildID()).isHomRef() && vc.getGenotype(trio.getChildID()).getGQ() >= HI_GQ_THRESHOLD;

                final boolean childIsHighDepth = vc.getGenotype(trio.getChildID()).getDP() >= HI_DP_THRESHOLD;
                final boolean momIsHighDepth = vc.getGenotype(trio.getChildID()).getDP() >= HI_DP_THRESHOLD;
                final boolean dadIsHighDepth = vc.getGenotype(trio.getChildID()).getDP() >= HI_DP_THRESHOLD;

                if (childIsHighDepth && momIsHighDepth && dadIsHighDepth &&
                        vc.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0) == 2 &&
                        childIsHighGQHet && oneParentHasAllele && momIsHighGQ && dadIsHighGQ) {
                        transmittedSingletonParent.add(matchingParentId);
                }
                //TODO: This only works for trios (not quads or other more complicated family structures that would effect number of singletons for parents or transmission to multiple kids)
                if (childIsHighDepth && momIsHighDepth && dadIsHighDepth &&
                        vc.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0) == 1 &&
                        childIsHighGQHomRef && momIsHighGQ && dadIsHighGQ) {
                        nonTransmittedSingletonParent.add(matchingParentId);
                }
            }
        }

        final Map<String, Object> attributeMap = new LinkedHashMap<>(1);
        if (!transmittedSingletonParent.isEmpty()) {
            attributeMap.put(GATKVCFConstants.TRANSMITTED_SINGLETON, transmittedSingletonParent);
        }
        if (!nonTransmittedSingletonParent.isEmpty()) {
            attributeMap.put(GATKVCFConstants.NON_TRANSMITTED_SINGLETON, nonTransmittedSingletonParent);
        }
        return attributeMap;
    }
}
