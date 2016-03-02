package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Compendium of utils to work in GENOTYPE_GIVEN_ALLELES mode.
 */
public final class GenotypingGivenAllelesUtils {

    /**
     * Composes the given allele variant-context providing information about the rods and reference location.
     * @param tracker the meta data tracker.
     * @param loc the query location.
     * @param snpsOnly whether we only should consider SNP variation.
     * @param logger where to output warnings.
     * @param allelesBinding the target variation context binding containing the given alleles.
     * @return never {@code null}
     */
    public static VariantContext composeGivenAllelesVariantContextFromRod(final FeatureContext tracker,
                                                                          final Locatable loc,
                                                                          final boolean snpsOnly,
                                                                          final Logger logger,
                                                                          final FeatureInput<VariantContext> allelesBinding) {
        if ( tracker == null ) throw new IllegalArgumentException("the tracker cannot be null");
        if ( loc == null ) throw new IllegalArgumentException("the location cannot be null");
        if ( allelesBinding == null ) throw new IllegalArgumentException("the alleles binding cannot be null");

        VariantContext vc = null;

        // search for usable record
        for ( final VariantContext rodVc : tracker.getValues(allelesBinding, new SimpleInterval(loc)) ) {
            if ( rodVc != null && ! rodVc.isFiltered() && (! snpsOnly || rodVc.isSNP() )) {
                if ( vc == null ) {
                    vc = rodVc;
                } else {
                    if (logger != null) {
                        logger.warn("Multiple valid VCF records detected in the alleles input file at site "
                                + loc + ", only considering the first record");
                    }
                }
            }
        }

        return vc;
    }

    /**
     * Create the list of artificial GGA-mode haplotypes by injecting each of the provided alternate alleles into the reference haplotype
     *
     * @param refHaplotype the reference haplotype
     * @param givenHaplotypes the list of alternate alleles in VariantContexts
     * @param activeRegionWindow the window containing the reference haplotype
     *
     * @return a non-null list of haplotypes
     */
    public static List<Haplotype> composeGivenHaplotypes(final Haplotype refHaplotype, final List<VariantContext> givenHaplotypes, final GenomeLoc activeRegionWindow) {
        if (refHaplotype == null) throw new IllegalArgumentException("the reference haplotype cannot be null");
        if (givenHaplotypes == null) throw new IllegalArgumentException("given haplotypes cannot be null");
        if (activeRegionWindow == null) throw new IllegalArgumentException("active region window cannot be null");
        if (activeRegionWindow.size() != refHaplotype.length()) throw new IllegalArgumentException("inconsistent reference haplotype and active region window");

        final Set<Haplotype> returnHaplotypes = new LinkedHashSet<>();
        final int activeRegionStart = refHaplotype.getAlignmentStartHapwrtRef();

        for( final VariantContext compVC : givenHaplotypes ) {
            if (!GATKVariantContextUtils.overlapsRegion(compVC, activeRegionWindow)) {
                throw new IllegalArgumentException(" some variant provided does not overlap with active region window");
            }
            for( final Allele compAltAllele : compVC.getAlternateAlleles() ) {
                final Haplotype insertedRefHaplotype = refHaplotype.insertAllele(compVC.getReference(), compAltAllele, activeRegionStart + compVC.getStart() - activeRegionWindow.getStart(), compVC.getStart());
                if( insertedRefHaplotype != null ) { // can be null if the requested allele can't be inserted into the haplotype
                    returnHaplotypes.add(insertedRefHaplotype);
                }
            }
        }

        return new ArrayList<>(returnHaplotypes);
    }
}
