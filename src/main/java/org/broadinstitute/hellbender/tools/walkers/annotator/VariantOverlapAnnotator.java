package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Annotate the ID field and attribute overlap FLAGs for a VariantContext against a RefMetaDataTracker or a list
 * of VariantContexts
 */
public final class VariantOverlapAnnotator {
    private final FeatureInput<VariantContext> dbSNPBinding;
    private final Map<FeatureInput<VariantContext>, String> overlapBindings;

    /**
     * Create a new VariantOverlapAnnotator without overall bindings
     *
     * @see #VariantOverlapAnnotator(org.broadinstitute.gatk.utils.commandline.RodBinding, java.util.Map, org.broadinstitute.gatk.utils.GenomeLocParser)
     */
    public VariantOverlapAnnotator(final FeatureInput<VariantContext> dbSNPBinding) {
        this(dbSNPBinding, Collections.<FeatureInput<VariantContext>, String>emptyMap());
    }

    /**
     * Create a new VariantOverlapAnnotator
     *
     * @param dbSNPBinding the RodBinding to use for updating ID field values, or null if that behavior isn't desired
     * @param overlapBindings a map of RodBindings / name to use for overlap annotation.  Each binding will be used to
     *                        add name => true for variants that overlap with variants found to a
     *                        RefMetaDataTracker at each location.  Can be empty but not null
     */
    public VariantOverlapAnnotator(final FeatureInput<VariantContext> dbSNPBinding, final Map<FeatureInput<VariantContext>, String> overlapBindings) {
        Utils.nonNull(overlapBindings, "overlapBindings cannot be null");

        this.dbSNPBinding = dbSNPBinding;
        this.overlapBindings = overlapBindings;
    }

    /**
     * Update rsID in vcToAnnotate with rsIDs from dbSNPBinding fetched from tracker
     * @see #annotateOverlap(java.util.List, String, htsjdk.variant.variantcontext.VariantContext)
     *
     * @param tracker non-null tracker, which we will use to update the rsID of vcToAnnotate
     *                for VariantContexts bound to dbSNPBinding that start at vcToAnnotate
     * @param vcToAnnotate a variant context to annotate
     * @return a VariantContext (may be == to vcToAnnotate) with updated rsID value
     */
    public VariantContext annotateRsID(final FeatureContext tracker, final VariantContext vcToAnnotate) {
        if ( dbSNPBinding != null ) {
            final SimpleInterval loc = new SimpleInterval(vcToAnnotate);
            return annotateRsID(tracker.getValues(dbSNPBinding, loc.getStart()), vcToAnnotate);
        } else {
            return vcToAnnotate;
        }
    }

    /**
     * Update rsID of vcToAnnotate with rsID match found in vcsAtLoc, if one exists
     *
     * @param vcsAtLoc a list of variant contexts starting at this location to use as sources for rsID values
     * @param vcToAnnotate a variant context to annotate
     * @return a VariantContext (may be == to vcToAnnotate) with updated rsID value
     */
    public static VariantContext annotateRsID(final List<VariantContext> vcsAtLoc, final VariantContext vcToAnnotate) {
        final String rsID = getRsID(vcsAtLoc, vcToAnnotate);

        // add the ID if appropriate
        if ( rsID != null ) {
            final VariantContextBuilder vcb = new VariantContextBuilder(vcToAnnotate);

            if ( ! vcToAnnotate.hasID() ) {
                return vcb.id(rsID).make();
            } else if ( ! vcToAnnotate.getID().contains(rsID) ) {
                return vcb.id(vcToAnnotate.getID() + VCFConstants.ID_FIELD_SEPARATOR + rsID).make();
            } // falling through to return VC lower down
        }

        // nothing to do, just return vc
        return vcToAnnotate;
    }

    /**
     * Add overlap attributes to vcToAnnotate against all overlapBindings in tracker
     *
     * @see #annotateOverlap(java.util.List, String, htsjdk.variant.variantcontext.VariantContext)
     * for more information
     *
     * @param tracker non-null tracker, which we will use to update the rsID of vcToAnnotate
     *                for VariantContexts bound to dbSNPBinding that start at vcToAnnotate
     * @param vcToAnnotate a variant context to annotate
     * @return a VariantContext (may be == to vcToAnnotate) with updated overlaps update fields value
     */
    public VariantContext annotateOverlaps(final FeatureContext tracker, final VariantContext vcToAnnotate) {
        if ( overlapBindings.isEmpty() ) {
            return vcToAnnotate;
        }

        VariantContext annotated = vcToAnnotate;
        final SimpleInterval loc = new SimpleInterval(vcToAnnotate);
        for ( final Map.Entry<FeatureInput<VariantContext>, String> overlapBinding : overlapBindings.entrySet() ) {
            final FeatureInput<VariantContext> fi = overlapBinding.getKey();
            final List<VariantContext> vcs = tracker.getValues(fi, loc.getStart());
            annotated = annotateOverlap(vcs, overlapBinding.getValue(), annotated);
        }

        return annotated;
    }

    /**
     * Add overlaps flag attributes to vcToAnnotate binding overlapTestVCs.getSource() => true if
     * an overlapping variant context can be found in overlapTestVCs with vcToAnnotate
     *
     * Overlaps here means that the reference alleles are the same and at least one alt
     * allele in vcToAnnotate is equals to one of the alt alleles in overlapTestVCs
     *
     * @param overlapTestVCs a non-null list of potential overlaps that start at vcToAnnotate
     * @param attributeKey the key to set to true in the attribute map for vcToAnnotate if it overlaps
     * @param vcToAnnotate a non-null VariantContext to annotate
     * @return
     */
    public VariantContext annotateOverlap(final List<VariantContext> overlapTestVCs, final String attributeKey, final VariantContext vcToAnnotate) {
        if ( overlapBindings.isEmpty() ) {
            return vcToAnnotate;
        }

        final boolean overlaps = overlaps(overlapTestVCs, vcToAnnotate);
        if ( overlaps ) {
            return new VariantContextBuilder(vcToAnnotate).attribute(attributeKey, true).make();
        } else {
            return vcToAnnotate;
        }
    }

    /**
     * Returns the ID field of the first VariantContext in rsIDSourceVCs that has the same reference allele
     * as vcToAnnotate and all of the alternative alleles in vcToAnnotate.
     *
     * Doesn't require vcToAnnotate to be a complete match, so
     *
     * A/C/G in VC in rsIDSourceVCs
     *
     * would match the a VC with A/C but not A/T.  Also we don't require all alleles to match
     * so we would also match A/C/T to A/C/G.
     *
     * Will only match rsIDSourceVCs that aren't failing filters.
     *
     * @param rsIDSourceVCs a non-null list of potential overlaps that start at vcToAnnotate
     * @param vcToAnnotate a non-null VariantContext to annotate
     * @return a String to use for the rsID from rsIDSourceVCs if one matches, or null if none matches
     */
    private static String getRsID(final List<VariantContext> rsIDSourceVCs, final VariantContext vcToAnnotate) {
        Utils.nonNull(rsIDSourceVCs, "rsIDSourceVCs cannot be null");
        Utils.nonNull(vcToAnnotate, "vcToAnnotate cannot be null");

        for ( final VariantContext vcComp : rsIDSourceVCs ) {
            if ( vcComp.isFiltered() ) {
                continue; // don't process any failed VCs
            }

            if ( ! vcComp.getContig().equals(vcToAnnotate.getContig()) || vcComp.getStart() != vcToAnnotate.getStart() ) {
                throw new IllegalArgumentException("source rsID VariantContext " + vcComp + " doesn't start at the same position as vcToAnnotate " + vcToAnnotate);
            }

            if ( vcToAnnotate.getReference().equals(vcComp.getReference()) ) {
                for ( final Allele allele : vcToAnnotate.getAlternateAlleles() ) {
                    if ( vcComp.getAlternateAlleles().contains(allele) ) {
                        return vcComp.getID();
                    }
                }
            }
        }

        return null;
    }

    /**
     * Does vcToAnnotate overlap with any of the records in potentialOverlaps?
     *
     * @param potentialOverlaps a non-null list of potential overlaps that start at vcToAnnotate
     * @param vcToAnnotate a non-null VariantContext to annotate
     * @return true if vcToAnnotate overlaps (position and all alt alleles) with some variant in potentialOverlaps
     */
    private static boolean overlaps(final List<VariantContext> potentialOverlaps, final VariantContext vcToAnnotate) {
        return getRsID(potentialOverlaps, vcToAnnotate) != null;
    }

    /**
     * Get the collection of the RodBinding names for those being used for overlap detection
     * @return a non-null collection of Strings
     */
    public Collection<String> getOverlapNames() {
        return overlapBindings.values();
    }
}
