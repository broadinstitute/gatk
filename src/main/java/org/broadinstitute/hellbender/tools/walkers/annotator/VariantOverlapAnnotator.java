package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Annotate the ID field and attribute overlap FLAGs for a VariantContext against a FeatureContext or a list
 * of VariantContexts. It uses dbSNP to look up RSIDs for variants.  Note that while attempts are made to account
 * for variant representation issues, in rare cases differences in variant representation can lead to a variant
 * not being annotated when it optimally should be.
 */
public final class VariantOverlapAnnotator {
    private final FeatureInput<VariantContext> dbSNP;
    private final Map<FeatureInput<VariantContext>, String> overlaps;

    /**
     * Create a new VariantOverlapAnnotator
     */
    public VariantOverlapAnnotator(final FeatureInput<VariantContext> dbSNP) {
        this(dbSNP, Collections.emptyMap());
    }

    /**
     * Create a new VariantOverlapAnnotator
     *
     * @param dbSNP used for updating ID field values, or null if that behavior isn't desired
     * @param overlaps a map of FeatureInput -> name to use for overlap annotation.  Each feature input will be used to
     *                        add name => true for variants that overlap it each location.
     *                        Can be empty but not null
     */
    public VariantOverlapAnnotator(final FeatureInput<VariantContext> dbSNP, final Map<FeatureInput<VariantContext>, String> overlaps) {
        Utils.nonNull(overlaps, "overlaps cannot be null");

        this.dbSNP = dbSNP;
        this.overlaps = overlaps;
    }

    /**
     * Update rsID in vcToAnnotate with rsIDs from dbSNP fetched from featureContext
     * @see #annotateOverlap(java.util.List, String, htsjdk.variant.variantcontext.VariantContext)
     *
     * @param featureContext non-null featureContext, which we will use to update the rsID of vcToAnnotate
     * @param vcToAnnotate a variant context to annotate
     * @return a VariantContext (may be == to vcToAnnotate) with updated rsID value
     */
    public VariantContext annotateRsID(final FeatureContext featureContext, final VariantContext vcToAnnotate) {
        if ( dbSNP != null ) {
            final SimpleInterval loc = new SimpleInterval(vcToAnnotate);
            return annotateRsID(featureContext.getValues(dbSNP, loc.getStart()), vcToAnnotate);
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
     * Add overlap attributes to vcToAnnotate against all overlaps in featureContext
     *
     * @see #annotateOverlap(java.util.List, String, htsjdk.variant.variantcontext.VariantContext)
     * for more information
     *
     * @param featureContext non-null featureContext, which we will use to update the rsID of vcToAnnotate
     * @param vcToAnnotate a variant context to annotate
     * @return a VariantContext (may be == to vcToAnnotate) with updated overlaps update fields value
     */
    public VariantContext annotateOverlaps(final FeatureContext featureContext, final VariantContext vcToAnnotate) {
        if ( overlaps.isEmpty() ) {
            return vcToAnnotate;
        }

        VariantContext annotated = vcToAnnotate;
        final SimpleInterval loc = new SimpleInterval(vcToAnnotate);
        for ( final Map.Entry<FeatureInput<VariantContext>, String> overlap : overlaps.entrySet() ) {
            final FeatureInput<VariantContext> fi = overlap.getKey();
            final List<VariantContext> vcs = featureContext.getValues(fi, loc.getStart());
            annotated = annotateOverlap(vcs, overlap.getValue(), annotated);
        }

        return annotated;
    }

    /**
     * Add overlaps flag attributes to vcToAnnotate overlapTestVCs.getSource() => true if
     * an overlapping variant context can be found in overlapTestVCs with vcToAnnotate
     *
     * Overlaps here means that the reference alleles are the same and at least one alt
     * allele in vcToAnnotate is equals to one of the alt alleles in overlapTestVCs
     *
     * @param overlapTestVCs a non-null list of potential overlaps that start at vcToAnnotate
     * @param attributeKey the key to set to true in the attribute map for vcToAnnotate if it overlaps
     * @param vcToAnnotate a non-null VariantContext to annotate
     */
    public VariantContext annotateOverlap(final List<VariantContext> overlapTestVCs, final String attributeKey, final VariantContext vcToAnnotate) {
        if ( overlaps.isEmpty() ) {
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
     * Returns the ID fields (separated by semicolons) of all VariantContexts in rsIDSourceVCs that have the same reference and alternate alleles
     * as vcToAnnotate, after splitting to bi-allelics and trimming both the rsIDSourceVCs and vcToAnnotate
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
     * @return a String to use for the rsID from rsIDSourceVCs if any match, or null if none match
     */
    private static String getRsID(final List<VariantContext> rsIDSourceVCs, final VariantContext vcToAnnotate) {
        Utils.nonNull(rsIDSourceVCs, "rsIDSourceVCs cannot be null");
        Utils.nonNull(vcToAnnotate, "vcToAnnotate cannot be null");
        final List<String> rsids = new ArrayList<>();

        final List<VariantContext> vcAnnotateList = GATKVariantContextUtils.splitVariantContextToBiallelics(vcToAnnotate, true,
                GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS, true);

        for ( final VariantContext vcCompSource : rsIDSourceVCs ) {
            if ( vcCompSource.isFiltered() ) {
                continue; // don't process any failed VCs
            }

            if (!vcCompSource.getContig().equals(vcToAnnotate.getContig())) {
                throw new IllegalArgumentException("source rsID VariantContext " + vcCompSource + " is not on same chromosome as vcToAnnotate " + vcToAnnotate);
            }

            final List<VariantContext> vcCompList = GATKVariantContextUtils.splitVariantContextToBiallelics(vcCompSource, true,
                    GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS, true);
            boolean addThisID = false;
            for (final VariantContext vcComp : vcCompList) {
                for (final VariantContext vcToAnnotateBi : vcAnnotateList) {
                    if (vcComp.getStart() == vcToAnnotateBi.getStart() && vcToAnnotateBi.getReference().equals(vcComp.getReference()) && vcComp.getAlternateAlleles().equals(vcToAnnotateBi.getAlternateAlleles())) {
                        addThisID = true;
                        break;
                    }
                }

                if (addThisID) {
                    rsids.add(vcCompSource.getID());
                    break;
                }
            }
        }

        return rsids.isEmpty()? null : String.join(";", rsids);
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
     * Get the collection of names of overlap sets ued by this annotator.
     */
    public Collection<String> getOverlapNames() {
        return overlaps.values();
    }
}
