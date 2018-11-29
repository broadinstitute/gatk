package org.broadinstitute.hellbender.engine.filters;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Set;
import java.util.LinkedHashSet;

/**
 * Keep only variants with any of these IDs.
 * Matching is done by case-sensitive exact match.
 */
public final class VariantIDsVariantFilter implements VariantFilter {
    private final static long serialVersionUID = 1L;

    private final Set<String> includeIDs = new LinkedHashSet<>();

    public VariantIDsVariantFilter(Set<String> keepIDs) {
        Utils.nonNull(keepIDs);
        includeIDs.addAll(keepIDs);
    }

    @Override
    public boolean test(final VariantContext vc) {
        if (vc.getID().indexOf(';') > 0) {
            String[] vc_ids = vc.getID().split(";");
            for (String vc_id : vc_ids) {
                if (includeIDs.contains(vc_id)) {
                    return true;
                }
            }
            return false;
        } else {
            return includeIDs.contains(vc.getID());
        }
    }
}
