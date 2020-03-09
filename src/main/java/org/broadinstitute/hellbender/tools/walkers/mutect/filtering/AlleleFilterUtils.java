package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.ListIterator;
import java.util.stream.Collectors;

public class AlleleFilterUtils {

    public static List<List<String>> decodeASFilters(VariantContext vc) {
        return AnnotationUtils.decodeAnyASListWithRawDelim(vc.getCommonInfo().getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY, "")).stream()
                .map(filters -> AnnotationUtils.decodeAnyASList(filters).stream().map(String::trim).collect(Collectors.toList()))
                .collect(Collectors.toList());
    }

    public static String encodeASFilters(List<List<String>> filters) {
        return AnnotationUtils.encodeAnyASListWithRawDelim(filters.stream().map(alleleFilters -> AnnotationUtils.encodeStringList(alleleFilters)).collect(Collectors.toList()));
    }

    public static String getMergedASFilterString(VariantContext vc, List<Boolean> isFiltered, String filterName) {
        List<List<String>> alleleFilters = decodeASFilters(vc);
        Utils.validateArg(isFiltered.size() == alleleFilters.size(), "lists are not the same size");
        ListIterator<Boolean> isFilteredIt = isFiltered.listIterator();

        List<List<String>> updatedFilters = alleleFilters.stream().map(filters -> {
            Boolean filtered = isFilteredIt.next();
            if (filtered) {
                return addFilter(filters, filterName);
            }
            else return filters;
        }).collect(Collectors.toList());
        return encodeASFilters(updatedFilters);
    }

    public static List<String> addFilter(List<String> currentFilters, String newFilter) {
        if (currentFilters.size() == 1 && currentFilters.contains(GATKVCFConstants.SITE_LEVEL_FILTERS)) {
            return Collections.singletonList(newFilter);
        } else {
            List<String> updated = new ArrayList<>();
            updated.addAll(currentFilters);
            updated.add(newFilter);
            return updated;
        }
    }
}
