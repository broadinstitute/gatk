package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Helps read and set allele specific filters in the INFO field. Also, helps with updating the Filter field
 */
public class AlleleFilterUtils {

    /**
     * Decode the AS_FilterStatus INFO attribute. It is important to trim the strings since the
     * parser puts spaces in the coded string. SITE can be returned in the list (i.e. it is not
     * removed during the processing).
     * @param vc the variant context to read the AS_FilterStatus attribute from
     * @return A list for each alt allele which contains a list of the filters that apply
     */
    public static List<List<String>> decodeASFilters(VariantContext vc) {
        return AnnotationUtils.decodeAnyASListWithRawDelim(vc.getCommonInfo().getAttributeAsString(GATKVCFConstants.AS_FILTER_STATUS_KEY, "")).stream()
                .map(filters -> AnnotationUtils.decodeAnyASList(filters).stream().map(String::trim).collect(Collectors.toList()))
                .collect(Collectors.toList());
    }

    /**
     * Create the encoded string for AS_FilterStatus from the list of filters for each alt allele.
     * The method assumes that SITE has been inserted for any empty filter lists
     * @param filters a list for each alt allele that contains a list of the filters that apply to it
     * @return the encoded string
     */
    public static String encodeASFilters(List<List<String>> filters) {
        return AnnotationUtils.encodeAnyASListWithRawDelim(filters.stream().map(alleleFilters -> AnnotationUtils.encodeStringList(alleleFilters)).collect(Collectors.toList()));
    }

    /**
     * Takes a list of boolean values that indicate whether the filter applies to each of the alternate alleles
     * @param vc variant context to use to get existing allele filters
     * @param isFiltered list for alternate alleles of whether the specified filter should apply
     * @param filterName the name of the filter to apply
     * @return the new encoded string for the AS_FilterStatus INFO attribute
     */
    public static String getMergedASFilterString(VariantContext vc, List<Boolean> isFiltered, String filterName) {
        List<List<String>> alleleFilters = decodeASFilters(vc);
        Utils.validateArg(isFiltered.size() == alleleFilters.size(), "lists are not the same size");
        ListIterator<Boolean> isFilteredIt = isFiltered.listIterator();

        List<List<String>> updatedFilters = alleleFilters.stream().map(filters -> {
            Boolean filtered = isFilteredIt.next();
            if (filtered) {
                return addAlleleFilters(filters, Collections.singletonList(filterName));
            }
            else return filters;
        }).collect(Collectors.toList());
        return encodeASFilters(updatedFilters);
    }

    /**
     * Adds the new filters to the list of current filters. Takes care of replacing the SITE keyword
     * if there were no previous filters
     * @param currentAlleleFilters the current list of filter for the allele
     * @param newFilters the new filters to add
     * @return the updated list of filters
     */
    protected static List<String> addAlleleFilters(List<String> currentAlleleFilters, List<String> newFilters) {
        if (newFilters.isEmpty()) {
            return currentAlleleFilters;
        } else if (currentAlleleFilters.isEmpty() || (currentAlleleFilters.size() == 1 && currentAlleleFilters.contains(GATKVCFConstants.SITE_LEVEL_FILTERS))) {
            // new filters is not empty and there are no filters currently set for the allele
            return newFilters;
        } else {
            LinkedHashSet<String> updated = new LinkedHashSet<>();
            updated.addAll(currentAlleleFilters);
            updated.addAll(newFilters);
            return updated.stream().collect(Collectors.toList());
        }

    }

    /**
     * Adds the new allele filters to the existing allele filters in the vc. Computes whether there are
     * new site filters and updates the filter in the vc. If there are no site filters, sets filters to pass
     * Sets the filters for each allele and calculates the intersection of the allele filters to set on the variant.
     * PASS if the intersection is empty.
     * @param vc The variant context to add the filters to, both at the allele and site level
     * @param newAlleleFilters filters to be applied to each allele, the intersection of these filters are applied at the site level
     * @param invalidatePreviousFilters whether existing filters should be removed
     * @return The updated variant context
     */
    public static VariantContext addAlleleAndSiteFilters(VariantContext vc, List<Set<String>> newAlleleFilters, boolean invalidatePreviousFilters) {
        if (newAlleleFilters.isEmpty()) {
            return vc;
        }
        List<List<String>> currentAlleleFilters = decodeASFilters(vc);
        if (!currentAlleleFilters.isEmpty() && newAlleleFilters.size() != currentAlleleFilters.size()) {
            // log an error
            return vc;
        }

        if (currentAlleleFilters.isEmpty() || invalidatePreviousFilters) {
            currentAlleleFilters = new ArrayList<>(Collections.nCopies(newAlleleFilters.size(), Collections.singletonList(GATKVCFConstants.SITE_LEVEL_FILTERS)));
        }
        ListIterator<List<String>> currentAlleleFiltersIt = currentAlleleFilters.listIterator();
        List<List<String>> updatedAlleleFilters = newAlleleFilters.stream().map(newfilters -> addAlleleFilters(currentAlleleFiltersIt.next(), newfilters.stream().collect(Collectors.toList()))).collect(Collectors.toList());
        String encodedFilters = encodeASFilters(updatedAlleleFilters);
        VariantContextBuilder vcb = new VariantContextBuilder(vc).attribute(GATKVCFConstants.AS_FILTER_STATUS_KEY, encodedFilters);

        if (invalidatePreviousFilters) {
            vcb.unfiltered();
        }
        Set<String> siteFiltersToAdd = newAlleleFilters.stream().skip(1)
                .collect(()->new HashSet<>(newAlleleFilters.get(0)), Set::retainAll, Set::retainAll);
        siteFiltersToAdd.stream().forEach(filter -> vcb.filter(filter));
        if ((vcb.getFilters() == null || vcb.getFilters().isEmpty()) && !invalidatePreviousFilters) {
            vcb.passFilters();
        }
        return vcb.make();
    }
}
