package org.broadinstitute.hellbender.tools.walkers.mutect.filtering;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import shaded.cloud_nio.com.google.errorprone.annotations.Var;

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
                return addFilter(filters, filterName);
            }
            else return filters;
        }).collect(Collectors.toList());
        return encodeASFilters(updatedFilters);
    }

    /**
     * Adds the new filter to the list of current filters. Takes care of replacing the SITE keyword
     * if there were no previous filters
     * @param currentFilters the current list of filter for the allele
     * @param newFilter the new filter to add
     * @return the new list of filters
     */
    protected static List<String> addFilter(List<String> currentFilters, String newFilter) {
        if (currentFilters.size() == 1 && currentFilters.contains(GATKVCFConstants.SITE_LEVEL_FILTERS)) {
            return Collections.singletonList(newFilter);
        } else {
            List<String> updated = new ArrayList<>();
            updated.addAll(currentFilters);
            updated.add(newFilter);
            return updated;
        }
    }

    /**
     * Sets the filters for each allele and calculates the intersection of the allele filters to set on the variant.
     * PASS if the intersection is empty.
     * @param vc The variant context to build from, however it assumes all relevant filters are set in the alleleFilters collection
     * @param alleleFilters filters to be applied to each allele, the intersection of these filters are applied at the site level
     * @return The updated variant context
     */
    public static VariantContext addAlleleAndComputeSiteFilters(VariantContext vc, List<Set<String>> alleleFilters) {
        String encodedFilters = AlleleFilterUtils.encodeASFilters(alleleFilters.stream().map(
                af -> af.isEmpty() ? Collections.singletonList(GATKVCFConstants.SITE_LEVEL_FILTERS) : af.stream().collect(Collectors.toList())).collect(Collectors.toList()));
        VariantContextBuilder vcb = new VariantContextBuilder(vc).attribute(GATKVCFConstants.AS_FILTER_STATUS_KEY, encodedFilters);

        Set<String> siteFilters = alleleFilters.stream().skip(1)
                .collect(()->new HashSet<>(alleleFilters.get(0)), Set::retainAll, Set::retainAll);

        if (!siteFilters.isEmpty()) {
            vcb.filters(siteFilters);
        } else {
            vcb.passFilters();
        }
        return vcb.make();
    }
}
