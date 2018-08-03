package org.broadinstitute.hellbender.tools.walkers.mutect;

import java.util.*;

public class FilterResult {
    private Set<String> filtersApplied = new HashSet<>();
    private Map<String, Object> newAnnotations = new HashMap<>();

    private double readOrientationPosterior;

    public void setReadOrientationPosterior(final double posterior){
        readOrientationPosterior = posterior;
    }

    public double getReadOrientationPosterior() {
        return readOrientationPosterior;
    }

    public void addFilter(final String filterName) {
        filtersApplied.add(filterName);
    }

    public void addAttribute(final String attributeName, final Object value){
        newAnnotations.put(attributeName, value);
    }

    public Set<String> getFilters(){
        return filtersApplied;
    }

    public Map<String, Object> getAttributes(){
        return newAnnotations;
    }
}

