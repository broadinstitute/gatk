package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.hellbender.cmdline.Argument;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Argument collection representing custom, dynamically discovered read filters.
 */
public class CustomReadFilterArgumentCollection implements Serializable {

    private static final long serialVersionUID = 1L;

    @Argument(fullName="customFilterName", optional=true)
    List<String> customReadFilters = new ArrayList<>();
}
