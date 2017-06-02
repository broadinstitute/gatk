package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Contig class annotations (see {@link ContigGermlinePloidyAnnotation}).
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public enum ContigClass {
    /**
     * Annotation for autosomal (non-sex chromosome) contigs
     */
    AUTOSOMAL,

    /**
     * Annotation for allosomal (sex chromosome) contigs
     */
    ALLOSOMAL;

    public static final Set<String> CONTIG_CLASS_NAMES_SET = Collections.unmodifiableSet(new HashSet<>(
            Arrays.asList(ContigClass.values()).stream().map(ContigClass::name).collect(Collectors.toList())));
}
