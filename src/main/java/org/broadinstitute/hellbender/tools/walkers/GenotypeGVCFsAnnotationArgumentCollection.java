package org.broadinstitute.hellbender.tools.walkers;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.DefaultGATKVariantAnnotationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class GenotypeGVCFsAnnotationArgumentCollection extends DefaultGATKVariantAnnotationArgumentCollection {
    private static final long serialVersionUID = 1L;
    
    public static final String KEEP_SPECIFIED_RAW_COMBINED_ANNOTATION_LONG_NAME = "keep-specific-combined-raw-annotation";
    public static final String KEEP_SPECIFIED_RAW_COMBINED_ANNOTATION_SHORT_NAME = "keep-specific-combined";

    /**
     * Keep only the specific combined raw annotations specified (removing the other raw annotations if keep-combined-raw-annotations is not set). See {@link ReducibleAnnotation}
     */
    @Argument(fullName= KEEP_SPECIFIED_RAW_COMBINED_ANNOTATION_LONG_NAME, shortName = KEEP_SPECIFIED_RAW_COMBINED_ANNOTATION_SHORT_NAME, optional = true,
            mutex = {GenotypeGVCFs.KEEP_COMBINED_LONG_NAME},
            doc="Keep only the specific combined raw annotations specified (removing the other raw annotations). Duplicate values will be ignored.")
    protected List<String> keepSpecifiedCombined = new ArrayList<>();

    @Override
    public List<String> getKeepSpecifiedCombinedAnnotationNames() {return Collections.unmodifiableList(keepSpecifiedCombined);}
}
