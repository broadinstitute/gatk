package org.broadinstitute.hellbender.cmdline.GATKPlugin;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.List;

/**
 * Arguments for requesting VariantContext annotations to be processed by {@link org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine}
 * for tools that process variants objects.
 */
public class DefaultGATKVariantAnnotationArgumentCollection extends GATKAnnotationArgumentCollection {
    private static final long serialVersionUID = 1L;

//    /**
//     * Argument collection constructor which defines tool default annotation arguments when they are not overridden by the user
//     *
//     * NOTE: This is a temporary solution to allow configurable annotations before the annotation plugin descriptor has been
//     * hooked up to the engine as the solution for adding default annotations. This constructor will soon not exist.
//     * Be aware that the new path will involve passing default annotations or groups directly to the GATKAnnotationPluginDescriptor.
//     *
//     * @param defaultGroups               List annotation group names to be used by default
//     * @param defaultAnnotations          List of annotation class names to be used by default
//     * @param defaultAnnotationsToExclude List of annotation class names to exclude by default. These override the default annotations and annotation groups.
//     */
//    public DefaultGATKVariantAnnotationArgumentCollection(final List<String> defaultGroups, final List<String> defaultAnnotations, final List<String> defaultAnnotationsToExclude) {
//        Utils.nonNull(defaultGroups);
//        Utils.nonNull(defaultAnnotations);
//        Utils.nonNull(defaultAnnotationsToExclude);
//
//        annotationGroupsToUse = new ArrayList<>(defaultGroups);
//        annotationsToUse = new ArrayList<>(defaultAnnotations);
//        annotationsToExclude = new ArrayList<>(defaultAnnotationsToExclude);
//    }

    /**
     * Which annotations to include in variant calls in the output. These supplement annotations provided by annotation groups.
     */
    @Argument(fullName = StandardArgumentDefinitions.ANNOTATION_LONG_NAME, shortName = StandardArgumentDefinitions.ANNOTATION_SHORT_NAME, doc = "One or more specific annotations to add to variant calls", optional = true)
    public List<String> annotationsToUse = new ArrayList<>();

    /**
     * Which annotations to exclude from output in the variant calls.  Note that this argument has higher priority than the
     * -A or -G arguments, so these annotations will be excluded even if they are explicitly included with the other
     * options.
     */
    @Argument(fullName = StandardArgumentDefinitions.ANNOTATIONS_TO_EXCLUDE_LONG_NAME, shortName = StandardArgumentDefinitions.ANNOTATIONS_TO_EXCLUDE_SHORT_NAME, doc = "One or more specific annotations to exclude from variant calls", optional = true)
    public List<String> annotationsToExclude = new ArrayList<>();

    /**
     * Which groups of annotations to add to the output variant calls.
     * Any requirements that are not met (e.g. failing to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Argument(fullName = StandardArgumentDefinitions.ANNOTATION_GROUP_LONG_NAME, shortName = StandardArgumentDefinitions.ANNOTATION_GROUP_SHORT_NAME, doc = "One or more groups of annotations to apply to variant calls", optional = true)
    public List<String> annotationGroupsToUse = new ArrayList<>();

    /**
     * Hook allowing for the user to remove default annotations from the tool
     */
    @Argument(fullName = StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS, shortName = StandardArgumentDefinitions.DISABLE_TOOL_DEFAULT_ANNOTATIONS, doc = "Disable all tool default annotations", optional = true)
    public boolean disableToolDefaultAnnotations = false;

    /**
     * You can use the -AX argument in combination with this one to exclude specific annotations. Note that some
     * annotations may not be actually applied if they are not applicable to the data provided or if they are
     * unavailable to the tool (e.g. there are several annotations that are currently not hooked up to
     * HaplotypeCaller). At present no error or warning message will be provided, the annotation will simply be
     * skipped silently. You can check the output VCF header to see which annotations were actually applied (although
     * this does not guarantee that the annotation was applied to all records in the VCF, since some annotations have
     * additional requirements, e.g. minimum number of samples or heterozygous sites only -- see the documentation
     * for individual annotations' requirements).
     */
    @Argument(fullName="use-all-annotations", doc="Use all possible annotations (not for the faint of heart)", optional=true)
    protected Boolean USE_ALL_ANNOTATIONS = false;

    @Override
    public List<String> getUserEnabledAnnotationNames() {
        return annotationsToUse;
    }

    @Override
    public List<String> getUserEnabledAnnotationGroups() {
        return annotationGroupsToUse;
    }

    @Override
    public List<String> getUserDisabledAnnotationNames() {
        return annotationsToExclude;
    }

    @Override
    public boolean getDisableToolDefaultAnnotations() {
        return disableToolDefaultAnnotations;
    }

    @Override
    public boolean getEnableAllAnnotations() {
        return USE_ALL_ANNOTATIONS;
    }
}
