package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Tool to annotate targets.
 * <p>
 *     This tools adds new columns to an existing target table with various target {@link TargetAnnotation annotations}.
 * </p>
 * <p>
 *     Each annotation may require some source of information. For example {@link TargetAnnotation#GC_CONTENT}
 *     requires a reference.
 * </p>
 * <p>
 *     Only annotation for which all data source dependencies are satisfied will be added to the ones present in the input
 *     target table into the output target table.
 * </p>
 * <p>
 *     There must be at least one annotation whose dependencies are satisfied, otherwise the
 *     tool will return an error.
 * </p>
 * <p>
 *     These are the currently supported annotations and their dependencies:
 *     <dl>
 *         <dt>{@link TargetAnnotation#GC_CONTENT GC_CONTENT}</dt>
 *         <dd><p>GC content in the target expressed as the fraction of unambiguous
 *         bases (A, C, G, T) that are G or C in the reference sequence that
 *         overlap the target.</p>
 *         <p> It can take values from 0 to 1
 *         or NaN in the case special case that there is no unambiguous base
 *         in the target region.</p>
 *         <p>
 *         This annotation requires a reference that the user can pass using the
 *         {@value StandardArgumentDefinitions#REFERENCE_LONG_NAME} argument.
 *         </p>
 *         </dd>
 *     </dl>
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        oneLineSummary = "Annotate targets with various properties, such as GC content",
        summary = "Annotate targets with various properties, such as GC content",
        programGroup = CopyNumberProgramGroup.class
)
public class AnnotateTargets extends TargetWalker {

    @Argument(
            doc = "Output annotated target file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional  = false
    )
    public File outputFile;

    private TargetWriter outputWriter;

    private Map<TargetAnnotation, TargetAnnotator> annotators;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        annotators = new LinkedHashMap<>(2);
        if (hasReference()) {
            annotators.put(TargetAnnotation.GC_CONTENT, new GCContentAnnotator());
            logger.info(String.format("Adding the %s annotation to the output; a reference has been provided", TargetAnnotation.GC_CONTENT));
        } else {
            logger.info(String.format("Omitting the %s annotation to the output; no reference was provided", TargetAnnotation.GC_CONTENT));
        }


        if (annotators.isEmpty()) {
            throw new UserException.BadInput("Resources needed to perform annotation are missing.");
        }
        try {
            outputWriter = new TargetWriter(outputFile, annotators.keySet());
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try {
            outputWriter.close();
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "problems closing the output file");
        }
        return super.onTraversalSuccess();
    }

    @Override
    public void apply(final Target target, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final TargetAnnotationCollection outputAnnotations = new TargetAnnotationCollection();
        // Add the input annotations:
        final TargetAnnotationCollection inputAnnotations = target.getAnnotations();
        for (final TargetAnnotation annotation : inputAnnotations.annotationSet()) {
            outputAnnotations.put(annotation, inputAnnotations.get(annotation));
        }
        // Calculate and add the requested annotations:
        for (final Map.Entry<TargetAnnotation, TargetAnnotator> annotatorsEntry : annotators.entrySet()) {
            outputAnnotations.put(annotatorsEntry.getKey(),
                    String.valueOf(annotatorsEntry.getValue().apply(target, readsContext, referenceContext, featureContext)));
        }
        // Compose the new target annotations, the new target and write it to the output.
        final Target newTarget = new Target(target.getName(), target.getInterval(), outputAnnotations);
        try {
            outputWriter.writeRecord(newTarget);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "problem writing a target in the output");
        }
    }

    private interface TargetAnnotator {
        Object apply(final Target target, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext);
    }

    private class GCContentAnnotator implements TargetAnnotator {
        @Override
        public Double apply(final Target target, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
            final Nucleotide.Counter counter = new Nucleotide.Counter();
            counter.addAll(referenceContext.getBases());
            final long gcCount = counter.get(Nucleotide.C) + counter.get(Nucleotide.G);
            final long atCount = counter.get(Nucleotide.A) + counter.get(Nucleotide.T);
            final long totalCount = gcCount + atCount;
            return totalCount == 0 ? Double.NaN : gcCount / (double) totalCount;
        }
    }
}
