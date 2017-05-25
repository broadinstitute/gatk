package org.broadinstitute.hellbender.tools.exome;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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
 * <h3>Examples</h3>

 * <pre>
 * java -Xmx4g -jar $gatk_jar AnnotateTargets \
 *   --targets targets.tsv \
 *   --reference ref_fasta.fa \
 *   --output targets.annotated.tsv
 * </pre>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        oneLineSummary = "Annotate targets with various properties, such as GC content",
        summary = "Annotate targets with various properties, such as GC content",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public class AnnotateTargets extends TargetWalker {

    @Argument(
            doc = "Output annotated target file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional  = false
    )
    protected File outputFile;

    @Argument(
            doc = "Bait file. Should be formatted as a tab-separated table (.tsv) with the following header columns:" +
                    " contig, start, stop, name.",
            fullName = ExomeStandardArgumentDefinitions.BAITS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.BAITS_FILE_SHORT_NAME,
            optional = true
    )
    protected FeatureInput<Target> baits;

    private EnumMap<TargetAnnotation, TargetAnnotator> annotators;

    /**
     * A target list containing traversed targets along with "immediate" (i.e. without requiring post-processing)
     * annotations
     */
    private List<Target> immediatelyAnnotatedTargetList;

    /**
     * Raw data for lazy annotators -- further processing is required after the traversal is complete
     */
    private List<EnumMap<TargetAnnotation, Object>> lazyAnnotationsRawData;


    private boolean hasLazyAnnotators;

    private TargetWriter outputWriter;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        annotators = new EnumMap<>(TargetAnnotation.class);

        /* GC content annotation */
        if (hasReference()) {
            annotators.put(TargetAnnotation.GC_CONTENT, new GCContentAnnotator());
            logger.info(String.format("Adding the %s annotation to the output; a reference has been provided",
                    TargetAnnotation.GC_CONTENT));
        } else {
            logger.info(String.format("Omitting the %s annotation to the output; no reference was provided",
                    TargetAnnotation.GC_CONTENT));
        }

        /* Bait count annotation */
        if (baits != null) {
            annotators.put(TargetAnnotation.BAIT_COUNT, new BaitCountAnnotator());
            logger.info(String.format("Adding the %s annotation to the output; a bait table has been provided",
                    TargetAnnotation.BAIT_COUNT));
        } else {
            logger.info(String.format("Omitting the %s annotation to the output; no bait table was provided",
                    TargetAnnotation.BAIT_COUNT));
        }

        if (annotators.isEmpty()) {
            throw new UserException.BadInput("Resources needed to perform annotation are missing.");
        }

        if (annotators.values().stream().anyMatch(TargetAnnotator::isLazy)) {
            lazyAnnotationsRawData = new ArrayList<>();
            hasLazyAnnotators = true;
        } else {
            hasLazyAnnotators = false;
        }

        try {
            outputWriter = new TargetWriter(outputFile, annotators.keySet());
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }

        /* A container of traversed targets along with their immediate (non-lazy) annotations */
        immediatelyAnnotatedTargetList = new ArrayList<>();
    }

    @Override
    public void apply(final Target target, final ReadsContext readsContext, final ReferenceContext referenceContext,
                      final FeatureContext featureContext) {
        // Calculate and add the requested annotations
        final TargetAnnotationCollection newAnnotations = new TargetAnnotationCollection();
        final EnumMap<TargetAnnotation, Object> lazyAnnotations = new EnumMap<>(TargetAnnotation.class);
        for (final Map.Entry<TargetAnnotation, TargetAnnotator> annotatorsEntry : annotators.entrySet()) {
            final TargetAnnotation annotation = annotatorsEntry.getKey();
            final TargetAnnotator annotator = annotatorsEntry.getValue();
            final Object newAnnotation = annotator.apply(target, readsContext, referenceContext, featureContext);
            if (!annotator.isLazy()) {
                newAnnotations.put(annotation, String.valueOf(newAnnotation));
            } else {
                lazyAnnotations.put(annotation, newAnnotation);
            }
        }
        /* Compose the new target with non-lazy annotations and add to the list */
        final Target newTarget = new Target(target.getName(), target.getInterval(),
                mergeAnnotations(target.getAnnotations(), newAnnotations));
        immediatelyAnnotatedTargetList.add(newTarget);
        if (hasLazyAnnotators) {
            lazyAnnotationsRawData.add(lazyAnnotations);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        /* write new targets */
        logger.info("Writing annotated targets to file...");
        try {
            outputWriter.writeAllRecords(generateFullyAnnotatedTargets());
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "problem writing targets to the output file");
        }

        try {
            outputWriter.close();
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "problems closing the output file");
        }
        return super.onTraversalSuccess();
    }

    /**
     * Merges two {@link TargetAnnotationCollection}s. If both collection intersect on a same set of annotations,
     * the second instance({@code other}) _overrides_ the first instance (@code{original}).
     *
     * @param original first collection
     * @param other second collection
     * @return merged collection
     */
    private static TargetAnnotationCollection mergeAnnotations(final TargetAnnotationCollection original,
                                                               final TargetAnnotationCollection other) {
        final TargetAnnotationCollection outputAnnotations = new TargetAnnotationCollection();
        for (final TargetAnnotation annotation : original.annotationSet()) {
            outputAnnotations.put(annotation, original.get(annotation));
        }
        for (final TargetAnnotation annotation : other.annotationSet()) {
            outputAnnotations.put(annotation, other.get(annotation));
        }
        return outputAnnotations;
    }

    /**
     * Performs lazy annotating operations and generates the fully annotated target list
     *
     * @return list of fully annotated targets ready to be saved
     */
    private List<Target> generateFullyAnnotatedTargets() {
        if (!hasLazyAnnotators) {
            return immediatelyAnnotatedTargetList;
        }
        final EnumMap<TargetAnnotation, List<String>> lazyAnnotationsMap = new EnumMap<>(TargetAnnotation.class);
        for (final TargetAnnotation targetAnnotation : annotators.keySet()) {
            if (annotators.get(targetAnnotation).isLazy()) {
                logger.info(String.format("Finalizing %s annotations...", targetAnnotation.name()));
                final List<Object> rawData = lazyAnnotationsRawData.stream()
                        .map(entry -> entry.get(targetAnnotation))
                        .collect(Collectors.toList());
                lazyAnnotationsMap.put(targetAnnotation, annotators.get(targetAnnotation).processAnnotations(
                        immediatelyAnnotatedTargetList, rawData));
            }
        }
        return IntStream.range(0, immediatelyAnnotatedTargetList.size())
                .mapToObj(targetIndex -> {
                    final Target target = immediatelyAnnotatedTargetList.get(targetIndex);
                    final TargetAnnotationCollection immediateAnnotations = target.getAnnotations();
                    final TargetAnnotationCollection lazyAnnotations = new TargetAnnotationCollection();
                    lazyAnnotationsMap.entrySet().forEach(entry ->
                            lazyAnnotations.put(entry.getKey(), entry.getValue().get(targetIndex)));
                    return new Target(target.getName(), target.getInterval(), mergeAnnotations(immediateAnnotations, lazyAnnotations));
                }).collect(Collectors.toList());
    }

    private interface TargetAnnotator {
        Object apply(final Target target, final ReadsContext readsContext, final ReferenceContext referenceContext,
                     final FeatureContext featureContext);

        /**
         * Indicates whether the annotator is lazy (requires post-processing) or can perform on-the-fly
         */
        boolean isLazy();

        /**
         * Lazy annotators must implement this function
         *
         * @param targets a list of traversed targets
         * @param rawAnnotations a list of collected "raw" annotations
         * @return processed annotations
         */
        default List<String> processAnnotations(final List<Target> targets, final List<Object> rawAnnotations) {
            throw new UnsupportedOperationException("Either this annotator is not lazy, or the annotation post-processing" +
                    " method is not implemented yet");
        }
    }

    private class GCContentAnnotator implements TargetAnnotator {
        @Override
        public Double apply(final Target target, final ReadsContext readsContext, final ReferenceContext referenceContext,
                            final FeatureContext featureContext) {
            final Nucleotide.Counter counter = new Nucleotide.Counter();
            counter.addAll(referenceContext.getBases());
            final long gcCount = counter.get(Nucleotide.C) + counter.get(Nucleotide.G);
            final long atCount = counter.get(Nucleotide.A) + counter.get(Nucleotide.T);
            final long totalCount = gcCount + atCount;
            return totalCount == 0 ? Double.NaN : gcCount / (double) totalCount;
        }

        @Override
        public boolean isLazy() {
            return false;
        }
    }

    private class BaitCountAnnotator implements TargetAnnotator {
        @Override
        public Set<SimpleInterval> apply(final Target target, final ReadsContext readsContext, final ReferenceContext referenceContext,
                                         final FeatureContext featureContext) {
            return featureContext.getValues(baits).stream().map(Target::getInterval).collect(Collectors.toSet());
        }

        @Override
        public boolean isLazy() {
            return true;
        }

        @Override
        @SuppressWarnings("unchecked")
        public List<String> processAnnotations(final List<Target> targets, final List<Object> rawAnnotations) {
            if (targets.isEmpty()) {
                return new ArrayList<>();
            }
            final List<Set<SimpleInterval>> overlappingBaitSetPerTarget = rawAnnotations.stream()
                    .map(rawAnnotation -> (Set<SimpleInterval>) rawAnnotation)
                    .collect(Collectors.toList());
            final double[] baitCounts = calculateBaitCounts(targets, overlappingBaitSetPerTarget);
            return Arrays.stream(baitCounts).mapToObj(String::valueOf).collect(Collectors.toList());
        }
    }

    @VisibleForTesting
    static double[] calculateBaitCounts(List<Target> targets, List<Set<SimpleInterval>> overlappingBaitSetPerTarget) {
        final int numTargets = targets.size();
        Utils.validateArg(overlappingBaitSetPerTarget.size() == numTargets, "The size of the list of overlapping baits" +
                " per target does not match the size of the list of targets");
        final Map<SimpleInterval, Set<Integer>> baitsToTargetIndexMap = new HashMap<>();
        for (int targetIndex = 0; targetIndex < numTargets; targetIndex++) {
            final Set<SimpleInterval> overlappingBaitsSet = overlappingBaitSetPerTarget.get(targetIndex);
            for (final SimpleInterval bait : overlappingBaitsSet) {
                if (!baitsToTargetIndexMap.containsKey(bait)) {
                    baitsToTargetIndexMap.put(bait, new HashSet<>());
                }
                baitsToTargetIndexMap.get(bait).add(targetIndex);
            }
        }

        /* if a bait overlaps with N targets, each target gets a 1/n share of the bait */
        final double[] baitCounts = new double[numTargets];
        baitsToTargetIndexMap.values().forEach(overlappingTargetIndices -> {
            final double sharePerTarget = 1.0 / overlappingTargetIndices.size();
            overlappingTargetIndices.forEach(targetIndex -> baitCounts[targetIndex] += sharePerTarget);
        });
        return baitCounts;
    }
}
