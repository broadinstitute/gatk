package org.broadinstitute.hellbender.tools.walkers.readorientation;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Learn the prior probability of read orientation artifact from the output of {@link CollectF1R2Counts}
 * Details of the model may be found in docs/mutect/mutect.pdf.
 *
 *
 * <h3>Usage Examples</h3>
 *
 * gatk LearnReadOrientationModel \
 *   -alt-table my-tumor-sample-alt.tsv \
 *   -ref-hist my-tumor-sample-ref.metrics \
 *   -alt-hist my-tumor-sample-alt-depth1.metrics \
 *   -O my-tumor-sample-artifact-prior.tsv
 */
@CommandLineProgramProperties(
        summary = "Get the maximum likelihood estimates of artifact prior probabilities in the orientation bias mixture model filter",
        oneLineSummary = "Get the maximum likelihood estimates of artifact prior probabilities in the orientation bias mixture model filter",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class LearnReadOrientationModel extends CommandLineProgram {
    public static final double DEFAULT_CONVERGENCE_THRESHOLD = 1e-4;
    public static final int DEFAULT_MAX_ITERATIONS = 20;

    public static final String EM_CONVERGENCE_THRESHOLD_LONG_NAME = "convergence-threshold";
    public static final String MAX_EM_ITERATIONS_LONG_NAME = "num-em-iterations";
    public static final String MAX_DEPTH_LONG_NAME = "max-depth";

    @Argument(fullName = CollectF1R2Counts.REF_SITE_METRICS_LONG_NAME, doc = "histograms of depths over ref sites for each reference context")
    private File refHistogramTable;

    @Argument(fullName = CollectF1R2Counts.ALT_DATA_TABLE_LONG_NAME,  doc = "a table of F1R2 and depth counts")
    private File altDataTable;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "table of artifact priors")
    private File output;

    @Argument(fullName = CollectF1R2Counts.ALT_DEPTH1_HISTOGRAM_LONG_NAME, doc = "histograms of depth 1 alt sites", optional = true)
    private File altHistogramTable = null;

    @Argument(fullName = EM_CONVERGENCE_THRESHOLD_LONG_NAME, doc = "Stop the EM when the distance between parameters between iterations falls below this value", optional = true)
    private double converagenceThreshold = DEFAULT_CONVERGENCE_THRESHOLD;

    @Argument(fullName = MAX_EM_ITERATIONS_LONG_NAME, doc = "give up on EM after this many iterations", optional = true)
    private int maxEMIterations = DEFAULT_MAX_ITERATIONS;

    @Argument(fullName = MAX_DEPTH_LONG_NAME, doc = "sites with depth higher than this value will be grouped", optional = true)
    private int maxDepth = F1R2FilterConstants.DEFAULT_MAX_DEPTH;

    List<Histogram<Integer>> refHistograms;

    List<Histogram<Integer>> altHistograms;

    final ArtifactPriorCollection artifactPriorCollection = new ArtifactPriorCollection();;

    @Override
    protected void onStartup(){
        final MetricsFile<?, Integer> referenceSiteMetrics = readMetricsFile(refHistogramTable);
        refHistograms = referenceSiteMetrics.getAllHistograms();

        if (altHistogramTable != null) {
            final MetricsFile<?, Integer> altSiteMetrics = readMetricsFile(altHistogramTable);
            altHistograms = altSiteMetrics.getAllHistograms();
        } else {
            altHistograms = Collections.emptyList();
        }
    }

    @Override
    public Object doWork(){
        final int defaultInitialListSize = 1_000_000;

        final Map<String, List<AltSiteRecord>> altDesignMatrixByContext =
                AltSiteRecord.readAltSiteRecords(altDataTable, defaultInitialListSize).stream()
                        .collect(Collectors.groupingBy(AltSiteRecord::getReferenceContext));

        // Since e.g. G->T under AGT F1R2 is equivalent to C->A under ACT F2R1, combine the data
        for (final String refContext : F1R2FilterConstants.CANONICAL_KMERS){
            final String reverseComplement = SequenceUtil.reverseComplement(refContext);

            // Merge ref histograms
            final Histogram<Integer> refHistogram = refHistograms.stream()
                    .filter(h -> h.getValueLabel().equals(refContext))
                    .findFirst().orElseGet(() -> F1R2FilterUtils.createRefHistogram(refContext, maxDepth));
            final Histogram<Integer> refHistogramRevComp = refHistograms.stream()
                    .filter(h -> h.getValueLabel().equals(reverseComplement))
                    .findFirst().orElseGet(() -> F1R2FilterUtils.createRefHistogram(reverseComplement, maxDepth));
            final Histogram<Integer> combinedRefHistograms = combineRefHistogramWithRC(refContext, refHistogram, refHistogramRevComp, maxDepth);


            // Merge alt depth=1 histograms
            final List<Histogram<Integer>> altDepthOneHistogramsForContext = altHistograms.stream()
                    .filter(h -> h.getValueLabel().startsWith(refContext))
                    .collect(Collectors.toList());
            final List<Histogram<Integer>> altDepthOneHistogramsRevComp = altHistograms.stream()
                    .filter(h -> h.getValueLabel().startsWith(reverseComplement))
                    .collect(Collectors.toList());
            final List<Histogram<Integer>> combinedAltHistograms = combineAltDepthOneHistogramWithRC(altDepthOneHistogramsForContext, altDepthOneHistogramsRevComp, maxDepth);

            // Finally, merge the rest of alt records
            final List<AltSiteRecord> altDesignMatrix = altDesignMatrixByContext.getOrDefault(refContext, new ArrayList<>()); // Cannot use Collections.emptyList() here because the input list must be mutable
            final List<AltSiteRecord> altDesignMatrixRevComp = altDesignMatrixByContext.getOrDefault(reverseComplement, Collections.emptyList());
            // Warning: the below method will mutate the content of {@link altDesignMatrixRevComp} and append to {@code altDesignMatrix}
            mergeDesignMatrices(altDesignMatrix, altDesignMatrixRevComp);


            if (combinedRefHistograms.getSumOfValues() == 0 || altDesignMatrix.isEmpty()) {
                logger.info(String.format("Skipping the reference context %s as we didn't find either the ref or alt table for the context", refContext));
                continue;
            }

            final LearnReadOrientationModelEngine engine = new LearnReadOrientationModelEngine(
                    combinedRefHistograms,
                    combinedAltHistograms,
                    altDesignMatrix,
                    converagenceThreshold,
                    maxEMIterations,
                    maxDepth,
                    logger);
            final ArtifactPrior artifactPrior = engine.learnPriorForArtifactStates();
            artifactPriorCollection.set(artifactPrior);
        }

        artifactPriorCollection.writeArtifactPriors(output);
        return "SUCCESS";
    }

    @VisibleForTesting
    public static Histogram<Integer> combineRefHistogramWithRC(final String refContext,
                                                               final Histogram<Integer> refHistogram,
                                                               final Histogram<Integer> refHistogramRevComp,
                                                               final int maxDepth){
        Utils.validateArg(refHistogram.getValueLabel()
                .equals(SequenceUtil.reverseComplement(refHistogramRevComp.getValueLabel())),
                "ref context = " + refHistogram.getValueLabel() + ", rev comp = " + refHistogramRevComp.getValueLabel());
        Utils.validateArg(refHistogram.getValueLabel().equals(refContext), "this better match");

        final Histogram<Integer> combinedRefHistogram = F1R2FilterUtils.createRefHistogram(refContext, maxDepth);

        for (final Integer depth : refHistogram.keySet()){
            final double newCount = refHistogram.get(depth).getValue() + refHistogramRevComp.get(depth).getValue();
            combinedRefHistogram.increment(depth, newCount);
        }

        return combinedRefHistogram;
    }

    @VisibleForTesting
    public static List<Histogram<Integer>> combineAltDepthOneHistogramWithRC(final List<Histogram<Integer>> altHistograms,
                                                                             final List<Histogram<Integer>> altHistogramsRevComp,
                                                                             final int maxDepth){
        if (altHistograms.isEmpty() && altHistogramsRevComp.isEmpty()){
            return Collections.emptyList();
        }

        final String refContext = ! altHistograms.isEmpty() ?
                F1R2FilterUtils.labelToTriplet(altHistograms.get(0).getValueLabel()).getLeft() :
                SequenceUtil.reverseComplement(F1R2FilterUtils.labelToTriplet(altHistogramsRevComp.get(0).getValueLabel()).getLeft());

        // Contract: altHistogram must be of the canonical representation of the kmer
        Utils.validateArg(F1R2FilterConstants.CANONICAL_KMERS.contains(refContext), "refContext must be the canonical representation but got " + refContext);

        final List<Histogram<Integer>> combinedHistograms = new ArrayList<>(F1R2FilterConstants.numAltHistogramsPerContext);

        for (Nucleotide altAllele : Nucleotide.STANDARD_BASES){
            // Skip when the alt base is the ref base, which doesn't make sense because this is a histogram of alt sites
            if (altAllele == F1R2FilterUtils.getMiddleBase(refContext)){
                continue;
            }

            final String reverseComplement = SequenceUtil.reverseComplement(refContext);
            final Nucleotide altAlleleRevComp = Nucleotide.valueOf(SequenceUtil.reverseComplement(altAllele.toString()));

            for (ReadOrientation orientation : ReadOrientation.values()) {
                final ReadOrientation otherOrientation = ReadOrientation.getOtherOrientation(orientation);
                final Histogram<Integer> altHistogram = altHistograms.stream()
                        .filter(h -> h.getValueLabel().equals(F1R2FilterUtils.tripletToLabel(refContext, altAllele, orientation)))
                        .findFirst().orElseGet(() -> F1R2FilterUtils.createAltHistogram(refContext, altAllele, orientation, maxDepth));

                final Histogram<Integer> altHistogramRevComp = altHistogramsRevComp.stream()
                        .filter(h -> h.getValueLabel().equals(F1R2FilterUtils.tripletToLabel(reverseComplement, altAlleleRevComp, otherOrientation)))
                        .findFirst().orElseGet(() -> F1R2FilterUtils.createAltHistogram(reverseComplement, altAlleleRevComp, otherOrientation, maxDepth));

                final Histogram<Integer> combinedHistogram = F1R2FilterUtils.createAltHistogram(refContext, altAllele, orientation, maxDepth);

                // Add the histograms manually - I don't like the addHistogram() in htsjdk method because it does so with side-effect
                for (final Integer depth : altHistogram.keySet()){
                    final double newCount = altHistogram.get(depth).getValue() + altHistogramRevComp.get(depth).getValue();
                    combinedHistogram.increment(depth, newCount);
                }
                combinedHistograms.add(combinedHistogram);
            }
        }

        return combinedHistograms;
    }


    /**
     * Contract: this method must be called after grouping the design matrices by context.
     * That is, {@param altDesignMatrix} must be a list of {@link AltSiteRecord} of a single reference context
     * (which is in {@link F1R2FilterConstants.CANONICAL_KMERS}) and {@param altDesignRevComp} contains only
     * {@link AltSiteRecord} of its reverse complement.
     */
    @VisibleForTesting
    public static void mergeDesignMatrices(final List<AltSiteRecord> altDesignMatrix, List<AltSiteRecord> altDesignMatrixRevComp){
        if (altDesignMatrix.isEmpty() && altDesignMatrixRevComp.isEmpty()){
            return;
        }

        // Order matters here. Assumes that all elements in the list have the same reference context
        Utils.validateArg(altDesignMatrix.isEmpty() || F1R2FilterConstants.CANONICAL_KMERS.contains(altDesignMatrix.get(0).getReferenceContext()),
                "altDesignMatrix must have the canonical representation");

        final Optional<String> refContext = altDesignMatrix.isEmpty() ? Optional.empty() :
                Optional.of(altDesignMatrix.get(0).getReferenceContext());
        final Optional<String> revCompContext = altDesignMatrixRevComp.isEmpty() ? Optional.empty() :
                Optional.of(altDesignMatrixRevComp.get(0).getReferenceContext());

        // If the matrices aren't empty, their reference context much be the reverse complement of each other
        if (refContext.isPresent() && revCompContext.isPresent()){
            Utils.validateArg(refContext.get().equals(SequenceUtil.reverseComplement(revCompContext.get())),
                    "ref context and its rev comp don't match");
        }

        altDesignMatrix.addAll(altDesignMatrixRevComp.stream().map(AltSiteRecord::getReverseComplementOfRecord).collect(Collectors.toList()));
    }

    private MetricsFile<?, Integer> readMetricsFile(File file){
        final MetricsFile<?, Integer> metricsFile = new MetricsFile<>();
        final Reader reader = IOUtil.openFileForBufferedReading(file);
        metricsFile.read(reader);
        CloserUtil.close(reader);
        return metricsFile;
    }

}
