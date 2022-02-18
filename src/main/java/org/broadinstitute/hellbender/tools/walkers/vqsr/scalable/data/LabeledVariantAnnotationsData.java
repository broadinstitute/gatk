package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data;

import com.google.common.collect.ImmutableList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsModel;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.modeling.VariantAnnotationsScorer;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Represents a collection of {@link LabeledVariantAnnotationsDatum} as a list of lists of datums.
 * The outer list is always per-variant. In allele-specific mode, each datum in the inner lists
 * corresponds to a single allele; otherwise, each inner list trivially contains a single datum corresponding
 * to the variant.
 */
public final class LabeledVariantAnnotationsData {
    private static final Logger logger = LogManager.getLogger(LabeledVariantAnnotationsData.class);

    // chunk size in temporary annotation files
    private static final int CHUNK_DIVISOR = 16;
    private static final int MAXIMUM_CHUNK_SIZE = HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX / CHUNK_DIVISOR;

    private static final int INITIAL_SIZE = 10_000_000;

    // TODO make labels enum?
    public static final String TRAINING_LABEL = "training";
    public static final String TRUTH_LABEL = "truth";

    private final List<String> sortedAnnotationNames;
    final List<String> sortedLabels;

    private final List<List<LabeledVariantAnnotationsDatum>> data;
    private final boolean useASAnnotations;

    public LabeledVariantAnnotationsData(final Collection<String> annotationNames,
                                         final Collection<String> labels,
                                         final boolean useASAnnotations) {
        this.data = new ArrayList<>(INITIAL_SIZE);
        this.sortedAnnotationNames = ImmutableList.copyOf(annotationNames.stream().distinct().sorted().collect(Collectors.toList()));
        if (sortedAnnotationNames.size() != annotationNames.size()) {
            logger.warn(String.format("Ignoring duplicate annotations: %s.", Utils.getDuplicatedItems(annotationNames)));
        }
        this.sortedLabels = ImmutableList.copyOf(labels.stream().distinct().sorted().collect(Collectors.toList()));
        if (sortedLabels.size() != labels.size()) {
            logger.warn(String.format("Ignoring duplicate labels: %s.", Utils.getDuplicatedItems(labels)));
        }
        this.useASAnnotations = useASAnnotations;
    }

    public List<String> getSortedAnnotationNames() {
        return sortedAnnotationNames;
    }

    public List<String> getSortedLabels() {
        return sortedLabels;
    }

    public int size() {
        return data.size();
    }

    public void clear() {
        data.clear();
    }

    public void add(final VariantContext vc,
                    final List<List<Allele>> altAllelesPerDatum,
                    final List<VariantType> variantTypePerDatum,
                    final List<TreeSet<String>> labelsPerDatum) {
        if (!useASAnnotations) {
            data.add(Collections.singletonList(new LabeledVariantAnnotationsDatum(
                    vc, altAllelesPerDatum.get(0), variantTypePerDatum.get(0), labelsPerDatum.get(0), sortedAnnotationNames, useASAnnotations)));
        } else {
            data.add(IntStream.range(0, altAllelesPerDatum.size()).boxed()
                    .map(i -> new LabeledVariantAnnotationsDatum(
                            vc, altAllelesPerDatum.get(i), variantTypePerDatum.get(i), labelsPerDatum.get(i), sortedAnnotationNames, useASAnnotations))
                    .collect(Collectors.toList()));
        }
    }

    public List<VariantType> getVariantTypeFlat() {
        return streamFlattenedData().map(datum -> datum.variantType).collect(Collectors.toList());
    }

    public List<Boolean> isLabelFlat(final String label) {
        return streamFlattenedData().map(datum -> datum.labels.contains(label)).collect(Collectors.toList());
    }

    private Stream<LabeledVariantAnnotationsDatum> streamFlattenedData() {
        return data.stream().flatMap(List::stream);
    }

    // TODO remove this
    static void logHeapUsage(final Logger logger,
                             final String message) {
        final int mb = 1024 * 1024;
        final Runtime runtime = Runtime.getRuntime();
        logger.debug("Used memory [MB]: " + (runtime.totalMemory() - runtime.freeMemory()) / mb);
        logger.debug(message);
    }

    /**
     * Writes a representation of the collection to an HDF5 file with the following directory structure:
     *
     *   ├── alleles
     *   │   ├── alt
     *   │   └── ref
     *   ├── annotations
     *   │   ├── chunk_0
     *   │   ├── ...
     *   │   ├── chunk_{num_chunks - 1}
     *   │   ├── names
     *   │   ├── num_chunks
     *   │   ├── num_columns
     *   │   └── num_rows
     *   ├── intervals
     *   │   ├── indexed_contig_names
     *   │   └── transposed_index_start_end
     *   └── labels
     *       ├── snp
     *       ├── ... (e.g., training, truth, etc.)
     *       └── ...
     *
     * TODO extract constants for paths here and in methods below
     *
     * See the methods {@link HDF5Utils#writeChunkedDoubleMatrix} and {@link HDF5Utils#writeIntervals} for additional details.
     * @param omitAllelesInHDF5 string arrays containing ref/alt alleles can be large, so we allow the option of omitting them
     */
    public void writeHDF5(final File outputFile,
                          final boolean omitAllelesInHDF5) {
        // TODO validate
        try (final HDF5File outputHDF5File = new HDF5File(outputFile, HDF5File.OpenMode.CREATE)) {
            IOUtils.canReadFile(outputHDF5File.getFile());
            logHeapUsage(logger, "writing intervals");
            HDF5Utils.writeIntervals(outputHDF5File, "/intervals",
                    streamFlattenedData().map(datum -> datum.interval).collect(Collectors.toList()));
            if (!omitAllelesInHDF5) {
                logHeapUsage(logger, "writing ref alleles");
                outputHDF5File.makeStringArray("/alleles/ref",
                        streamFlattenedData().map(datum -> datum.refAllele.getDisplayString()).toArray(String[]::new));
                if (!useASAnnotations) {
                    logHeapUsage(logger, "writing alt alleles");
                    outputHDF5File.makeStringArray("/alleles/alt",
                            streamFlattenedData()
                                    .map(datum -> datum.altAlleles.stream().map(Allele::getDisplayString).collect(Collectors.joining(",")))
                                    .toArray(String[]::new));
                } else {
                    logHeapUsage(logger, "writing alt alleles");
                    outputHDF5File.makeStringArray("/alleles/alt",
                            streamFlattenedData().map(datum -> datum.altAlleles.get(0).getDisplayString()).toArray(String[]::new));
                }
            }
            outputHDF5File.makeStringArray("/annotations/names", sortedAnnotationNames.toArray(new String[0]));
            logHeapUsage(logger, "writing annotations");
            HDF5Utils.writeChunkedDoubleMatrix(outputHDF5File, "/annotations",
                    streamFlattenedData().map(datum -> datum.annotations).toArray(double[][]::new), MAXIMUM_CHUNK_SIZE);
            logHeapUsage(logger, "writing snp");
            outputHDF5File.makeDoubleArray("/labels/snp",
                    streamFlattenedData().mapToDouble(datum -> datum.variantType == VariantType.SNP ? 1 : 0).toArray());
            logHeapUsage(logger, "writing labels");
            for (final String label : sortedLabels) {
                outputHDF5File.makeDoubleArray(String.format("/labels/%s", label),
                        streamFlattenedData().mapToDouble(datum -> datum.labels.contains(label) ? 1 : 0).toArray());
            }
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of annotations and metadata (%s). Output file at %s may be in a bad state.",
                    exception, outputFile.getAbsolutePath()));
        }
    }

    public static List<String> readAnnotationNames(final File annotationsFile) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return Arrays.asList(annotationsHDF5File.readStringArray("/annotations/names"));
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of annotation names from %s: %s",
                    annotationsFile.getAbsolutePath(), exception));
        }
    }

    public static double[][] readAnnotations(final File annotationsFile) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return HDF5Utils.readChunkedDoubleMatrix(annotationsHDF5File, "/annotations");
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of annotations from %s: %s",
                    annotationsFile.getAbsolutePath(), exception));
        }
    }

    public static List<Boolean> readLabel(final File annotationsFile,
                                          final String label) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return Arrays.stream(annotationsHDF5File.readDoubleArray(String.format("/labels/%s", label))).boxed().map(d -> d == 1).collect(Collectors.toList());
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of label %s from %s: %s",
                    label, annotationsFile.getAbsolutePath(), exception));
        }
    }

    /**
     * Subsets annotation data according to a boolean filter and writes a limited representation to a temporary HDF5 file.
     * Intended for passing annotations via the file interfaces of {@link VariantAnnotationsModel} and {@link VariantAnnotationsScorer}.
     */
    public static File subsetAnnotationsToTemporaryFile(final List<String> annotationNames,
                                                         final double[][] allAnnotations,
                                                         final List<Boolean> isSubset) {
        final double[][] subsetData = IntStream.range(0, isSubset.size()).boxed().filter(isSubset::get).map(i -> allAnnotations[i]).toArray(double[][]::new);
        final File subsetAnnotationsFile = IOUtils.createTempFile("subset.annot", ".hdf5");
        try (final HDF5File subsetAnnotationsHDF5File = new HDF5File(subsetAnnotationsFile, HDF5File.OpenMode.CREATE)) {
            subsetAnnotationsHDF5File.makeStringArray("/annotations/names", annotationNames.toArray(new String[0]));
            HDF5Utils.writeChunkedDoubleMatrix(subsetAnnotationsHDF5File, "/annotations", subsetData, MAXIMUM_CHUNK_SIZE);
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of annotations (%s). Output file at %s may be in a bad state.",
                    exception, subsetAnnotationsFile.getAbsolutePath()));
        }
        return subsetAnnotationsFile;
    }
}
