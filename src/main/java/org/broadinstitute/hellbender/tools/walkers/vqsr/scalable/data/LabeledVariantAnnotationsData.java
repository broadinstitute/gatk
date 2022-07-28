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
    // TODO this could be exposed
    private static final int CHUNK_DIVISOR = 16;
    private static final int MAXIMUM_CHUNK_SIZE = HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX / CHUNK_DIVISOR;

    private static final int INITIAL_SIZE = 10_000_000;

    public static final String TRAINING_LABEL = "training";
    public static final String CALIBRATION_LABEL = "calibration";
    public static final String SNP_LABEL = "snp";

    public static final String INTERVALS_PATH = "/intervals";
    public static final String ALLELES_REF_PATH = "/alleles/ref";
    public static final String ALLELES_ALT_PATH = "/alleles/alt";
    public static final String ANNOTATIONS_NAMES_PATH = "/annotations/names";
    public static final String ANNOTATIONS_PATH = "/annotations";
    public static final String LABELS_PATH = "/labels";
    public static final String LABELS_SNP_PATH = LABELS_PATH + "/snp";

    private final List<String> sortedAnnotationNames;
    final List<String> sortedLabels;

    private final List<List<LabeledVariantAnnotationsDatum>> data;
    private final boolean useASAnnotations;

    public LabeledVariantAnnotationsData(final Collection<String> annotationNames,
                                         final Collection<String> labels,
                                         final boolean useASAnnotations,
                                         final int initialSize) {
        data = new ArrayList<>(initialSize);
        sortedAnnotationNames = ImmutableList.copyOf(annotationNames.stream().distinct().sorted().collect(Collectors.toList()));
        Utils.validateArg(sortedAnnotationNames.size() > 0, "Number of annotation names must be positive.");
        if (sortedAnnotationNames.size() != annotationNames.size()) {
            logger.warn(String.format("Ignoring duplicate annotations: %s.", Utils.getDuplicatedItems(annotationNames)));
        }
        sortedLabels = ImmutableList.copyOf(labels.stream().distinct().sorted().collect(Collectors.toList()));
        if (sortedLabels.size() != labels.size()) {
            logger.warn(String.format("Ignoring duplicate labels: %s.", Utils.getDuplicatedItems(labels)));
        }
        this.useASAnnotations = useASAnnotations;
    }

    public LabeledVariantAnnotationsData(final Collection<String> annotationNames,
                                         final Collection<String> labels,
                                         final boolean useASAnnotations) {
        this(annotationNames, labels, useASAnnotations, INITIAL_SIZE);
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

    /**
     * Adds an element to the underlying {@link #data} collection.
     */
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

    /**
     * Sets the element at a specified index in the underlying {@link #data} collection.
     */
    public void set(final int index,
                    final VariantContext vc,
                    final List<List<Allele>> altAllelesPerDatum,
                    final List<VariantType> variantTypePerDatum,
                    final List<TreeSet<String>> labelsPerDatum) {
        if (!useASAnnotations) {
            data.set(index, Collections.singletonList(new LabeledVariantAnnotationsDatum(
                    vc, altAllelesPerDatum.get(0), variantTypePerDatum.get(0), labelsPerDatum.get(0), sortedAnnotationNames, useASAnnotations)));
        } else {
            data.set(index, IntStream.range(0, altAllelesPerDatum.size()).boxed()
                    .map(i -> new LabeledVariantAnnotationsDatum(
                            vc, altAllelesPerDatum.get(i), variantTypePerDatum.get(i), labelsPerDatum.get(i), sortedAnnotationNames, useASAnnotations))
                    .collect(Collectors.toList()));
        }
    }

    /**
     * @return  list of {@link VariantType} indicators, with length given by the number of corresponding sites
     */
    public List<VariantType> getVariantTypeFlat() {
        return streamFlattenedData().map(datum -> datum.variantType).collect(Collectors.toList());
    }

    /**
     * @return  list of boolean label indicators, with length given by the number of sites;
     *          an element in the list will be true if the corresponding site is assigned to the specified label
     */
    public List<Boolean> isLabelFlat(final String label) {
        return streamFlattenedData().map(datum -> datum.labels.contains(label)).collect(Collectors.toList());
    }

    private Stream<LabeledVariantAnnotationsDatum> streamFlattenedData() {
        return data.stream().flatMap(List::stream);
    }

    /**
     * Writes a representation of the collection to an HDF5 file with the following directory structure:
     *
     *      <p>
     *          |--- alleles<br>
     *          │    |--- alt<br>
     *          │    |--- ref<br>
     *          |--- annotations<br>
     *          │    |--- chunk_0<br>
     *          │    |--- ...<br>
     *          │    |--- chunk_{num_chunks - 1}<br>
     *          │    |--- names<br>
     *          │    |--- num_chunks<br>
     *          │    |--- num_columns<br>
     *          │    |--- num_rows<br>
     *          |--- intervals<br>
     *          │    |--- indexed_contig_names<br>
     *          │    |--- transposed_index_start_end<br>
     *          |--- labels<br>
     *          │    |--- snp<br>
     *          │    |--- ... (e.g., training, calibration, etc.)<br>
     *          │    |--- ...<br>
     *      </p>
     *
     * Here, each chunk is a double matrix, with dimensions given by (number of sites in the chunk) x (number of annotations).
     * See the methods {@link HDF5Utils#writeChunkedDoubleMatrix} and {@link HDF5Utils#writeIntervals} for additional details.
     *
     * @param omitAllelesInHDF5 string arrays containing ref/alt alleles can be large, so we allow the option of omitting them
     */
    public void writeHDF5(final File outputFile,
                          final boolean omitAllelesInHDF5) {

        try (final HDF5File outputHDF5File = new HDF5File(outputFile, HDF5File.OpenMode.CREATE)) {
            IOUtils.canReadFile(outputHDF5File.getFile());
            HDF5Utils.writeIntervals(outputHDF5File, INTERVALS_PATH,
                    streamFlattenedData().map(datum -> datum.interval).collect(Collectors.toList()));
            if (!omitAllelesInHDF5) {
                outputHDF5File.makeStringArray(ALLELES_REF_PATH,
                        streamFlattenedData().map(datum -> datum.refAllele.getDisplayString()).toArray(String[]::new));
                if (!useASAnnotations) {
                    outputHDF5File.makeStringArray(ALLELES_ALT_PATH,
                            streamFlattenedData()
                                    .map(datum -> datum.altAlleles.stream().map(Allele::getDisplayString).collect(Collectors.joining(",")))
                                    .toArray(String[]::new));
                } else {
                    outputHDF5File.makeStringArray(ALLELES_ALT_PATH,
                            streamFlattenedData().map(datum -> datum.altAlleles.get(0).getDisplayString()).toArray(String[]::new));
                }
            }
            outputHDF5File.makeStringArray(ANNOTATIONS_NAMES_PATH, sortedAnnotationNames.toArray(new String[0]));
            HDF5Utils.writeChunkedDoubleMatrix(outputHDF5File, ANNOTATIONS_PATH,
                    streamFlattenedData().map(datum -> datum.annotations).toArray(double[][]::new), MAXIMUM_CHUNK_SIZE);
            outputHDF5File.makeDoubleArray(LABELS_SNP_PATH,
                    streamFlattenedData().mapToDouble(datum -> datum.variantType == VariantType.SNP ? 1 : 0).toArray());
            for (final String label : sortedLabels) {
                outputHDF5File.makeDoubleArray(String.format("%s/%s", LABELS_PATH, label),
                        streamFlattenedData().mapToDouble(datum -> datum.labels.contains(label) ? 1 : 0).toArray());
            }
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of annotations and metadata (%s). Output file at %s may be in a bad state.",
                    exception, outputFile.getAbsolutePath()));
        }
    }

    /**
     * @return  list of annotation names, with length given by the number of annotations, read from the specified file
     */
    public static List<String> readAnnotationNames(final File annotationsFile) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return Arrays.asList(annotationsHDF5File.readStringArray(ANNOTATIONS_NAMES_PATH));
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of annotation names from %s: %s",
                    annotationsFile.getAbsolutePath(), exception));
        }
    }

    /**
     * @return  matrix with dimensions (number of sites) x (number of annotations), read from the specified file
     */
    public static double[][] readAnnotations(final File annotationsFile) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return HDF5Utils.readChunkedDoubleMatrix(annotationsHDF5File, ANNOTATIONS_PATH);
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of annotations from %s: %s",
                    annotationsFile.getAbsolutePath(), exception));
        }
    }

    /**
     * @return  list of boolean label indicators, with length given by the number of corresponding sites, read from the specified file;
     *          an element in the list will be true if the corresponding site is assigned to the specified label
     */
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
        Utils.validateArg(annotationNames.size() > 0, "Number of annotation names must be positive.");
        Utils.validateArg(allAnnotations.length > 0, "Number of annotation data points must be positive.");
        Utils.validateArg(annotationNames.size() == allAnnotations[0].length,
                "Number of annotation names must match number of features in annotation data.");
        final double[][] subsetData = IntStream.range(0, isSubset.size()).boxed().filter(isSubset::get).map(i -> allAnnotations[i]).toArray(double[][]::new);
        final File subsetAnnotationsFile = IOUtils.createTempFile("subset.annot", ".hdf5");
        try (final HDF5File subsetAnnotationsHDF5File = new HDF5File(subsetAnnotationsFile, HDF5File.OpenMode.CREATE)) {
            subsetAnnotationsHDF5File.makeStringArray(ANNOTATIONS_NAMES_PATH, annotationNames.toArray(new String[0]));
            HDF5Utils.writeChunkedDoubleMatrix(subsetAnnotationsHDF5File, ANNOTATIONS_PATH, subsetData, MAXIMUM_CHUNK_SIZE);
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of annotations (%s). Output file at %s may be in a bad state.",
                    exception, subsetAnnotationsFile.getAbsolutePath()));
        }
        return subsetAnnotationsFile;
    }
}
