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
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/*
 * TODO this whole class needs refactoring. it's cleaned up significantly from VQSR version,
 *  but still has a long way to go. we should decide how strongly to couple to the tranche code and refactor
 *  that at the same time
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

    private final List<LabeledVariantAnnotationsDatum> data;
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
                    final List<Allele> altAlleles,
                    final VariantType variantType,
                    final Set<String> labels) {
        data.add(new LabeledVariantAnnotationsDatum(vc, altAlleles, variantType, labels, sortedAnnotationNames, useASAnnotations));
    }

    private Stream<LabeledVariantAnnotationsDatum> streamFlattenedData() {
        return data.stream();
    }

    public void writeHDF5(final File outputFile,
                          final boolean omitAllelesInHDF5) {
        // TODO validate
        try (final HDF5File outputHDF5File = new HDF5File(outputFile, HDF5File.OpenMode.CREATE)) {
            IOUtils.canReadFile(outputHDF5File.getFile());
            HDF5Utils.writeIntervals(outputHDF5File, "/intervals",
                    streamFlattenedData().map(datum -> datum.interval).collect(Collectors.toList()));
            if (!omitAllelesInHDF5) {
                outputHDF5File.makeStringArray("/alleles/ref",
                        streamFlattenedData().map(datum -> datum.refAllele.getDisplayString()).toArray(String[]::new));
                if (!useASAnnotations) {
                    outputHDF5File.makeStringArray("/alleles/alt",
                            streamFlattenedData()
                                    .map(datum -> datum.altAlleles.stream().map(Allele::getDisplayString).collect(Collectors.joining(",")))
                                    .toArray(String[]::new));
                } else {
                    outputHDF5File.makeStringArray("/alleles/alt",
                            streamFlattenedData().map(datum -> datum.altAlleles.get(0).getDisplayString()).toArray(String[]::new));
                }
            }
            outputHDF5File.makeStringArray("/annotations/names", sortedAnnotationNames.toArray(new String[0]));
            HDF5Utils.writeChunkedDoubleMatrix(outputHDF5File, "/annotations",
                    streamFlattenedData().map(datum -> datum.annotations).toArray(double[][]::new), MAXIMUM_CHUNK_SIZE);
            outputHDF5File.makeDoubleArray("/labels/snp",
                    streamFlattenedData().mapToDouble(datum -> datum.variantType == VariantType.SNP ? 1 : 0).toArray());
            for (final String label : sortedLabels) {
                outputHDF5File.makeDoubleArray(String.format("/labels/%s", label),
                        streamFlattenedData().mapToDouble(datum -> datum.labels.contains(label) ? 1 : 0).toArray());
            }
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of annotations and metadata (%s). Output file at %s may be in a bad state.",
                    exception, outputFile.getAbsolutePath()));
        }
        logger.debug(String.format("Annotations and metadata written to %s.", outputFile.getAbsolutePath()));
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
