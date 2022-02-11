package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

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
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/*
 * TODO this whole class needs refactoring. it's cleaned up significantly from VQSR version,
 *  but still has a long way to go. we should decide how strongly to couple to the tranche code and refactor
 *  that at the same time
 */
final class LabeledVariantAnnotationsData {
    private static final Logger logger = LogManager.getLogger(LabeledVariantAnnotationsData.class);

    private static final int INITIAL_SIZE = 1000000;

    // TODO make labels enum?
    static final String TRAINING_LABEL = "training";
    static final String TRUTH_LABEL = "truth";

    private final List<String> sortedAnnotationKeys;
    final List<String> sortedLabels;
    final boolean useASAnnotations;

    private final List<LabeledVariantAnnotationsDatum> data;
    final List<List<Allele>> alternateAlleles;

    LabeledVariantAnnotationsData(final Collection<String> annotationKeys,
                                  final Collection<String> labels,
                                  final boolean useASAnnotations) {
        this.data = new ArrayList<>(INITIAL_SIZE);
        this.alternateAlleles = new ArrayList<>(INITIAL_SIZE);
        this.sortedAnnotationKeys = ImmutableList.copyOf(annotationKeys.stream().distinct().sorted().collect(Collectors.toList()));
        if (sortedAnnotationKeys.size() != annotationKeys.size()) {
            logger.warn(String.format("Ignoring duplicate annotations: %s.", Utils.getDuplicatedItems(annotationKeys)));
        }
        this.sortedLabels = ImmutableList.copyOf(labels.stream().distinct().sorted().collect(Collectors.toList()));
        if (sortedLabels.size() != labels.size()) {
            logger.warn(String.format("Ignoring duplicate labels: %s.", Utils.getDuplicatedItems(labels)));
        }
        this.useASAnnotations = useASAnnotations;
    }

    void addDatum(final VariantContext vc,
                  final Allele refAllele, 
                  final Allele altAllele,
                  final Set<String> labels) {
        final LabeledVariantAnnotationsDatum datum = new LabeledVariantAnnotationsDatum(vc, refAllele, altAllele, labels, sortedAnnotationKeys, useASAnnotations);
        data.add(datum);
    }

    // TODO clean this up; need to store alleles to enable passing of training sites-only VCF to score tool for marking of training sites
    void addAlternateAlleles(final List<Allele> alternateAlleles) {
        this.alternateAlleles.add(alternateAlleles);
    }

    public List<LabeledVariantAnnotationsDatum> getData() {
        return data;
    }

    List<String> getSortedAnnotationKeys() {
        return sortedAnnotationKeys;
    }

    static void setScores(final List<LabeledVariantAnnotationsDatum> data,
                          final double[] scores) {
        IntStream.range(0, data.size()).forEach(i -> data.get(i).score = scores[i]);
    }

    void writeAnnotationsHDF5(final File outputFile,
                              final int maximumChunkSize) {
        // TODO validate
        try (final HDF5File outputHDF5File = new HDF5File(outputFile, HDF5File.OpenMode.CREATE)) {
            IOUtils.canReadFile(outputHDF5File.getFile());

            outputHDF5File.makeStringArray("/data/annotation_names", sortedAnnotationKeys.toArray(new String[0]));
            HDF5Utils.writeChunkedDoubleMatrix(outputHDF5File, "/data/annotations",
                    data.stream().map(vd -> vd.annotations.stream().mapToDouble(x -> x).toArray()).toArray(double[][]::new),
                    maximumChunkSize);
            for (final String label : sortedLabels) {
                outputHDF5File.makeDoubleArray(String.format("/data/is_%s", label),
                        data.stream().mapToDouble(vd -> vd.labels.contains(label) ? 1 : 0).toArray());
            }
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of annotations (%s). Output file at %s may be in a bad state.",
                    exception, outputFile.getAbsolutePath()));
        }
        logger.info(String.format("Annotations written to %s.", outputFile.getAbsolutePath()));
    }

    static double[][] readData(final File annotationsFile) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return HDF5Utils.readChunkedDoubleMatrix(annotationsHDF5File, "/data/annotations");
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of annotations from %s: %s",
                    annotationsFile.getAbsolutePath(), exception));
        }
    }

    static List<Boolean> readLabel(final File annotationsFile,
                               final String label) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return readBooleanList(annotationsHDF5File, String.format("/data/is_%s", label));
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of label %s from %s: %s",
                    label, annotationsFile.getAbsolutePath(), exception));
        }
    }

    static String[] readAnnotationNames(final File annotationsFile) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return annotationsHDF5File.readStringArray("/data/annotation_names");
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of annotation names from %s: %s",
                    annotationsFile.getAbsolutePath(), exception));
        }
    }

    private static List<Boolean> readBooleanList(final HDF5File hdf5File,
                                                 final String path) {
        return Arrays.stream(hdf5File.readDoubleArray(path))
                .mapToObj(d -> (d == 1))
                .collect(Collectors.toList());
    }
}
