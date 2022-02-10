package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.ImmutableList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

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
final class VariantDataCollection {
    private static final Logger logger = LogManager.getLogger(VariantDataCollection.class);

    private static final int INITIAL_SIZE = 1000000;

    // TODO make labels enum?
    static final String TRAINING_LABEL = "training";
    static final String TRUTH_LABEL = "truth";

    static final String BIALLELIC_SNP_LABEL = "biallelic_snp";
    static final String TRANSITION_LABEL = "transition";
    static final List<String> COMPUTED_LABELS = ImmutableList.of(BIALLELIC_SNP_LABEL, TRANSITION_LABEL);

    private final List<String> sortedAnnotationKeys;
    final List<String> sortedLabels;
    final boolean useASAnnotations;

    private final List<VariantDatum> data;
    final List<List<Allele>> alternateAlleles;

    VariantDataCollection(final Collection<String> annotationKeys,
                          final Collection<String> resourceLabels,
                          final boolean useASAnnotations) {
        this.data = new ArrayList<>(INITIAL_SIZE);
        this.alternateAlleles = new ArrayList<>(INITIAL_SIZE);
        this.sortedAnnotationKeys = ImmutableList.copyOf(annotationKeys.stream().distinct().sorted().collect(Collectors.toList()));
        if (sortedAnnotationKeys.size() != annotationKeys.size()) {
            logger.warn(String.format("Ignoring duplicate annotations: %s.", Utils.getDuplicatedItems(annotationKeys)));
        }
        if (!Collections.disjoint(resourceLabels, COMPUTED_LABELS)) {
            throw new UserException(String.format("Resource labels derived from provided tracks cannot contain the computed labels: %s", COMPUTED_LABELS));
        }
        this.sortedLabels = ImmutableList.copyOf(Stream.concat(resourceLabels.stream(), COMPUTED_LABELS.stream()).distinct().sorted().collect(Collectors.toList()));
        if (sortedLabels.size() != resourceLabels.size() + COMPUTED_LABELS.size()) {
            logger.warn(String.format("Ignoring duplicate resource labels: %s.", Utils.getDuplicatedItems(resourceLabels)));
        }
        this.useASAnnotations = useASAnnotations;
    }

    void addDatum(final VariantContext vc,
                  final Allele refAllele, 
                  final Allele altAllele,
                  final Set<String> labels) {
        final VariantDatum datum = new VariantDatum();

        // Populate the datum with lots of fields from the VariantContext, unfortunately the VC is too big so we just
        // pull in only the things we absolutely need.
        datum.referenceAllele = refAllele;
        datum.alternateAllele = altAllele;
        decodeAnnotations(datum, vc);

        // non-deterministic because order of calls depends on load of machine TODO SL: not sure what this means?
        datum.loc = new SimpleInterval(vc);

        datum.labels = labels;
        if (vc.isSNP() && vc.isBiallelic()) {
            datum.labels.add(BIALLELIC_SNP_LABEL);
            if (GATKVariantContextUtils.isTransition(vc)) {
                datum.labels.add(TRANSITION_LABEL);
            }
        }

        data.add(datum);
    }

    // TODO clean this up; need to store alleles to enable passing of training sites-only VCF to score tool for marking of training sites
    void addAlternateAlleles(final List<Allele> alternateAlleles) {
        this.alternateAlleles.add(alternateAlleles);
    }

    public List<VariantDatum> getData() {
        return data;
    }

    List<String> getSortedAnnotationKeys() {
        return sortedAnnotationKeys;
    }

    private void decodeAnnotations(final VariantDatum datum,
                                   final VariantContext vc) {
        final double[] annotations = new double[sortedAnnotationKeys.size()];
        int iii = 0;
        for(final String key : sortedAnnotationKeys) {
            annotations[iii] = decodeAnnotation(key, vc, useASAnnotations, datum);
            iii++;
        }
        datum.annotations = annotations;
    }

    private static double decodeAnnotation(final String annotationKey,
                                           final VariantContext vc,
                                           final boolean useASannotations,
                                           final VariantDatum datum ) {
        double value;

        try {
            //if we're in allele-specific mode and an allele-specific annotation has been requested, parse the appropriate value from the list
            if (useASannotations && annotationKey.startsWith(GATKVCFConstants.ALLELE_SPECIFIC_PREFIX)) {
                final List<Object> valueList = vc.getAttributeAsList(annotationKey);
                //FIXME: we need to look at the ref allele here too
                if (vc.hasAllele(datum.alternateAllele)) {
                    final int altIndex = vc.getAlleleIndex(datum.alternateAllele) - 1; //- 1 is to convert the index from all alleles (including reference) to just alternate alleles
                    value = Double.parseDouble((String) valueList.get(altIndex));
                } else {
                    //if somehow our alleles got mixed up
                    throw new IllegalStateException("ExtractAnnotationsVariantDatum allele " + datum.alternateAllele + " is not contained in the input VariantContext.");
                }
            } else {
                value = vc.getAttributeAsDouble(annotationKey, Double.NaN);
            }
            if (Double.isInfinite(value)) {
                value = Double.NaN;
            }
        } catch (final NumberFormatException e) {
            value = Double.NaN;
        }
        return value;
    }

    static void setScores(final List<VariantDatum> data,
                          final double[] scores) {
        IntStream.range(0, data.size()).forEach(i -> data.get(i).score = scores[i]);
    }

    void writeAnnotationsHDF5(final File outputFile,
                              final int maximumChunkSize) {
        // TODO validate
        try (final HDF5File outputHDF5File = new HDF5File(outputFile, HDF5File.OpenMode.CREATE)) {
            IOUtils.canReadFile(outputHDF5File.getFile());

            outputHDF5File.makeStringArray("/data/annotation_names", sortedAnnotationKeys.toArray(new String[0]));
            HDF5Utils.writeChunkedDoubleMatrix(outputHDF5File, "/data/annotations", data.stream().map(vd -> vd.annotations).toArray(double[][]::new), maximumChunkSize);
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
