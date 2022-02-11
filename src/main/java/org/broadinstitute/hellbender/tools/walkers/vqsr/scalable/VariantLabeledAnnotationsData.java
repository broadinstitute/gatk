package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/*
 * TODO this whole class needs refactoring. it's cleaned up significantly from VQSR version,
 *  but still has a long way to go. we should decide how strongly to couple to the tranche code and refactor
 *  that at the same time
 */
final class VariantLabeledAnnotationsData {
    private static final Logger logger = LogManager.getLogger(VariantLabeledAnnotationsData.class);

    private static final int INITIAL_SIZE = 1000000;

    // TODO make labels enum?
    static final String TRAINING_LABEL = "training";
    static final String TRUTH_LABEL = "truth";

    private final List<String> sortedAnnotationKeys;
    final List<String> sortedLabels;

    private final List<VariantContext> variantContexts;
    private final List<List<AlleleLabeledAnnotationsDatum>> data;

    VariantLabeledAnnotationsData(final Collection<String> annotationKeys,
                                  final Collection<String> labels,
                                  final int initialSize) {
        this.variantContexts = new ArrayList<>(initialSize);
        this.data = new ArrayList<>(initialSize);
        this.sortedAnnotationKeys = ImmutableList.copyOf(annotationKeys.stream().distinct().sorted().collect(Collectors.toList()));
        if (sortedAnnotationKeys.size() != annotationKeys.size()) {
            logger.warn(String.format("Ignoring duplicate annotations: %s.", Utils.getDuplicatedItems(annotationKeys)));
        }
        this.sortedLabels = ImmutableList.copyOf(labels.stream().distinct().sorted().collect(Collectors.toList()));
        if (sortedLabels.size() != labels.size()) {
            logger.warn(String.format("Ignoring duplicate labels: %s.", Utils.getDuplicatedItems(labels)));
        }
    }

    public List<VariantContext> getVariantContexts() {
        return Collections.unmodifiableList(variantContexts);
    }

    public List<List<AlleleLabeledAnnotationsDatum>> getData() {
        return Collections.unmodifiableList(data);
    }

    VariantLabeledAnnotationsData(final Collection<String> annotationKeys,
                                  final Collection<String> labels) {
        this(annotationKeys, labels, INITIAL_SIZE);
    }

    int size() {
        return data.size();
    }

    void clear() {
        variantContexts.clear();
        data.clear();
    }

    void add(final VariantContext vc,
             final List<Allele> altAlleles,
             final List<VariantType> variantTypesPerAltAllele,
             final List<Set<String>> labelsPerAltAllele) {
        variantContexts.add(vc);
        if (altAlleles == null) {
            // if altAllele is null we are not in AS mode
            data.add(Collections.singletonList(new AlleleLabeledAnnotationsDatum(vc,
                    null, variantTypesPerAltAllele.get(0), labelsPerAltAllele.get(0), sortedAnnotationKeys)));
        } else {
            // AS mode
            data.add(IntStream.range(0, altAlleles.size()).boxed()
                    .map(i -> new AlleleLabeledAnnotationsDatum(vc,
                            altAlleles.get(i), variantTypesPerAltAllele.get(i), labelsPerAltAllele.get(i), sortedAnnotationKeys))
                    .collect(Collectors.toList()));
        }
    }

    List<String> getSortedAnnotationKeys() {
        return sortedAnnotationKeys;
    }

    void writeLabeledAnnotationsBatchToHDF5(final File outputFile,
                                            final int batchIndex) {
        // TODO validate
        try (final HDF5File outputHDF5File = new HDF5File(outputFile, HDF5File.OpenMode.CREATE)) {
            IOUtils.canReadFile(outputHDF5File.getFile());

            outputHDF5File.makeStringArray("/data/annotation_names", sortedAnnotationKeys.toArray(new String[0]));
            outputHDF5File.makeDoubleMatrix(String.format("/data/annotations/chunk_%d", batchIndex),
                    data.stream().flatMap(List::stream).map(vd -> vd.annotations).toArray(double[][]::new));
            for (final String label : sortedLabels) {
                outputHDF5File.makeDoubleArray(String.format("/data/is_%s/chunk_%d", label, batchIndex),
                        data.stream().flatMap(List::stream).mapToDouble(vd -> vd.labels.contains(label) ? 1 : 0).toArray());
            }
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of labels and annotations to (chunk_%d) (%s). Output file at %s may be in a bad state.",
                    batchIndex, exception, outputFile.getAbsolutePath()));
        }
        logger.debug(String.format("Labels and annotations (chunk_%d) written to %s.", batchIndex, outputFile.getAbsolutePath()));
    }

    static double[][] readAnnotationsBatch(final File annotationsFile,
                                           final int batchIndex) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return HDF5Utils.readChunkedDoubleMatrix(annotationsHDF5File, String.format("/data/annotations/chunk_%d", batchIndex));
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of annotations (chunk_%d) from %s: %s",
                    batchIndex, annotationsFile.getAbsolutePath(), exception));
        }
    }

    static List<Boolean> readLabelBatch(final File annotationsFile,
                                        final String label,
                                        final int batchIndex) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return readBooleanList(annotationsHDF5File, String.format("/data/is_%s/chunk_%d", label, batchIndex));
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of label %s (chunk_%d) from %s: %s",
                    label, batchIndex, annotationsFile.getAbsolutePath(), exception));
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

//    void writeBatchToVCF(final VariantContextWriter vcfWriter,
//                         final boolean writeAlleles,
//                         final boolean writeTrainingOnly,
//                         final boolean writeScores) {
//
//
//        final List<AlleleLabeledAnnotationsDatum> data = this.dataBatch.getFlattenedData();
//
//        // we need to sort in coordinate order in order to produce a valid VCF
//        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
//        data.sort((vd1, vd2) -> IntervalUtils.compareLocatables(vd1.loc, vd2.loc, sequenceDictionary));
//
//        // create dummy alleles to be used
//        List<Allele> alleles = Arrays.asList(Allele.create("N", true), Allele.create(DUMMY_ALLELE, false));
//
//        for (int i = 0; i < data.size(); i++) {
//            final AlleleLabeledAnnotationsDatum datum = data.get(i);
//            if (writeTrainingOnly && !datum.labels.contains(VariantLabeledAnnotationsData.TRAINING_LABEL)) {
//                continue;
//            }
//            if (useASAnnotations) {
//                alleles = Arrays.asList(datum.refAllele, datum.altAllele); //use the alleles to distinguish between multiallelics in AS mode
//            } else if (writeAlleles) {
//                final List<Allele> allelesToWrite = this.dataBatch.alternateAlleles.get(i);
//                allelesToWrite.add(0, datum.refAllele);
//                alleles = allelesToWrite;
//            }
//            final VariantContextBuilder builder = new VariantContextBuilder(SCORE_KEY, datum.loc.getContig(), datum.loc.getStart(), datum.loc.getEnd(), alleles);
//            builder.attribute(VCFConstants.END_KEY, datum.loc.getEnd());
//
//            if (writeScores) {
//                builder.attribute(SCORE_KEY, String.format("%.4f", datum.score));
//            }
//
//            if (datum.labels.contains(VariantLabeledAnnotationsData.TRAINING_LABEL)) {
//                builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
//            }
//
//            vcfWriter.add(builder.make());
//        }
//        vcfWriter.close();
//        logger.info(String.format("Recalibration VCF written to %s.", outputVCFFile.getAbsolutePath()));
//    }
}
