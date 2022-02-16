package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.LinkedListMultimap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5LibException;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
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
public final class LabeledVariantAnnotationsData {
    private static final Logger logger = LogManager.getLogger(LabeledVariantAnnotationsData.class);

    // chunk size in temporary annotation files
    private static final int CHUNK_DIVISOR = 16;
    private static final int MAXIMUM_CHUNK_SIZE = HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX / CHUNK_DIVISOR;

    // TODO make labels enum?
    public static final String TRAINING_LABEL = "training";
    public static final String TRUTH_LABEL = "truth";

    private final List<String> sortedAnnotationNames;
    final List<String> sortedLabels;

    private final LinkedListMultimap<VariantContext, LabeledVariantAnnotationsDatum> variantContextsToDataMultimap;
    private final boolean useASAnnotations;

    public LabeledVariantAnnotationsData(final Collection<String> annotationNames,
                                         final Collection<String> labels,
                                         final int initialSize,
                                         final boolean useASAnnotations) {
        this.variantContextsToDataMultimap = LinkedListMultimap.create(initialSize);
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

    public int size() {
        return variantContextsToDataMultimap.size();
    }

    public void clear() {
        variantContextsToDataMultimap.clear();
    }

    public void add(final VariantContext vc,
                    final List<Allele> altAllelesPerDatum,
                    final List<VariantType> variantTypePerDatum,
                    final List<Set<String>> labelsPerDatum) {
        if (!useASAnnotations) {
            // if altAllele is null we are not in AS mode
            variantContextsToDataMultimap.put(vc,
                    new LabeledVariantAnnotationsDatum(vc,
                            altAllelesPerDatum,
                            variantTypePerDatum.get(0),
                            labelsPerDatum.get(0),
                            sortedAnnotationNames,
                            useASAnnotations));
        } else {
            // AS mode
            variantContextsToDataMultimap.putAll(vc,
                    IntStream.range(0, altAllelesPerDatum.size()).boxed()
                            .map(i -> new LabeledVariantAnnotationsDatum(vc,
                                    Collections.singletonList(altAllelesPerDatum.get(i)),
                                    variantTypePerDatum.get(i),
                                    labelsPerDatum.get(i),
                                    sortedAnnotationNames,
                                    useASAnnotations))
                            .collect(Collectors.toList()));
        }
    }

    private List<LabeledVariantAnnotationsDatum> getFlattenedData() {
        return Collections.unmodifiableList(variantContextsToDataMultimap.values());
    }

    public void writeBatchToHDF5(final File outputFile,
                                 final int batchIndex,
                                 final boolean omitAllelesInHDF5) {
        // TODO validate
        try (final HDF5File outputHDF5File = new HDF5File(outputFile, batchIndex == 0 ? HDF5File.OpenMode.CREATE : HDF5File.OpenMode.READ_WRITE)) {
            IOUtils.canReadFile(outputHDF5File.getFile());

            HDF5Utils.writeIntervals(outputHDF5File, String.format("/intervals/chunk_%d", batchIndex),
                    getFlattenedData().stream().map(datum -> datum.interval).collect(Collectors.toList()));
            if (!omitAllelesInHDF5) {
                outputHDF5File.makeStringArray(String.format("/ref/chunk_%d", batchIndex),
                        getFlattenedData().stream().map(datum -> datum.refAllele.getDisplayString()).toArray(String[]::new));
                if (!useASAnnotations) {
                    outputHDF5File.makeStringArray(String.format("/alt/chunk_%d", batchIndex),
                            getFlattenedData().stream()
                                    .map(datum -> datum.altAlleles.stream().map(Allele::getDisplayString).collect(Collectors.joining(",")))
                                    .toArray(String[]::new));
                } else {
                    outputHDF5File.makeStringArray(String.format("/alt/chunk_%d", batchIndex),
                            getFlattenedData().stream().map(datum -> datum.altAlleles.get(0).getDisplayString()).toArray(String[]::new));
                }
            }
            outputHDF5File.makeStringArray("/annotations/names", sortedAnnotationNames.toArray(new String[0]));

            // TODO clean this up
            // for the below quantities, we also write the additional fields num_chunks, num_rows, and num_columns
            // so that HDF5Utils.readChunkedDoubleMatrix can be used to read the result
            final double[][] annotations = getFlattenedData().stream().map(datum -> datum.annotations).toArray(double[][]::new);
            outputHDF5File.makeDoubleMatrix(String.format("/annotations/chunk_%d", batchIndex), annotations);
            outputHDF5File.makeDouble("/annotations" + HDF5Utils.NUMBER_OF_ROWS_SUB_PATH,
                    batchIndex == 0 ? annotations.length : (int) outputHDF5File.readDouble("/annotations" + HDF5Utils.NUMBER_OF_ROWS_SUB_PATH) + annotations.length);
            outputHDF5File.makeDouble("/annotations" + HDF5Utils.NUMBER_OF_COLUMNS_SUB_PATH, annotations[0].length);
            outputHDF5File.makeDouble("/annotations" + HDF5Utils.NUMBER_OF_CHUNKS_SUB_PATH, batchIndex + 1);

            final double[] isSNP = getFlattenedData().stream().mapToDouble(datum -> datum.variantType == VariantType.SNP ? 1 : 0).toArray();
            outputHDF5File.makeDoubleArray(String.format("/labels/snp/chunk_%d", batchIndex), isSNP);
            outputHDF5File.makeDouble("/labels/snp" + HDF5Utils.NUMBER_OF_ROWS_SUB_PATH,
                    batchIndex == 0 ? isSNP.length : (int) outputHDF5File.readDouble("/labels/snp" + HDF5Utils.NUMBER_OF_ROWS_SUB_PATH) + isSNP.length);
            outputHDF5File.makeDouble("/labels/snp" + HDF5Utils.NUMBER_OF_COLUMNS_SUB_PATH, 1);
            outputHDF5File.makeDouble("/labels/snp" + HDF5Utils.NUMBER_OF_CHUNKS_SUB_PATH, batchIndex + 1);

            for (final String label : sortedLabels) {
                final double[] isLabel = getFlattenedData().stream().mapToDouble(datum -> datum.labels.contains(label) ? 1 : 0).toArray();
                outputHDF5File.makeDoubleArray(String.format("/labels/%s/chunk_%d", label, batchIndex), isLabel);
                outputHDF5File.makeDouble(String.format("/labels/%s/%s", label, HDF5Utils.NUMBER_OF_ROWS_SUB_PATH),
                        batchIndex == 0 ? isSNP.length : (int) outputHDF5File.readDouble(String.format("/labels/%s/%s", label, HDF5Utils.NUMBER_OF_ROWS_SUB_PATH)) + isSNP.length);
                outputHDF5File.makeDouble(String.format("/labels/%s/%s", label, HDF5Utils.NUMBER_OF_COLUMNS_SUB_PATH), 1);
                outputHDF5File.makeDouble(String.format("/labels/%s/%s", label, HDF5Utils.NUMBER_OF_CHUNKS_SUB_PATH), batchIndex + 1);
            }
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during writing of annotations and metadata to (chunk_%d) (%s). Output file at %s may be in a bad state.",
                    batchIndex, exception, outputFile.getAbsolutePath()));
        }
        logger.debug(String.format("Annotations and metadata (chunk_%d) written to %s.", batchIndex, outputFile.getAbsolutePath()));
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

    public static double[][] readAnnotationsBatch(final File annotationsFile,
                                                  final int batchIndex) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return annotationsHDF5File.readDoubleMatrix(String.format("/annotations/chunk_%d", batchIndex));
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of annotations (chunk_%d) from %s: %s",
                    batchIndex, annotationsFile.getAbsolutePath(), exception));
        }
    }

    public static List<Boolean> readLabel(final File annotationsFile,
                                          final String label) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return Arrays.stream(readChunkedDoubleArray(annotationsHDF5File, String.format("/labels/%s", label))).boxed().map(d -> d == 1).collect(Collectors.toList());
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of label %s from %s: %s",
                    label, annotationsFile.getAbsolutePath(), exception));
        }
    }

    public static List<Boolean> readLabelBatch(final File annotationsFile,
                                               final String label,
                                               final int batchIndex) {
        try (final HDF5File annotationsHDF5File = new HDF5File(annotationsFile, HDF5File.OpenMode.READ_ONLY)) {
            IOUtils.canReadFile(annotationsHDF5File.getFile());
            return Arrays.stream(annotationsHDF5File.readDoubleArray(String.format("/labels/%s/chunk_%d", label, batchIndex))).boxed().map(d -> d == 1).collect(Collectors.toList());
        } catch (final HDF5LibException exception) {
            throw new GATKException(String.format("Exception encountered during reading of label %s from %s: %s",
                    label, annotationsFile.getAbsolutePath(), exception));
        }
    }

    public static double[] readChunkedDoubleArray(final HDF5File file,
                                                  final String path) {
        Utils.nonNull(file);
        IOUtils.canReadFile(file.getFile());
        Utils.nonNull(path);

        final String numRowsPath = path + HDF5Utils.NUMBER_OF_ROWS_SUB_PATH;
        final String numColumnsPath = path + HDF5Utils.NUMBER_OF_COLUMNS_SUB_PATH;
        final String numChunksPath = path + HDF5Utils.NUMBER_OF_CHUNKS_SUB_PATH;
        Utils.validateArg(file.isPresent(numRowsPath) && file.isPresent(numColumnsPath) && file.isPresent(numChunksPath),
                String.format("HDF5 file %s does not contain a chunked array in path %s.", file.getFile().getAbsolutePath(), path));

        final int numRows = (int) file.readDouble(numRowsPath);
        final int numColumns = (int) file.readDouble(numColumnsPath);
        final int numChunks = (int) file.readDouble(numChunksPath);

        Utils.validateArg(numColumns == 1,
                String.format("HDF5 file %s does not contain a chunked array with a single column in path %s.", file.getFile().getAbsolutePath(), path));

        final double[] fullArray = new double[numRows];
        int numRowsRead = 0;
        for (int chunkIndex = 0; chunkIndex < numChunks; chunkIndex++) {
            final double[] arrayChunk = file.readDoubleArray(path + HDF5Utils.CHUNK_INDEX_PATH_SUFFIX + chunkIndex);
            if (numRowsRead + arrayChunk.length > numRows) {
                throw new UserException.BadInput("Array chunk contains too many rows.");
            }
            System.arraycopy(arrayChunk, 0, fullArray, numRowsRead, arrayChunk.length);
            numRowsRead += arrayChunk.length;
        }
        if (numRowsRead != numRows) {
            throw new UserException.BadInput("Array chunks do not contain expected total number of rows.");
        }
        return fullArray;
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

//    void writeBatchToVCF(final VariantContextWriter vcfWriter,
//                         final boolean writeAlleles,
//                         final boolean writeTrainingOnly,
//                         final boolean writeScores) {
//        // TODO validate
//
////        TODO shouldn't this already be sorted?
//        // we need to sort in coordinate order in order to produce a valid VCF
////        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
////        dataBatch.sort((vd1, vd2) -> IntervalUtils.compareLocatables(vd1.loc, vd2.loc, sequenceDictionary));
//
//        List<Allele> alleles;
//
//        for (int i = 0; i < size(); i++) {
//            final AlleleLabeledAnnotationsDatum datum = dataBatch.get(i);
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
