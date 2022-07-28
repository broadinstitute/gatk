package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.tools.walkers.vqsr.VariantRecalibrator;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsData;
import org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

/**
 * Extracts site-level variant annotations, labels, and other metadata from a VCF file to HDF5 files.
 *
 * <p>
 *     This tool is intended to be used as the first step in a variant-filtering workflow that supersedes the
 *     {@link VariantRecalibrator} workflow. This tool extracts site-level annotations, labels, and other relevant metadata
 *     from variant sites (or alleles, in allele-specific mode) that are or are not present in specified labeled
 *     resource VCFs (e.g., training or calibration VCFs). The former, present sites are considered labeled; each site
 *     can have multiple labels. The latter sites are considered unlabeled and can be randomly downsampled using
 *     reservoir sampling; extraction of these is optional. The outputs of the tool are HDF5 files containing the
 *     extracted data for labeled and (optional) unlabeled variant sets, as well as a sites-only indexed VCF containing
 *     the labeled variants.
 * </p>
 * 
 * <p>
 *     The extracted sets can be provided as input to the {@link TrainVariantAnnotationsModel} tool
 *     to produce an annotation-based model for scoring variant calls. This model can in turn be provided
 *     along with a VCF file to the {@link ScoreVariantAnnotations} tool, which assigns a score to each call
 *     (with a lower score indicating that a call is more likely to be an artifact and should perhaps be filtered).
 *     Each score can also be converted to a corresponding sensitivity to a calibration set, if the latter is available.
 * </p>
 *
 * <p>
 *     Note that annotations and metadata are collected in memory during traversal until they are written to HDF5 files
 *     upon completion of the traversal. Memory requirements thus roughly scale linearly with both the number of sites
 *     extracted and the number of annotations.
 * </p>
 *
 * <p>
 *     Note that HDF5 files may be viewed using <a href="https://support.hdfgroup.org/products/java/hdfview/">hdfview</a>
 *     or loaded in Python using <a href="http://www.pytables.org/">PyTables</a> or <a href="http://www.h5py.org/">h5py</a>.
 * </p>
 * 
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Input VCF file. Site-level annotations will be extracted from the contained variants (or alleles, 
 *         if the {@value USE_ALLELE_SPECIFIC_ANNOTATIONS_LONG_NAME} argument is specified).
 *     </li>
 *     <li>
 *         Annotations to extract.
 *     </li>
 *     <li>
 *         Variant types (i.e., SNP and/or INDEL) to extract. Logic for determining variant type was retained from
 *         {@link VariantRecalibrator}; see {@link VariantType}. Extracting SNPs and INDELs separately in two runs of
 *         this tool can be useful if one wishes to extract different sets of annotations for each variant type,
 *         for example.
 *     </li>
 *     <li>
 *         (Optional) Resource VCF file(s). Each resource should be tagged with a label, which will be assigned to
 *         extracted sites that are present in the resource. In typical use, the {@value LabeledVariantAnnotationsData#TRAINING_LABEL}
 *         and {@value LabeledVariantAnnotationsData#CALIBRATION_LABEL} labels should be used to tag at least one resource
 *         apiece. The resulting sets of sites will be used for model training and conversion of scores to
 *         calibration-set sensitivity, respectively; the trustworthiness of the respective resources should be
 *         taken into account accordingly. The {@value LabeledVariantAnnotationsData#SNP_LABEL} label is
 *         reserved by the tool, as it is used to label sites determined to be SNPs, and thus it cannot be used to tag
 *         provided resources.
 *     </li>
 *     <li>
 *         (Optional) Maximum number of unlabeled variants (or alleles) to randomly sample with reservoir downsampling.
 *         If nonzero, annotations will also be extracted from unlabeled sites (i.e., those that are not present
 *         in the labeled resources).
 *     </li>
 *     <li>
 *         Output prefix.
 *         This is used as the basename for output files.
 *     </li>
 * </ul>
 *
 * <h3>Outputs</h3>
 *
 * <ul>
 *     <li>
 *         (Optional) Labeled-annotations HDF5 file (.annot.hdf5). Annotation data and metadata for those sites that
 *         are present in labeled resources are stored in the following HDF5 directory structure:
 *
 *         <p>
 *           |--- alleles<br>
 *           |    |--- alt<br>
 *           |    |--- ref<br>
 *           |--- annotations<br>
 *           |    |--- chunk_0<br>
 *           |    |--- ...<br>
 *           |    |--- chunk_{num_chunks - 1}<br>
 *           |    |--- names<br>
 *           |    |--- num_chunks<br>
 *           |    |--- num_columns<br>
 *           |    |--- num_rows<br>
 *           |--- intervals<br>
 *           |    |--- indexed_contig_names<br>
 *           |    |--- transposed_index_start_end<br>
 *           |--- labels<br>
 *           |    |--- snp<br>
 *           |    |--- ... (e.g., training, calibration, etc.)<br>
 *           |    |--- ...<br>
 *         </p>
 *
 *         <p>
 *             Here, each chunk is a double matrix, with dimensions given by (number of sites in the chunk) x (number of annotations).
 *             See the methods {@link HDF5Utils#writeChunkedDoubleMatrix} and {@link HDF5Utils#writeIntervals} for additional details.
 *             If {@value USE_ALLELE_SPECIFIC_ANNOTATIONS_LONG_NAME} is specified, each record corresponds to an individual allele;
 *             otherwise, each record corresponds to a variant site, which may contain multiple alleles.
 *             Storage of alleles can be omitted using the {@value OMIT_ALLELES_IN_HDF5_LONG_NAME} argument, which will reduce
 *             the size of the file. This file will only be produced if resources are provided and the number of extracted
 *             labeled sites is nonzero.
 *         </p>
 *
 *     </li>
 *     <li>
 *         Labeled sites-only VCF file and index. The VCF will not be gzipped if the {@value DO_NOT_GZIP_VCF_OUTPUT_LONG_NAME}
 *         argument is set to true. The VCF can be provided as a resource in subsequent runs of
 *         {@link ScoreVariantAnnotations} and used to indicate labeled sites that were extracted.
 *         This can be useful if the {@value StandardArgumentDefinitions#INTERVALS_LONG_NAME} argument was used to
 *         subset sites in training or calibration resources for extraction; this may occur when setting up
 *         training/validation/test splits, for example. Note that records for the random sample of unlabeled sites are
 *         currently not included in the VCF.
 *     </li>
 *     <li>
 *         (Optional) Unlabeled-annotations HDF5 file. This will have the same directory structure as in the
 *         labeled-annotations HDF5 file. However, note that records are currently written in the order they
 *         appear in the downsampling reservoir after random sampling, and hence, are not in genomic order.
 *         This file will only be produced if a nonzero value of the {@value MAXIMUM_NUMBER_OF_UNLABELED_VARIANTS_LONG_NAME}
 *         argument is provided.
 *     </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <p>
 *     Extract annotations from training/calibration SNP/INDEL sites, producing the outputs
 *     1) {@code extract.annot.hdf5}, 2) {@code extract.vcf.gz}, and 3) {@code extract.vcf.gz.tbi}.
 *     The HDF5 file can then be provided to {@link TrainVariantAnnotationsModel}
 *     to train a model using a positive-only approach.
 *
 * <pre>
 *     gatk ExtractVariantAnnotations \
 *          -V input.vcf \
 *          -A annotation_1 \
 *          ...
 *          -A annotation_N \
 *          --mode SNP \
 *          --resource snp-training,training=true snp-training.vcf \
 *          --resource snp-calibration,calibration=true snp-calibration.vcf \
 *          --mode INDEL \
 *          --resource indel-training,training=true indel-training.vcf \
 *          --resource indel-calibration,calibration=true indel-calibration.vcf \
 *          -O extract
 * </pre>
 * </p>
 *
 * <p>
 *     Extract annotations from both training/calibration SNP/INDEL sites and a random sample of
 *     1000000 unlabeled (i.e., non-training/calibration) sites, producing the outputs
 *     1) {@code extract.annot.hdf5}, 2) {@code extract.unlabeled.annot.hdf5}, 3) {@code extract.vcf.gz},
 *     and 4) {@code extract.vcf.gz.tbi}. The HDF5 files can then be provided to {@link TrainVariantAnnotationsModel}
 *     to train a model using a positive-negative approach (similar to that used in {@link VariantRecalibrator}).
 *
 * <pre>
 *     gatk ExtractVariantAnnotations \
 *          -V input.vcf \
 *          -A annotation_1 \
 *          ...
 *          -A annotation_N \
 *          --mode SNP \
 *          --resource snp-training,training=true snp-training.vcf \
 *          --resource snp-calibration,calibration=true snp-calibration.vcf \
 *          --mode INDEL \
 *          --resource indel-training,training=true indel-training.vcf \
 *          --resource indel-calibration,calibration=true indel-calibration.vcf \
 *          --maximum-number-of-unlableled-variants 1000000
 *          -O extract
 * </pre>
 * </p>
 *
 * <p>
 *     In the (atypical) event that resource VCFs are unavailable, one can still extract annotations from a random sample of
 *     unlabeled sites, producing the outputs 1) {@code extract.unlabeled.annot.hdf5},
 *     2) {@code extract.vcf.gz} (which will contain no records), and 3) {@code extract.vcf.gz.tbi}.
 *     This random sample cannot be used by {@link TrainVariantAnnotationsModel}, but may still be useful for
 *     exploratory analyses.
 *
 * <pre>
 *     gatk ExtractVariantAnnotations \
 *          -V input.vcf \
 *          -A annotation_1 \
 *          ...
 *          -A annotation_N \
 *          --mode SNP \
 *          --mode INDEL \
 *          --maximum-number-of-unlableled-variants 1000000
 *          -O extract
 * </pre>
 * </p>
 *
 * DEVELOPER NOTE: See documentation in {@link LabeledVariantAnnotationsWalker}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Extracts site-level variant annotations, labels, and other metadata from a VCF file to HDF5 files.",
        oneLineSummary = "Extracts site-level variant annotations, labels, and other metadata from a VCF file to HDF5 files",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class ExtractVariantAnnotations extends LabeledVariantAnnotationsWalker {

    public static final String MAXIMUM_NUMBER_OF_UNLABELED_VARIANTS_LONG_NAME = "maximum-number-of-unlabeled-variants";
    public static final String RESERVOIR_SAMPLING_RANDOM_SEED_LONG_NAME = "reservoir-sampling-random-seed";

    public static final String UNLABELED_TAG = ".unlabeled";

    @Argument(
            fullName = MAXIMUM_NUMBER_OF_UNLABELED_VARIANTS_LONG_NAME,
            doc = "Maximum number of unlabeled variants to extract. " +
                    "If greater than zero, reservoir sampling will be used to randomly sample this number " +
                    "of sites from input sites that are not present in the specified resources.",
            minValue = 0)
    private int maximumNumberOfUnlabeledVariants = 0;

    @Argument(
            fullName = RESERVOIR_SAMPLING_RANDOM_SEED_LONG_NAME,
            doc = "Random seed to use for reservoir sampling of unlabeled variants.")
    private int reservoirSamplingRandomSeed = 0;

    private RandomGenerator rng;
    private LabeledVariantAnnotationsData unlabeledDataReservoir; // will not be sorted in genomic order
    private int unlabeledIndex = 0;

    @Override
    public void afterOnTraversalStart() {
        if (!resourceLabels.contains(LabeledVariantAnnotationsData.TRAINING_LABEL)) {
            logger.warn("No training set found! If you are using the downstream TrainVariantAnnotationsModel and ScoreVariantAnnotations tools, " +
                    "provide sets of known polymorphic loci marked with the training=true feature input tag. " +
                    "For example, --resource:hapmap,training=true hapmap.vcf");
        }
        if (!resourceLabels.contains(LabeledVariantAnnotationsData.CALIBRATION_LABEL)) {
            logger.warn("No calibration set found! If you are using the downstream TrainVariantAnnotationsModel and ScoreVariantAnnotations tools " +
                    "and wish to convert scores to sensitivity to a calibration set of variants, " +
                    "provide sets of known polymorphic loci marked with the calibration=true feature input tag. " +
                    "For example, --resource:hapmap,calibration=true hapmap.vcf");
        }

        rng = RandomGeneratorFactory.createRandomGenerator(new Random(reservoirSamplingRandomSeed));
        unlabeledDataReservoir = maximumNumberOfUnlabeledVariants == 0
                ? null
                : new LabeledVariantAnnotationsData(annotationNames, resourceLabels, useASAnnotations, maximumNumberOfUnlabeledVariants);
    }

    @Override
    protected void nthPassApply(final VariantContext variant,
                                final ReadsContext readsContext,
                                final ReferenceContext referenceContext,
                                final FeatureContext featureContext,
                                final int n) {
        if (n == 0) {
            final List<Triple<List<Allele>, VariantType, TreeSet<String>>> metadata = extractVariantMetadata(
                    variant, featureContext, unlabeledDataReservoir != null);
            final boolean isVariantExtracted = !metadata.isEmpty();
            if (isVariantExtracted) {
                final boolean isUnlabeled = metadata.stream().map(Triple::getRight).allMatch(Set::isEmpty);
                if (!isUnlabeled) {
                    addExtractedVariantToData(data, variant, metadata);
                    writeExtractedVariantToVCF(variant, metadata);
                } else {
                    // Algorithm R for reservoir sampling: https://en.wikipedia.org/wiki/Reservoir_sampling#Simple_algorithm
                    if (unlabeledIndex < maximumNumberOfUnlabeledVariants) {
                        addExtractedVariantToData(unlabeledDataReservoir, variant, metadata);
                    } else {
                        final int j = rng.nextInt(unlabeledIndex);
                        if (j < maximumNumberOfUnlabeledVariants) {
                            setExtractedVariantInData(unlabeledDataReservoir, variant, metadata, j);
                        }
                    }
                    unlabeledIndex++;
                }
            }
        }
    }

    @Override
    protected void afterNthPass(final int n) {
        if (n == 0) {
            writeAnnotationsToHDF5();
            data.clear();
            if (unlabeledDataReservoir != null) {
                writeUnlabeledAnnotationsToHDF5();
                // TODO write extracted unlabeled variants to VCF, which can be used to mark extraction in scoring step
                unlabeledDataReservoir.clear();
            }
            if (vcfWriter != null) {
                vcfWriter.close();
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {

        logger.info(String.format("%s complete.", getClass().getSimpleName()));

        return null;
    }

    private static void setExtractedVariantInData(final LabeledVariantAnnotationsData data,
                                                  final VariantContext variant,
                                                  final List<Triple<List<Allele>, VariantType, TreeSet<String>>> metadata,
                                                  final int index) {
        data.set(index, variant,
                metadata.stream().map(Triple::getLeft).collect(Collectors.toList()),
                metadata.stream().map(Triple::getMiddle).collect(Collectors.toList()),
                metadata.stream().map(Triple::getRight).collect(Collectors.toList()));
    }

    private void writeUnlabeledAnnotationsToHDF5() {
        final File outputUnlabeledAnnotationsFile = new File(outputPrefix + UNLABELED_TAG + ANNOTATIONS_HDF5_SUFFIX);
        if (unlabeledDataReservoir.size() == 0) {
            throw new GATKException("No unlabeled variants were present in the input VCF.");
        }
        for (final VariantType variantType : variantTypesToExtract) {
            logger.info(String.format("Extracted unlabeled annotations for %d variants of type %s.",
                    unlabeledDataReservoir.getVariantTypeFlat().stream().mapToInt(t -> t == variantType ? 1 : 0).sum(), variantType));
        }
        logger.info(String.format("Extracted unlabeled annotations for %s total variants.", unlabeledDataReservoir.size()));

        logger.info("Writing unlabeled annotations...");
        // TODO coordinate sort
        unlabeledDataReservoir.writeHDF5(outputUnlabeledAnnotationsFile, omitAllelesInHDF5);
        logger.info(String.format("Unlabeled annotations and metadata written to %s.", outputUnlabeledAnnotationsFile.getAbsolutePath()));
    }
}