package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.HDF5Utils;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

/**
 * TODO
 */
@CommandLineProgramProperties(
        // TODO
        summary = "",
        oneLineSummary = "",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
public class LabeledVariantAnnotationsWalker extends MultiVariantWalker {

    private static final String ANNOTATIONS_HDF5_SUFFIX = ".annot.hdf5";
    private static final String DEFAULT_VCF_SUFFIX = ".vcf";

    private static final int DEFAULT_CHUNK_DIVISOR = 16;
    static final int DEFAULT_MAXIMUM_CHUNK_SIZE = HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX / DEFAULT_CHUNK_DIVISOR;

    private static final String SCORE_KEY = GATKVCFConstants.VQS_LOD_KEY;
    private static final String DUMMY_ALLELE = "<VQSR>";

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output prefix.")
    String outputPrefix;

    /**
     * Any set of VCF files to use as lists of training or truth sites.
     * Training - The program builds the model using input variants that overlap with these training sites.
     * Truth - The program uses these truth sites to determine where to set the cutoff in sensitivity.
     */
    @Argument(
            fullName = StandardArgumentDefinitions.RESOURCE_LONG_NAME,
            doc = "") // TODO
    private List<FeatureInput<VariantContext>> resources = new ArrayList<>(10);

    @Argument(
            fullName = "mode",
            shortName = "mode",
            doc = "Variant type to extract")
    public VariantTypeMode mode = VariantTypeMode.SNP;

    /**
     * Extract per-allele annotations.
     * Annotations should be specified using their full names with AS_ prefix.
     * Non-allele-specific (scalar) annotations will be applied to all alleles.
     */
    @Argument(
            fullName = "use-allele-specific-annotations",
            shortName = "AS",
            doc = "If specified, attempt to use the allele-specific versions of the specified annotations.",
            optional = true)
    public boolean useASAnnotations = false;

    /**
     * See the input VCF file's INFO field for a list of all available annotations.
     */
    @Argument(
            fullName = "use-annotation",
            shortName = "an",
            doc = "The names of the annotations to extract.",
            minElements = 1)
    private List<String> useAnnotations = new ArrayList<>();

    @Argument(
            fullName = "ignore-filter",
            doc = "If specified, use variants marked as filtered by the specified filter name in the input VCF file.",
            optional = true)
    private List<String> ignoreInputFilters = new ArrayList<>();

    @Argument(
            fullName = "ignore-all-filters",
            doc = "If specified, ignore all input filters.",
            optional = true)
    private boolean ignoreAllFilters = false;

    @Advanced
    @Argument(
            fullName = "trust-all-polymorphic",
            doc = "Trust that unfiltered records in the resources contain only polymorphic sites to decrease runtime.",
            optional = true)
    private boolean trustAllPolymorphic = false;

    @Advanced
    @Argument(
            doc = "Maximum HDF5 matrix chunk size.  Large matrices written to HDF5 are chunked into equally sized " +
                    "subsets of rows (plus a subset containing the remainder, if necessary) to avoid a hard limit in " +
                    "Java HDF5 on the number of elements in a matrix.  However, since a single row is not allowed to " +
                    "be split across multiple chunks, the number of columns must be less than the maximum number of " +
                    "values in each chunk.  Decreasing this number will reduce heap usage when writing chunks.",
            fullName = "maximum-chunk-size",
            minValue = 1,
            maxValue = HDF5Utils.MAX_NUMBER_OF_VALUES_PER_HDF5_MATRIX,
            optional = true
    )
    int maximumChunkSize = DEFAULT_MAXIMUM_CHUNK_SIZE;

    LabeledVariantAnnotationsData data;
    File outputAnnotationsFile;
    private File outputVCFFile;
    private final Set<String> ignoreInputFilterSet = new TreeSet<>();

    // TODO document
    public boolean isExtractAll() {
        return true;
    }

    public String getVCFSuffix() {
        return DEFAULT_VCF_SUFFIX;
    }

    public void beforeOnTraversalStart() {
        // override
    }

    @Override
    public void onTraversalStart() {

        beforeOnTraversalStart();

        final String vcfSuffix = getVCFSuffix();

        outputAnnotationsFile = new File(outputPrefix + ANNOTATIONS_HDF5_SUFFIX);
        outputVCFFile = new File(outputPrefix + vcfSuffix);

        for (final File outputFile : Arrays.asList(outputAnnotationsFile, outputVCFFile)) {
            if ((outputFile.exists() && !outputFile.canWrite()) ||
                    (!outputFile.exists() && !outputFile.getAbsoluteFile().getParentFile().canWrite())) {
                throw new UserException(String.format("Cannot create output file at %s.", outputFile));
            }
        }

        if (ignoreInputFilters != null) {
            ignoreInputFilterSet.addAll(ignoreInputFilters);
        }

        final Set<String> resourceLabels = new TreeSet<>();
        for (final FeatureInput<VariantContext> resource : resources) {
            final TreeSet<String> trackResourceLabels = resource.getTagAttributes().entrySet().stream()
                    .filter(e -> e.getValue().equals("true"))
                    .map(Map.Entry::getKey)
                    .sorted()
                    .collect(Collectors.toCollection(TreeSet::new));
            resourceLabels.addAll(trackResourceLabels);
            logger.info( String.format("Found %s track: labels = %s", resource.getName(), trackResourceLabels));
        }
        resourceLabels.forEach(String::intern);

        if (!resourceLabels.contains(LabeledVariantAnnotationsData.TRAINING_LABEL)) {
            throw new CommandLineException(
                    "No training set found! Please provide sets of known polymorphic loci marked with the training=true feature input tag. " +
                            "For example, --resource:hapmap,training=true,truth=true hapmapFile.vcf");
        }

        if (!resourceLabels.contains(LabeledVariantAnnotationsData.TRUTH_LABEL)) {
            throw new CommandLineException(
                    "No truth set found! Please provide sets of known polymorphic loci marked with the truth=true feature input tag. " +
                            "For example, --resource:hapmap,training=true,truth=true hapmapFile.vcf");
        }

        data = new LabeledVariantAnnotationsData(useAnnotations, resourceLabels, useASAnnotations);
    }

    @Override
    public void apply(final VariantContext vc,
                      final ReadsContext readsContext,
                      final ReferenceContext ref,
                      final FeatureContext featureContext) {
        addVariantDatum(vc, featureContext);
    }

    private void addVariantDatum(final VariantContext vc,
                                 final FeatureContext featureContext) {
        if (vc == null) {
            return;
        }
        if (!(ignoreAllFilters || vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters()))) {
            return;
        }
        if (VariantTypeMode.checkVariationClass(vc, mode) && !useASAnnotations) {
            final Set<String> overlappingResourceLabels = findOverlappingResourceLabels(vc, null, null, featureContext);
            if (isExtractAll() || !overlappingResourceLabels.isEmpty()) {
                data.addDatum(vc, vc.getReference(), null, overlappingResourceLabels);
                if (!isExtractAll()) {
                    data.addAlternateAlleles(vc.getAlternateAlleles());
                }
            }
        } else if (useASAnnotations) {
            for (final Allele altAllele : vc.getAlternateAlleles()) {
                if (!GATKVCFConstants.isSpanningDeletion(altAllele) && VariantTypeMode.checkVariationClass(vc, altAllele, mode)) {
                    final Set<String> overlappingResourceLabels = findOverlappingResourceLabels(vc, vc.getReference(), altAllele, featureContext);
                    //note that this may not be the minimal representation for the ref and alt allele
                    if (isExtractAll() || !overlappingResourceLabels.isEmpty()) {
                        data.addDatum(vc, vc.getReference(), altAllele, overlappingResourceLabels);
                        if (!isExtractAll()) {
                            data.addAlternateAlleles(Collections.singletonList(altAllele));
                        }
                    }
                }
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        return null;
    }

    private Set<String> findOverlappingResourceLabels(final VariantContext vc,
                                                      final Allele refAllele,
                                                      final Allele altAllele,
                                                      final FeatureContext featureContext) {
        final Set<String> overlappingResourceLabels = new TreeSet<>();
        for (final FeatureInput<VariantContext> resource : resources) {
            final List<VariantContext> resourceVCs = featureContext.getValues(resource, featureContext.getInterval().getStart());
            for (final VariantContext resourceVC : resourceVCs) {
                if (useASAnnotations && !doAllelesMatch(refAllele, altAllele, resourceVC)) {
                    continue;
                }
                if (isValidVariant(vc, resourceVC, trustAllPolymorphic)) {
                    overlappingResourceLabels.addAll(resource.getTagAttributes().entrySet().stream()
                            .filter(e -> e.getValue().equals("true"))
                            .map(Map.Entry::getKey)
                            .collect(Collectors.toSet()));
                }
            }
        }
        return overlappingResourceLabels;
    }

    private static boolean isValidVariant(final VariantContext vc,
                                          final VariantContext resourceVC,
                                          final boolean trustAllPolymorphic) {
        return resourceVC != null && resourceVC.isNotFiltered() && resourceVC.isVariant() && VariantTypeMode.checkVariationClass(vc, resourceVC) &&
                (trustAllPolymorphic || !resourceVC.hasGenotypes() || resourceVC.isPolymorphicInSamples());
    }

    private static boolean doAllelesMatch(final Allele refAllele,
                                          final Allele altAllele,
                                          final VariantContext resourceVC) {
        if (altAllele == null) {
            return true;
        }
        try {
            return GATKVariantContextUtils.isAlleleInList(refAllele, altAllele, resourceVC.getReference(), resourceVC.getAlternateAlleles());
        } catch (final IllegalStateException e) {
            throw new IllegalStateException("Reference allele mismatch at position " + resourceVC.getContig() + ':' + resourceVC.getStart() + " : ", e);
        }
    }

    void writeAnnotationsHDF5() {
        data.writeAnnotationsHDF5(outputAnnotationsFile, maximumChunkSize);
    }

    void writeVCF(final boolean writeAlleles,
                  final boolean writeTrainingOnly,
                  final boolean writeScores) {
        final VariantContextWriter vcfWriter = createVCFWriter(outputVCFFile);
        vcfWriter.writeHeader(constructVCFHeader());

        final List<LabeledVariantAnnotationsDatum> data = this.data.getData();

        // we need to sort in coordinate order in order to produce a valid VCF
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        data.sort((vd1, vd2) -> IntervalUtils.compareLocatables(vd1.loc, vd2.loc, sequenceDictionary));

        // create dummy alleles to be used
        List<Allele> alleles = Arrays.asList(Allele.create("N", true), Allele.create(DUMMY_ALLELE, false));

        for (int i = 0; i < data.size(); i++) {
            final LabeledVariantAnnotationsDatum datum = data.get(i);
            if (writeTrainingOnly && !datum.labels.contains(LabeledVariantAnnotationsData.TRAINING_LABEL)) {
                continue;
            }
            if (useASAnnotations) {
                alleles = Arrays.asList(datum.refAllele, datum.altAllele); //use the alleles to distinguish between multiallelics in AS mode
            } else if (writeAlleles) {
                final List<Allele> allelesToWrite = this.data.alternateAlleles.get(i);
                allelesToWrite.add(0, datum.refAllele);
                alleles = allelesToWrite;
            }
            final VariantContextBuilder builder = new VariantContextBuilder(SCORE_KEY, datum.loc.getContig(), datum.loc.getStart(), datum.loc.getEnd(), alleles);
            builder.attribute(VCFConstants.END_KEY, datum.loc.getEnd());

            if (writeScores) {
                builder.attribute(SCORE_KEY, String.format("%.4f", datum.score));
            }

            if (datum.labels.contains(LabeledVariantAnnotationsData.TRAINING_LABEL)) {
                builder.attribute(GATKVCFConstants.POSITIVE_LABEL_KEY, true);
            }

            vcfWriter.add(builder.make());
        }
        vcfWriter.close();
        logger.info(String.format("Recalibration VCF written to %s.", outputVCFFile.getAbsolutePath()));
    }

    private VCFHeader constructVCFHeader() {
        //TODO: this should be refactored/consolidated as part of
        // https://github.com/broadinstitute/gatk/issues/2112
        // https://github.com/broadinstitute/gatk/issues/121,
        // https://github.com/broadinstitute/gatk/issues/1116 and
        // Initialize VCF header lines
        Set<VCFHeaderLine> hInfo = getDefaultToolVCFHeaderLines();
        hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(SCORE_KEY));
        hInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.POSITIVE_LABEL_KEY));
        hInfo.add(GATKVCFHeaderLines.getFilterLine(VCFConstants.PASSES_FILTERS_v4));
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        if (hasReference()) {
            hInfo = VcfUtils.updateHeaderContigLines(
                    hInfo, referenceArguments.getReferencePath(), sequenceDictionary, true);
        } else if (null != sequenceDictionary) {
            hInfo = VcfUtils.updateHeaderContigLines(hInfo, null, sequenceDictionary, true);
        }
        return new VCFHeader(hInfo);
    }
}