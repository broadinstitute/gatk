package org.broadinstitute.hellbender.tools.walkers.sv;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.ClusteringParameters;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterWalker;
import org.broadinstitute.hellbender.tools.sv.cluster.StratifiedClusteringTableParser;
import org.broadinstitute.hellbender.tools.sv.concordance.ClosestSVFinder;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceAnnotator;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceLinkage;
import org.broadinstitute.hellbender.tools.sv.concordance.StratifiedConcordanceEngine;
import org.broadinstitute.hellbender.tools.sv.stratify.OptionalSVStratificationEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import picard.vcf.GenotypeConcordance;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>This tool calculates SV genotype concordance between an "evaluation" VCF and a "truth" VCF. For each evaluation
 * variant, a single truth variant is matched based on the following order of criteria:</p>
 *
 * <ol>
 *     <li>Total breakend distance</li>
 *     <li>Min breakend distance (among the two sides)</li>
 *     <li>Genotype concordance</li>
 * </ol>
 *
 * after meeting minimum overlap criteria. Evaluation VCF variants that are sucessfully matched are annotated with
 * genotype concordance metrics, including allele frequency of the truth variant. Concordance metrics are computed
 * on the intersection of sample sets of the two VCFs, but all other annotations including variant truth status
 * and allele frequency use all records and samples available. See output header for descriptions
 * of the specific fields. For multi-allelic CNVs, only a copy state concordance metric is
 * annotated. Allele frequencies will be recalculated automatically if unavailable in the provided VCFs.
 *
 * Minimum matching criteria can be specified as in {@link SVCluster} (e.g. reciprocal overlap). These serve to
 * improve computational efficiency, but may be set to relaxed values that can be iterated on post hoc. Note that
 * this set of parameters includes minimum sample overlap (Jaccard index of carrier samples), which should generally
 * be set to 0 for concordance analysis.
 *
 * This tool also allows supports stratification of the SVs into groups with specified matching criteria including SV type,
 * size range, and interval overlap. Please see the {@link GroupedSVCluster} tool documentation for further details
 * on how to specify stratification groups. Stratification only affects the criteria applied to each "eval" SV. In
 * other words, "truth" variants are not stratified and can match to an "eval" SV from any stratification group.
 *
 * Note that unlike {@link GroupedSVCluster}, this tool allows any variant to
 * match more than one stratification group. If this occurs, groups will be prioritized by their ordering in the input
 * stratification table, with groups appearing first receiving higher priority. While all matching groups will
 * be listed in the STRAT INFO field, the variant ID pertaining to the highest-priority group will be populated in
 * the TRUTH_VID field (groups with no matching variant are ignored). It is therefore recommended that the groups with
 * the most specific clustering criteria be listed as higher priority.
 *
 * The "default" stratification group, with clustering parameters specified directly through the clustering program
 * arguments (e.g. --depth-breakend-window, --pesr-interval-overlap, etc.), is always present and given lowest priority.
 *
 * <ul>
 *     <li>GENOTYPE_CONCORDANCE</li>
 *     <li>CNV_CONCORDANCE</li>
 *     <li>NON_REF_GENOTYPE_CONCORDANCE</li>
 *     <li>HET_PPV</li>
 *     <li>HET_SENSITIVITY</li>
 *     <li>HOMVAR_PPV</li>
 *     <li>HOMVAR_SENSITIVITY</li>
 *     <li>VAR_PPV</li>
 *     <li>VAR_SENSITIVITY</li>
 *     <li>VAR_SPECIFICITY</li>
 *     <li>TRUTH_VARIANT_ID</li>
 *     <li>TRUTH_AC</li>
 *     <li>TRUTH_AF</li>
 *     <li>TRUTH_AN</li>
 *     <li>TRUTH_RECIPROCAL_OVERLAP</li>
 *     <li>TRUTH_SIZE_SIMILARITY</li>
 *     <li>TRUTH_DISTANCE_START</li>
 *     <li>TRUTH_DISTANCE_END</li>
 * </ul>
 *
 * Allele frequency related fields (AF, AC, AN, TRUTH_AF, TRUTH_AC, TRUTH_AN) are passed through if already assigned in
 * the input VCFs and otherwise recalculated.
 *
 * Output records determined to be "true positives" are annotated with the following INFO fields according to overlap
 * with the highest priority matching stratification group:
 *
 * This tool performs a final sorting step on the emitted records using {@link SortingCollection}, which may
 * inflate memory usage and degrade performance on very large VCFs. Performance may be improved by reducing the
 * {@link #maxRecordsInRam} parameter.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Evaluation VCF
 *     </li>
 *     <li>
 *         Truth VCF
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         The evaluation VCF annotated with genotype concordance metrics
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVConcordance \
 *       --sequence-dictionary ref.dict \
 *       --eval evaluation.vcf.gz \
 *       --truth truth.vcf.gz \
 *       -O output.vcf.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Annotates structural variant genotype concordance",
        oneLineSummary = "Annotates structural variant genotype concordance",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class SVConcordance extends AbstractConcordanceWalker {

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected GATKPath outputFile;

    /**
     * Expected format is tab-delimited and contains columns NAME, RECIPROCAL_OVERLAP, SIZE_SIMILARITY, BREAKEND_WINDOW,
     * SAMPLE_OVERLAP. First line must be a header with column names. Comment lines starting with
     * {@link TableUtils#COMMENT_PREFIX} are ignored.
     */
    @Argument(
            doc = "Configuration file (.tsv) containing the clustering parameters for each group",
            fullName = GroupedSVCluster.CLUSTERING_CONFIG_FILE_LONG_NAME,
            optional = true
    )
    public GATKPath strataClusteringConfigFile;

    @Argument(fullName = SVClusterWalker.MAX_RECORDS_IN_RAM_LONG_NAME,
            doc = "When writing VCF files that need to be sorted, this will specify the number of records stored in " +
                    "RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a " +
                    "VCF file, and increases the amount of RAM needed.",
            optional=true)
    public int maxRecordsInRam = 1000;

    @ArgumentCollection
    protected final SVClusterEngineArgumentsCollection defaultClusteringArgs = new SVClusterEngineArgumentsCollection();
    @ArgumentCollection
    private final OptionalSVStratificationEngineArgumentsCollection stratArgs = new OptionalSVStratificationEngineArgumentsCollection();

    protected StratifiedConcordanceEngine engine;
    protected SAMSequenceDictionary dictionary;
    protected VariantContextWriter writer;

    // Output records are written in the order of a stable sort by (header contig index, start). When the eval header
    // orders contigs as the traversal does, records are held in an in-memory reorder buffer and written as soon as no
    // pending record can precede them; otherwise they all go through an on-disk sorting collection.
    protected SortingCollection<VariantContext> sortingBuffer;
    private PriorityQueue<PendingOutput> reorderBuffer;
    private Map<String, Integer> outputContigIndex;
    private long numOutputRecords = 0;
    private String lastAddedContig;
    private int lastAddedStart;


    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> true;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        // Use master sequence dictionary i.e. hg38 .dict file since the "best" dictionary is grabbed
        // from the VCF, which is sometimes out of order
        dictionary = getMasterSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }

        // Check that vcfs are sorted the same
        SequenceDictionaryUtils.validateDictionaries("eval", getEvalHeader().getSequenceDictionary(),
                "truth", getTruthHeader().getSequenceDictionary(), false, true);
        writer = createVCFWriter(outputFile);
        final VCFHeader header = createHeader(getEvalHeader());
        writer.writeHeader(header);
        outputContigIndex = new HashMap<>();
        for (final VCFContigHeaderLine contigLine : header.getContigLines()) {
            outputContigIndex.put(contigLine.getID(), contigLine.getContigIndex());
        }
        if (!outputContigIndex.isEmpty() && outputContigOrderMatchesTraversal()) {
            reorderBuffer = new PriorityQueue<>();
        } else {
            sortingBuffer = SortingCollection.newInstance(
                    VariantContext.class,
                    new VCFRecordCodec(header, true),
                    header.getVCFRecordComparator(),
                    maxRecordsInRam,
                    tmpDir.toPath());
        }

        // Concordance computations should be done on common samples only
        final Set<String> commonSamples = new HashSet<>(Sets.intersection(
                new HashSet<>(getEvalHeader().getGenotypeSamples()),
                new HashSet<>(getTruthHeader().getGenotypeSamples())));
        // Output used to pass through the sorting collection's VCF text round trip, which re-decodes genotypes when the
        // header's samples are not in sorted order and so drops FORMAT keys that no sample has a value for. Without the
        // round trip, reproduce that by not emitting the one key the annotator can leave valueless for every sample.
        final boolean omitUndeterminedCopyNumberEquality = reorderBuffer != null && !header.samplesWereAlreadySorted();
        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator(commonSamples, omitUndeterminedCopyNumberEquality);

        // Load stratification groups
        if ((stratArgs.configFile == null) ^ (strataClusteringConfigFile == null)) {
            throw new UserException.BadInput("Both --" + OptionalSVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME
                    + " and --" + GroupedSVCluster.CLUSTERING_CONFIG_FILE_LONG_NAME + " must be used together, but only one was specified.");
        }
        final Map<String, ClosestSVFinder> clusterEngineMap = new HashMap<>();
        if (strataClusteringConfigFile != null) {
            try (final TableReader<StratifiedClusteringTableParser.StratumParameters> tableReader = TableUtils.reader(strataClusteringConfigFile.toPath(), StratifiedClusteringTableParser::tableParser)) {
                for (final StratifiedClusteringTableParser.StratumParameters parameters : tableReader) {
                    // Identical parameters for each linkage type
                    final ClusteringParameters pesrParams = ClusteringParameters.createPesrParameters(parameters.reciprocalOverlap(), parameters.sizeSimilarity(), parameters.breakendWindow(), parameters.sampleOverlap());
                    final ClusteringParameters mixedParams = ClusteringParameters.createMixedParameters(parameters.reciprocalOverlap(), parameters.sizeSimilarity(), parameters.breakendWindow(), parameters.sampleOverlap());
                    final ClusteringParameters depthParams = ClusteringParameters.createDepthParameters(parameters.reciprocalOverlap(), parameters.sizeSimilarity(), parameters.breakendWindow(), parameters.sampleOverlap());
                    final SVConcordanceLinkage linkage = new SVConcordanceLinkage(dictionary);
                    linkage.setDepthOnlyParams(depthParams);
                    linkage.setMixedParams(mixedParams);
                    linkage.setEvidenceParams(pesrParams);
                    final ClosestSVFinder engine = new ClosestSVFinder(linkage, collapser::annotate, false, dictionary);
                    clusterEngineMap.put(parameters.name(), engine);
                }
            } catch (final IOException e) {
                throw new GATKException("IO error while reading config table", e);
            }
        }
        final SVStratificationEngine stratEngine = SVStratify.loadStratificationConfig(stratArgs.configFile, stratArgs, dictionary);
        engine = new StratifiedConcordanceEngine(clusterEngineMap, stratEngine, stratArgs, defaultClusteringArgs, collapser, dictionary);
    }


    /**
     * True if every contig in the output header is indexed in the same relative order as the traversal (sequence
     * dictionary) order, so that records from contigs not yet traversed can never sort before records already seen.
     */
    private boolean outputContigOrderMatchesTraversal() {
        int lastIndex = -1;
        for (final SAMSequenceRecord sequence : dictionary.getSequences()) {
            final Integer index = outputContigIndex.get(sequence.getSequenceName());
            if (index != null) {
                if (index <= lastIndex) {
                    return false;
                }
                lastIndex = index;
            }
        }
        return true;
    }

    @Override
    public Object onTraversalSuccess() {
        bufferOutput(engine.flush(true));
        if (!engine.isEmpty()) {
            throw new GATKException("Concordance engine is not empty, but it should be");
        }
        if (reorderBuffer != null) {
            while (!reorderBuffer.isEmpty()) {
                writer.add(reorderBuffer.poll().variant);
            }
        } else {
            for (final VariantContext variant : sortingBuffer) {
                writer.add(variant);
            }
        }
        return super.onTraversalSuccess();
    }

    private void bufferOutput(final Collection<VariantContext> variants) {
        for (final VariantContext variant : variants) {
            if (reorderBuffer != null) {
                reorderBuffer.add(new PendingOutput(variant, outputContigIndex.get(variant.getContig()), numOutputRecords++));
            } else {
                sortingBuffer.add(variant);
            }
        }
    }

    /**
     * Writes buffered records that sort at or before every record still to come. Records still to come are eval
     * variants active in the engine and variants not yet added, which start at or after the last added variant.
     * A pending record that ties with a later one keeps its place because it was produced first, which matches the
     * stable sort of the sorting collection.
     */
    private void writeReadyOutput() {
        final List<SimpleInterval> bounds = engine.getActiveEvalStarts();
        bounds.add(new SimpleInterval(lastAddedContig, lastAddedStart, lastAddedStart));
        while (!reorderBuffer.isEmpty() && precedesAll(reorderBuffer.peek(), bounds)) {
            writer.add(reorderBuffer.poll().variant);
        }
    }

    private boolean precedesAll(final PendingOutput pending, final List<SimpleInterval> bounds) {
        for (final SimpleInterval bound : bounds) {
            final int contigCompare = Integer.compare(pending.contigIndex, outputContigIndex.get(bound.getContig()));
            if (contigCompare > 0 || (contigCompare == 0 && pending.variant.getStart() > bound.getStart())) {
                return false;
            }
        }
        return true;
    }

    @Override
    public void closeTool() {
        if (sortingBuffer != null) {
            sortingBuffer.cleanup();
        }
        if (writer != null) {
            writer.close();
        }
        super.closeTool();
    }

    @Override
    public void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {
        if (truthVersusEval.hasTruth()) {
            final VariantContext truth = truthVersusEval.getTruth();
            final SVCallRecord site = SVCallRecordUtils.create(new VariantContextBuilder(truth).noGenotypes().make(), dictionary);
            engine.addTruthVariant(new LazyTruthRecord(site, truth, dictionary));
            lastAddedContig = site.getContigA();
            lastAddedStart = site.getPositionA();
        }
        if (truthVersusEval.hasEval()) {
            final SVCallRecord record = SVCallRecordUtils.create(truthVersusEval.getEval(), dictionary);
            engine.addEvalVariant(record);
            lastAddedContig = record.getContigA();
            lastAddedStart = record.getPositionA();
        }
        bufferOutput(engine.flush(false));
        if (reorderBuffer != null) {
            writeReadyOutput();
        }
    }

    /**
     * Strips all FORMAT fields except for the genotype and copy state
     */
    private static Genotype stripTruthGenotype(final Genotype genotype) {
        final GenotypeBuilder builder = new GenotypeBuilder(genotype.getSampleName()).alleles(genotype.getAlleles());
        if (genotype.hasExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT)) {
            builder.attribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT, genotype.getExtendedAttribute(GATKSVVCFConstants.COPY_NUMBER_FORMAT));
        }
        if (genotype.hasExtendedAttribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT)) {
            builder.attribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT, genotype.getExtendedAttribute(GATKSVVCFConstants.EXPECTED_COPY_NUMBER_FORMAT));
        }
        return builder.make();
    }

    /**
     * Output record waiting in the reorder buffer, ordered by output contig index, start, then production order.
     */
    private static final class PendingOutput implements Comparable<PendingOutput> {
        final VariantContext variant;
        final int contigIndex;
        final long order;

        PendingOutput(final VariantContext variant, final int contigIndex, final long order) {
            this.variant = variant;
            this.contigIndex = contigIndex;
            this.order = order;
        }

        @Override
        public int compareTo(final PendingOutput other) {
            int result = Integer.compare(contigIndex, other.contigIndex);
            if (result == 0) {
                result = Integer.compare(variant.getStart(), other.variant.getStart());
            }
            return result != 0 ? result : Long.compare(order, other.order);
        }
    }

    /**
     * Truth record whose genotypes are decoded only when first needed. Only truth variants chosen as the closest match,
     * compared in a genotype tie-break, or tested for sample overlap need genotypes; decoding every truth variant's
     * genotypes, all FORMAT fields included, dominated runtime and memory on large cohorts. On first access the
     * genotypes are decoded from the source variant and stripped to GT and copy state, as before.
     */
    private static final class LazyTruthRecord extends SVCallRecord {
        private VariantContext source; // released once genotypes are decoded
        private GenotypesContext genotypes;

        LazyTruthRecord(final SVCallRecord site, final VariantContext source, final SAMSequenceDictionary dictionary) {
            super(site.getId(), site.getContigA(), site.getPositionA(), site.getStrandA(), site.getContigB(),
                    site.getPositionB(), site.getStrandB(), site.getType(), site.getComplexSubtype(),
                    site.getComplexEventIntervals(), site.getLength(), site.getEvidence(), site.getAlgorithms(),
                    site.getAlleles(), Collections.emptyList(), site.getAttributes(), site.getFilters(),
                    site.getLog10PError(), dictionary);
            this.source = source;
        }

        @Override
        public GenotypesContext getGenotypes() {
            if (genotypes == null) {
                final GenotypesContext stripped = GenotypesContext.create();
                stripped.addAll(source.getGenotypes().stream().map(SVConcordance::stripTruthGenotype).collect(Collectors.toList()));
                genotypes = stripped;
                source = null;
            }
            return genotypes;
        }
    }

    @Override
    protected boolean areVariantsAtSameLocusConcordant(final VariantContext truth, final VariantContext eval) {
        return true;
    }

    protected VCFHeader createHeader(final VCFHeader header) {
        header.addMetaDataLine(new VCFFormatHeaderLine(GenotypeConcordance.CONTINGENCY_STATE_TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "The genotype concordance contingency state"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.TRUTH_CN_EQUAL_FORMAT, 1, VCFHeaderLineType.Integer, "Truth CNV copy state is equal (1=True, 0=False)"));
        header.addMetaDataLine(Concordance.TRUTH_STATUS_HEADER_LINE);
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.GENOTYPE_CONCORDANCE_INFO, 1, VCFHeaderLineType.Float, "Genotype concordance"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.COPY_NUMBER_CONCORDANCE_INFO, 1, VCFHeaderLineType.Float, "CNV copy number concordance"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.NON_REF_GENOTYPE_CONCORDANCE_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype concordance"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HET_PPV_INFO, 1, VCFHeaderLineType.Float, "Heterozygous genotype positive predictive value"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HET_SENSITIVITY_INFO, 1, VCFHeaderLineType.Float, "Heterozygous genotype sensitivity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HOMVAR_PPV_INFO, 1, VCFHeaderLineType.Float, "Homozygous genotype positive predictive value"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.HOMVAR_SENSITIVITY_INFO, 1, VCFHeaderLineType.Float, "Homozygous genotype sensitivity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.VAR_PPV_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype positive predictive value"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.VAR_SENSITIVITY_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype sensitivity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.VAR_SPECIFICITY_INFO, 1, VCFHeaderLineType.Float, "Non-ref genotype specificity"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, 1, VCFHeaderLineType.String, "Matching truth set variant id"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Truth set allele count"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Truth set allele number"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Truth set allele frequency"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_RECIPROCAL_OVERLAP_INFO, 1, VCFHeaderLineType.Float, "Reciprocal overlap with truth variant"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_SIZE_SIMILARITY_INFO, 1, VCFHeaderLineType.Float, "Size similarity with truth variant"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_DISTANCE_START_INFO, 1, VCFHeaderLineType.Integer, "Start coordinate distance in bp to truth variant's start"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_DISTANCE_END_INFO, 1, VCFHeaderLineType.Integer, "End coordinate distance in bp to truth variant's end"));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        SVStratify.addStratifyMetadata(header);
        return header;
    }
}
