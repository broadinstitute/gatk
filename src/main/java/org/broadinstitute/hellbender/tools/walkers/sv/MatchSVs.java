package org.broadinstitute.hellbender.tools.walkers.sv;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.collections4.Predicate;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
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
import org.broadinstitute.hellbender.tools.sv.concordance.*;
import org.broadinstitute.hellbender.tools.sv.stratify.OptionalSVStratificationEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngine;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.IOException;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Annotates site-level structural variant matches across VCFs",
        oneLineSummary = "Annotates site-level structural variant matches across VCFs",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)

public final class MatchSVs extends AbstractConcordanceWalker {
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
    protected SortingCollection<VariantContext> sortingBuffer;
    protected VariantContextWriter writer;


    @Override
    protected Predicate<VariantContext> makeTruthVariantFilter() {
        return vc -> true;
    }

    protected VCFHeader createHeader(final VCFHeader inputHeader) {
        final VCFHeader header = new VCFHeader(inputHeader.getMetaDataInInputOrder());
        header.setSequenceDictionary(dictionary);

        header.addMetaDataLine(Concordance.TRUTH_STATUS_HEADER_LINE);
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_VARIANT_ID_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Matching truth set variant id"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_ALLELE_COUNT_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Truth set allele count"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_ALLELE_NUMBER_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Truth set allele number"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_ALLELE_FREQUENCY_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Truth set allele frequency"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_RECIPROCAL_OVERLAP_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Reciprocal overlap with truth variant"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_SIZE_SIMILARITY_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Size similarity with truth variant"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_DISTANCE_START_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Start coordinate distance in bp to truth variant's start"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.TRUTH_DISTANCE_END_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "End coordinate distance in bp to truth variant's end"));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        SVStratify.addStratifyMetadata(header);
        return header;
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
        sortingBuffer = SortingCollection.newInstance(
                VariantContext.class,
                new VCFRecordCodec(header, true),
                header.getVCFRecordComparator(),
                maxRecordsInRam,
                tmpDir.toPath());

        // Load stratification groups
        if ((stratArgs.configFile == null) ^ (strataClusteringConfigFile == null)) {
            throw new UserException.BadInput("Both --" + OptionalSVStratificationEngineArgumentsCollection.STRATIFY_CONFIG_FILE_LONG_NAME
                    + " and --" + GroupedSVCluster.CLUSTERING_CONFIG_FILE_LONG_NAME + " must be used together, but only one was specified.");
        }
        final Map<String, CorrespondingSVSelector> clusterEngineMap = new HashMap<>();
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
                    CorrespondingSVSelector engine;
                    engine = new SVMatchesFinder(linkage);
                    clusterEngineMap.put(parameters.name(), engine);
                }
            } catch (final IOException e) {
                throw new GATKException("IO error while reading config table", e);
            }
        }

        final SVConcordanceLinkage defaultLinkage = new SVConcordanceLinkage(dictionary);
        defaultLinkage.setDepthOnlyParams(defaultClusteringArgs.getDepthParameters());
        defaultLinkage.setMixedParams(defaultClusteringArgs.getMixedParameters());
        defaultLinkage.setEvidenceParams(defaultClusteringArgs.getPESRParameters());
        CorrespondingSVSelector defaultEngine;
        defaultEngine = new SVMatchesFinder(defaultLinkage);

        final SVStratificationEngine stratEngine = SVStratify.loadStratificationConfig(stratArgs.configFile, stratArgs, dictionary);
        engine = new StratifiedConcordanceEngine(clusterEngineMap, defaultEngine, stratEngine, stratArgs);
    }

    @Override
    public Object onTraversalSuccess() {
        for (final VariantContext variant : engine.flush(true)) {
            sortingBuffer.add(variant);
        }
        if (!engine.isEmpty()) {
            throw new GATKException("Concordance engine is not empty, but it should be");
        }
        for (final VariantContext variant : sortingBuffer) {
            writer.add(variant);
        }
        return super.onTraversalSuccess();
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
    protected boolean areVariantsAtSameLocusConcordant(final VariantContext truth, final VariantContext eval) {
        return true;
    }

    @Override
    public void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {
        if (truthVersusEval.hasTruth()) {
            final SVCallRecord record = SVCallRecordUtils.create(truthVersusEval.getTruth(), true, false, dictionary);
            engine.addTruthVariant(record);
        }
        if (truthVersusEval.hasEval()) {
            final SVCallRecord record = SVCallRecordUtils.create(truthVersusEval.getEval(), true, false, dictionary);
            engine.addEvalVariant(record);
        }
        for (final VariantContext variant : engine.flush(false)) {
            sortingBuffer.add(variant);
        }
    }

}
