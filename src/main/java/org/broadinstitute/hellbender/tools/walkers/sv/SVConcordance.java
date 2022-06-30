package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
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
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.tools.sv.concordance.ClosestSVFinder;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceAnnotator;
import org.broadinstitute.hellbender.tools.sv.concordance.SVConcordanceLinkage;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;
import picard.vcf.GenotypeConcordance;

/**
 * <p>This tool calculates SV genotype concordance between an "evaluation" VCF and a "truth" VCF.</p>
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

    public static final String USE_TRUTH_AF_LONG_NAME = "use-truth-af";
    public static final String BIALLELIC_DUPLICATIONS_LONG_NAME = "force-biallelic-dups";

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFile;

    /**
     * By default, truth allele frequencies are calculated on the fly, using the evaluation record's allele number as the
     * denominator. This option forces the tool to use the allele frequency annotations (AF/AN/AC) of the closest-
     * matching truth variant record (by min distance to both breakpoints) for truth allele frequency annotations.
     */
    @Argument(
            doc = "Use allele frequency annotations from the truth vcf from the best-matching record.",
            fullName = USE_TRUTH_AF_LONG_NAME,
            optional = true
    )
    private boolean useTruthAf = false;

    /**
     * By default, the tool assumes that DUP records have multi-copy alleles (i.e. with a 1-copy allele, 2-copy allele,
     * etc.). In this case, a genotype's alleles cannot be inferred when the copy number greater than 1 + the
     * sample ploidy, and the tool will treat the genotype as a no-call. For example, a genotype with
     * ECN=2 and CN=4 will be treated as a no-call (./.) since this copy state correspond to either a homozygous 1-copy or
     * heterozygous 2-copy / reference genotype. This option overrides that behavior by assuming only a 1-copy allele
     * and can therefore infer the genotype allele for cases where the copy number is greater than 1 + the ploidy.
     * In the previous example, the genotype would be treated as homozygous (1/1). Note that in cases where the copy
     * number is more than twice the ploidy, the genotype is treated as homozygous.
     */
    @Argument(
            doc = "Interpret duplication alleles as single-copy when inferring genotypes from copy number",
            fullName = BIALLELIC_DUPLICATIONS_LONG_NAME,
            optional = true
    )
    private boolean biallelicDups = false;

    @ArgumentCollection
    private final SVClusterEngineArgumentsCollection clusterParameterArgs = new SVClusterEngineArgumentsCollection();

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private SVConcordanceLinkage linkage;
    private ClosestSVFinder engine;
    private String currentContig = null;

    @Override
    public void onTraversalStart() {
        // Use master sequence dictionary i.e. hg38 .dict file since the "best" dictionary is grabbed
        // from the VCF, which is sometimes out of order
        dictionary = getMasterSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }

        linkage = new SVConcordanceLinkage(dictionary);
        linkage.setDepthOnlyParams(clusterParameterArgs.getDepthParameters());
        linkage.setMixedParams(clusterParameterArgs.getMixedParameters());
        linkage.setEvidenceParams(clusterParameterArgs.getPESRParameters());

        final SVConcordanceAnnotator collapser = new SVConcordanceAnnotator(useTruthAf);
        engine = new ClosestSVFinder(linkage, collapser::annotate, dictionary);

        writer = createVCFWriter(outputFile);
        writer.writeHeader(createHeader(getEvalHeader()));
    }

    @Override
    public Object onTraversalSuccess() {
        flushClusters(true);
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (writer != null) {
            writer.close();
        }
    }

    @Override
    public void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {
        if (truthVersusEval.hasTruth()) {
            add(truthVersusEval.getTruth(), true);
        }
        if (truthVersusEval.hasEval()) {
            add(truthVersusEval.getEval(), false);
        }
    }

    private void add(final VariantContext variant, final boolean isTruth) {
        SVCallRecord record = SVCallRecordUtils.create(variant);
        if (!record.getContigA().equals(currentContig)) {
            flushClusters(true);
            currentContig = record.getContigA();
        }
        if (biallelicDups) {
            record = SVCallRecordUtils.convertToBiallelicDupGenotypes(record, this.logger);
        }
        engine.add(record, isTruth);
        flushClusters(false);
    }

    @Override
    protected boolean areVariantsAtSameLocusConcordant(final VariantContext truth, final VariantContext eval) {
        return true;
    }

    private void flushClusters(final boolean force) {
        engine.flush(force).stream()
                .map(SVCallRecordUtils::getVariantBuilder)
                .map(VariantContextBuilder::make)
                .forEach(writer::add);
    }

    private VCFHeader createHeader(final VCFHeader header) {
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
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        return header;
    }


}
