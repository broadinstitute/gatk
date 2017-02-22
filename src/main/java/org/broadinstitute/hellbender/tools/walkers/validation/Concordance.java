package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.collections4.Predicate;
import org.apache.commons.lang.mutable.MutableLong;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.variantutils.VariantsToTable;

import java.io.File;
import java.io.IOException;
import java.util.EnumMap;

/**
 *
 * This tool evaluates a vcf against a validated truth vcf. We assume that the truth vcf only contains PASS variants.
 * The summary statistics (# true positives, # false positives, # false negatives, sensitivity, precision)
 * are reported as a summary tsv (--summary).
 *
 * Optionally, the tool also produces vcfs of
 * 1) true positives and false negatives (i.e. all variants in the truth vcf)
 *      This mode is useful for calculating sensitivity
 * 2) true positives and false positives (i.e. all variants in the eval vcf)
 *      This mode is useful for obtaining a training data set for machine learning classifiers of artifacts.
 *
 * These vcfs can be passed to {@link VariantsToTable} to produce a tsv for statistical analysis in R or Python.
 *
 *
 * java -jar gatk.jar Concordance -eval na12878-eval.vcf --truth na12878-truth.vcf  --summary summary.tsv
 *
 * Created by Takuto Sato on 1/30/17.
 */
@CommandLineProgramProperties(
        summary = "Evaluate a vcf against a vcf of validated (true) variants",
        oneLineSummary = "Evaluate a vcf against a vcf of validated (true) variants",
        programGroup = VariantProgramGroup.class
)
public class Concordance extends AbstractConcordanceWalker {
    public static final String SUMMARY_LONG_NAME = "summary";
    public static final String SUMMARY_SHORT_NAME = "S";

    public static final String TRUE_POSITIVES_AND_FALSE_NEGATIVES_LONG_NAME = "truePositivesAndFalseNegatives";
    public static final String TRUE_POSITIVES_AND_FALSE_NEGATIVES_SHORT_NAME = "tpfn";
    public static final String TRUE_POSITIVES_AND_FALSE_POSITIVES_LONG_NAME = "truePositivesAndFalsePositives";
    public static final String TRUE_POSITIVES_AND_FALSE_POSITIVES_SHORT_NAME = "tpfp";

    public static final String TRUTH_STATUS_VCF_ATTRIBUTE = "STATUS";
    private static VCFInfoHeaderLine TRUTH_STATUS_HEADER_LINE =
            new VCFInfoHeaderLine(TRUTH_STATUS_VCF_ATTRIBUTE, 1,VCFHeaderLineType.String, "Truth status: TP/FP/FN for true positive/false positive/false negative.");

    @Argument(doc = "A table of summary statistics (true positives, sensitivity, etc.)",
            fullName = SUMMARY_LONG_NAME,
            shortName = SUMMARY_SHORT_NAME)
    protected File summary;

    @Argument(doc = "A vcf to write true positives and false negatives",
            fullName = TRUE_POSITIVES_AND_FALSE_NEGATIVES_LONG_NAME,
            shortName = TRUE_POSITIVES_AND_FALSE_NEGATIVES_SHORT_NAME,
            optional = true)
    protected File truePositivesAndFalseNegativesVcf = null;

    @Argument(doc = "A vcf to write true positives and false positives",
            fullName = TRUE_POSITIVES_AND_FALSE_POSITIVES_LONG_NAME,
            shortName = TRUE_POSITIVES_AND_FALSE_POSITIVES_SHORT_NAME,
            optional = true)
    protected File truePositivesAndFalsePositivesVcf = null;

    // we count true positives, false positives, false negatives for snps and indels
    private final EnumMap<ConcordanceState, MutableLong> snpCounts = new EnumMap<>(ConcordanceState.class);
    private final EnumMap<ConcordanceState, MutableLong> indelCounts = new EnumMap<>(ConcordanceState.class);
    private VariantContextWriter truePositivesAndFalseNegativesVcfWriter;
    private VariantContextWriter truePositivesAndFalsePositivesVcfWriter;

    @Override
    public void onTraversalStart() {
        for (final ConcordanceState state : ConcordanceState.values()) {
            snpCounts.put(state, new MutableLong(0));
            indelCounts.put(state, new MutableLong(0));
        }

        if (truePositivesAndFalseNegativesVcf != null) {
            truePositivesAndFalseNegativesVcfWriter = createVCFWriter(truePositivesAndFalseNegativesVcf);
            final VCFHeader truthHeader = getTruthHeader();
            truthHeader.addMetaDataLine(TRUTH_STATUS_HEADER_LINE);
            truePositivesAndFalseNegativesVcfWriter.writeHeader(truthHeader);
        }

        if (truePositivesAndFalsePositivesVcf != null) {
            truePositivesAndFalsePositivesVcfWriter = createVCFWriter(truePositivesAndFalsePositivesVcf);
            final VCFHeader evalHeader = getEvalHeader();
            evalHeader.addMetaDataLine(TRUTH_STATUS_HEADER_LINE);
            truePositivesAndFalsePositivesVcfWriter.writeHeader(evalHeader);
        }
    }

    @Override
    protected void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext, final ReferenceContext refContext) {
        final ConcordanceState concordanceState = truthVersusEval.getConcordance();
        if (truthVersusEval.getTruthIfPresentElseEval().isSNP()) {
            snpCounts.get(concordanceState).increment();
        } else {
            indelCounts.get(concordanceState).increment();
        }

        if (truePositivesAndFalseNegativesVcfWriter != null && concordanceState != ConcordanceState.FALSE_POSITIVE) {
            final VariantContext vc = annotateWithConcordanceState(truthVersusEval.getTruth(), concordanceState);
            truePositivesAndFalseNegativesVcfWriter.add(vc);
        }

        if (truePositivesAndFalsePositivesVcfWriter != null && concordanceState != ConcordanceState.FALSE_NEGATIVE) {
            final VariantContext vc = annotateWithConcordanceState(truthVersusEval.getEval(), concordanceState);
            truePositivesAndFalsePositivesVcfWriter.add(vc);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        try ( ConcordanceSummaryRecord.Writer concordanceSummaryWriter = ConcordanceSummaryRecord.getWriter(summary) ){
            concordanceSummaryWriter.writeRecord(new ConcordanceSummaryRecord(VariantContext.Type.SNP,
                    snpCounts.get(ConcordanceState.TRUE_POSITIVE).longValue(),
                    snpCounts.get(ConcordanceState.FALSE_POSITIVE).longValue(),
                    snpCounts.get(ConcordanceState.FALSE_NEGATIVE).longValue()));
            concordanceSummaryWriter.writeRecord(new ConcordanceSummaryRecord(VariantContext.Type.INDEL,
                    indelCounts.get(ConcordanceState.TRUE_POSITIVE).longValue(),
                    indelCounts.get(ConcordanceState.FALSE_POSITIVE).longValue(),
                    indelCounts.get(ConcordanceState.FALSE_NEGATIVE).longValue()));
        } catch (IOException e){
            throw new UserException("Encountered an IO exception writing the concordance summary table", e);
        }

        if (truePositivesAndFalsePositivesVcfWriter != null) {
            truePositivesAndFalsePositivesVcfWriter.close();
        }

        if (truePositivesAndFalseNegativesVcfWriter != null) {
            truePositivesAndFalseNegativesVcfWriter.close();
        }

        return "SUCCESS";
    }

    @Override
    protected boolean areVariantsAtSameLocusConcordant(final VariantContext truth, final VariantContext eval) {
        final boolean sameRefAllele = truth.getReference().equals(eval.getReference());
        // we assume that the truth has a single alternate allele
        final boolean containsAltAllele = eval.getAlternateAlleles().contains(truth.getAlternateAllele(0));

        return sameRefAllele && containsAltAllele;
    }

    @Override
    protected Predicate<VariantContext> makeVariantFilter() {
        return vc -> !vc.isFiltered() && ! vc.isSymbolicOrSV();
    }

    private VariantContext annotateWithConcordanceState(final VariantContext vc, final ConcordanceState state) {
        return new VariantContextBuilder(vc).attribute(TRUTH_STATUS_VCF_ATTRIBUTE, state.getAbbreviation()).make();
    }
}
