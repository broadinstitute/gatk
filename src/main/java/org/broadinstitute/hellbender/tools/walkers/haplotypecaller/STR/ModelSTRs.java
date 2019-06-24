package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by valentin on 11/3/16.
 */
public class ModelSTRs extends STRAnalysisWalker {

    @Argument(doc = "sites where the major allele balance is less than this fraction will be ignored", shortName = "mMAB", fullName = "minimumMajorAlleleBalance", optional = true)
    public double minimumMajorAlleleBalance = 0.6;

    @Argument(doc = "number of bp around the STR to be considered when calculating the GC content", shortName = "gcPad", fullName = "gcContentPadding", optional = true)
    public int gcContentPadding = 10;

    @Argument(doc = "output Context table", shortName = "O", fullName = "output")
    public File output;

    private PrintWriter outputWriter;

    @Argument(doc = "maximum repeat count diff reported in output", optional = true)
    public int maximumRepeatCountDiffReportedInOutput = 10;

    @Argument(doc = "minimum repeat count diff reportedd in output", optional = true)
    public int minimumRepeatCountDiffReportedInOutput = -10;

    @Argument(doc = "dont include non-variant blocks", optional = true)
    private boolean dontIncludeNonVariantsBlocks = false;

    private String lastOutputContig;
    private int lastOutputPosition;

    @Argument(doc = "minimum GQ to trust a hom/ref call", optional = true)
    private int minHomRefGQ = 30;
    @Argument(doc = "minimum GQ to trust truth calls", optional = true)
    private int minTruthGQ = 30;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        try {
            outputWriter = new PrintWriter(new FileWriter(output));
            final String[] fixedColumnNames = new String[]{
                    "CONTIG", "POS", "GC", "DEPTH", "UNIT", "U.LENGTH", "MA.RC", "T.LENGTH", "CORRECT", "ERRORS", "Minus", "Plus", "Minus.overflow"
            };
            final List<String> columnNames = new ArrayList<>();
            for (final String string : fixedColumnNames) {
                columnNames.add(string);
            }
            for (int i = minimumRepeatCountDiffReportedInOutput; i <= maximumRepeatCountDiffReportedInOutput; i++) {
                final String prefix = (i < 0) ? "Minus." : (i == 0) ? "" : "Plus.";
                columnNames.add(prefix + Math.abs(i));
            }
            columnNames.add("Plus.overflow");
            outputWriter.println(String.join("\t", columnNames));
            outputWriter.flush();
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(output, e);
        }
    }

    private STRContext map(final VariantContext call, final FeatureContext features, final ReferenceContext ref) {
        if (ref == null || call == null || features == null) {
            return null;
        }
        final VariantContext truth = findBestMatch(call, getTruthFeatures(features));

        if (isNonVariantBlock(call)) {
            if (dontIncludeNonVariantsBlocks) {
                return null;
            }
            final STRContext result = model.composeContext(ref, call);
            if (call.getGenotype(0).getGQ() < minHomRefGQ) {
                logger.debug("Discarded STR at non-variant block since hom/ref call is not trustworthy " + result);
                return null;
            }
            if (result == null) {
                return null;
            } else {
                if (truth != null && !truth.getGenotype(0).isHom()) { // ignore het sites.
                    logger.debug("Discarded due to het truth: " + call);
                    return null;
                } else if (truth != null && truthGQ(truth) < minTruthGQ) {
                    logger.debug("Discarded due to lack of confidence in truth call:" + call + " " + truth);
                    return null;
                } else if (truth != null && !truth.getGenotype(0).isHomRef()) {
                    logger.debug("Discarded due to lack of confidence in truth call hom-ref: " + truth);
                    return null;
                }

                result.setGCContent(calculateGCContent(ref, result));
                return result;
            }
        } else if (call.getStart() != ref.getInterval().getStart()) {
            return null;
        } else {
            // might return null if no repeat unit could be determined for the site.
            final STRContext result = model.composeContext(ref, call);
            if (result == null) {
                return null;
            } else {
                if (truth != null && !truth.getGenotype(0).isHom()) { // ignore het sites.
                    logger.debug("Discarded due to het truth: " + call);
                    return null;
                } else if (truth != null && truthGQ(truth) < minTruthGQ) {
                    logger.debug("Discarded due to lack of confidence in truth call:" + call + " " + truth);
                    return null;
                }
                final int truthAlleleIndex = truth == null ? -1 : result.getAlleles().alleleIndex(truth.getReference(), truth.getGenotype(0).getAllele(0));
                final STRAllele truthAllele = truthAlleleIndex == -1 ? null : result.getAlleles().get(truthAlleleIndex);
                if (truthAllele == null || result.getMajorAlleleRepeatCount() != truthAllele.repeatCount) {
                    logger.debug("Discarded due to mismatching truth: " + call);
                    return null; // skip cases where the truth allele is not even the major allele.
                }

                final double majorAlleleBalance = result.getMajorAlleleBalance();
                if (majorAlleleBalance < minimumMajorAlleleBalance) {
                    logger.debug("Discarded due to low MAB: " + call);
                    return null;
                    //} else if (majorAlleleBalance == 1.0) { // ignore sites where only one allele has AD != 0.
                    //    logger.debug("Discarded due to high MAB == 1.0: " + vc);
                    //    return null;
                } else {
                    result.setGCContent(calculateGCContent(ref, result));
                    return result;
                }
            }
        }
    }

    private int truthGQ(final VariantContext truth) {
        if (truth == null) {
            return 0;
        } else if (truth.getGenotype(0).hasGQ()) {
            return truth.getGenotype(0).getGQ();
        } else {
            return (int) Math.floor(truth.getPhredScaledQual());
        }
    }

    private static boolean isNonVariantBlock(final VariantContext vc) {
        if (vc == null) {
            return false;
        } else if (vc.getAlternateAlleles().size() == 0) {
            return true;
        } else if (vc.getAlternateAlleles().size() == 1) {
            return vc.getAlternateAllele(0).equals(Allele.NON_REF_ALLELE);
        } else {
            return false;
        }
    }

    private double calculateGCContent(final ReferenceContext ref, final STRContext context) {
        final int startOffset = -gcContentPadding;
        final int stopOffset = context.getAlleles().getRepeatUnitLength() * context.getMajorAlleleRepeatCount() + gcContentPadding - 1;
        final SimpleInterval window = ref.getWindow();
        final int startOffsetInWindow = Math.max(window.getStart(), ref.getInterval().getStart() + startOffset) - window.getStart();
        final int stopOffsetInWindow = Math.min(window.getEnd(), ref.getInterval().getStart() + stopOffset) - window.getStart();
        final byte[] windowBases = ref.getBases();
        int gcCount = 0;
        int atCount = 0;
        for (int i = startOffsetInWindow; i <= stopOffsetInWindow; i++) {
            switch (Character.toLowerCase(windowBases[i])) {
                case 'c':
                case 'g': gcCount++; break;
                case 'a':
                case 't': atCount++; break;
            }
        }
        if (gcCount == atCount) { // cover for the case that both are zero (i.e. all N?); in practice this would never happen.
            return 0.5;
        } else {
            return gcCount / ((double) atCount + gcCount);
        }
    }

    private void outputContext(final STRContext value) {
        if (value.getRepeatLocus().getContig().equals(lastOutputContig) && lastOutputPosition >= value.getRepeatLocus().getStart()) {
            return; // skip overlapping AssessSTRs.
        }
        lastOutputContig = value.getRepeatLocus().getContig();
        lastOutputPosition = value.getRepeatLocus().getStart() + value.getMaximumRepeatCount() * value.getAlleles().getRepeatUnitLength();
        final int[] diffADs = new int[maximumRepeatCountDiffReportedInOutput - minimumRepeatCountDiffReportedInOutput + 1];
        final List<String> fields = new ArrayList<>(diffADs.length + 20);
        Arrays.fill(diffADs, 0);
        fields.clear();
        final int[] ads = value.getAlleleDepths();
        final STRAlleleSet alleles = value.getAlleles();
        final int majorAlleleRepeatCount = value.getMajorAlleleRepeatCount();
        int errorMinus = 0;
        int errorPlus = 0;
        int minusOverflow = 0;
        int plusOverflow = 0;
        for (int i = 0; i < ads.length; i++) {
            final int diffADsIndex = alleles.get(i).repeatCount - majorAlleleRepeatCount - minimumRepeatCountDiffReportedInOutput;
            if (diffADsIndex < 0) {
                minusOverflow++;
            } else if (diffADsIndex >= diffADs.length) {
                plusOverflow++;
            } else {
                diffADs[diffADsIndex] = ads[i];
            }
            if (alleles.get(i).repeatCount < majorAlleleRepeatCount)
                errorMinus += ads[i];
            else if (alleles.get(i).repeatCount > majorAlleleRepeatCount)
                errorPlus += ads[i];
        }
        final long adSum = MathUtils.sum(value.getAlleleDepths());
        fields.add(value.getRepeatLocus().getContig());
        fields.add("" + value.getRepeatLocus().getStart());
        fields.add("" + value.getGCContent());
        fields.add("" + MathUtils.sum(value.getAlleleDepths()));
        fields.add("" + value.getAlleles().getRepeatUnitString());
        fields.add("" + value.getAlleles().getRepeatUnitLength());
        fields.add("" + value.getMajorAlleleRepeatCount());
        fields.add("" + value.getMajorAlleleRepeatCount() * value.getAlleles().getRepeatUnitLength());
        fields.add("" + Math.round(value.getMajorAlleleBalance() * adSum));
        fields.add("" + (adSum - Math.round(value.getMajorAlleleBalance() * adSum)));
        fields.add("" + errorMinus);
        fields.add("" + errorPlus);
        fields.add("" + minusOverflow);
        for (int i = 0; i < diffADs.length; i++) {
            fields.add("" + diffADs[i]);
        }
        fields.add("" + plusOverflow);
        outputWriter.println(String.join("\t", fields));
    }

    @Override
    public Object onTraversalSuccess() {
        outputWriter.flush();
        return "SUCCESS";
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (isCall(variant)) {
            final STRContext result = map(variant, featureContext, referenceContext);
            if (result != null) {
                outputContext(result);
            }
        }
    }
}
