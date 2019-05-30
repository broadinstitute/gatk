package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.STR;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.walkers.*;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.commandline.*;
import org.broadinstitute.hellbender.utils.contexts.AlignmentContext;
import org.broadinstitute.hellbender.utils.contexts.ReferenceContext;
import org.broadinstitute.hellbender.utils.exceptions.UserException;
import org.broadinstitute.hellbender.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

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
@Requires({DataSource.REFERENCE, DataSource.REFERENCE_ORDERED_DATA})
@Reference(window=@Window(start = -100, stop = 1000))
public class ModelSTRs extends RodWalker<STRContext, Integer> {

    @Argument(doc = "sites where the major allele balance is less than this fraction will be ignored", shortName = "mMAB", fullName = "minimumMajorAlleleBalance", required = false)
    public double minimumMajorAlleleBalance = 0.6;

    @Argument(doc = "number of bp around the STR to be considered when calculating the GC content", shortName = "gcPad", fullName = "gcContentPadding", required = false)
    public int gcContentPadding = 10;

    @Input(doc = "input VCF", shortName = "V", fullName = "variant")
    public RodBinding<VariantContext> input;

    @Input(doc = "truth VCF calls", shortName = "truth", fullName = "truthVariant", required = false)
    public RodBinding<VariantContext> truth;

    @ArgumentCollection
    public STRModel model = new STRModel();

    @Output(doc = "output Context table", shortName = "O", fullName = "output")
    public File output;

    private PrintWriter outputWriter;

    @Argument(doc = "maximum repeat count diff reported in output", required = false)
    public int maximumRepeatCountDiffReportedInOutput = 10;

    @Argument(doc = "minimum repeat count diff reportedd in output", required = false)
    public int minimumRepeatCountDiffReportedInOutput = -10;

    @Argument(doc = "dont include non-variant blocks", required = false)
    private boolean dontIncludeNonVariantsBlocks = false;

    private String lastOutputContig;
    private int lastOutputPosition;

    @Argument(doc = "minimum GQ to trust a hom/ref call", required = false)
    private int minHomRefGQ = 30;
    @Argument(doc = "minimum GQ to trust truth calls", required = false)
    private int minTruthGQ = 30;

    @Override
    public void initialize() {
        super.initialize();
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

    @Override
    public STRContext map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if (tracker == null || ref == null) {
            return null;
        }
        final VariantContext vc = tracker.getFirstValue(input);
        final VariantContext truthVC = findTruthVC(tracker.getValues(truth), vc);

        if (vc == null) {
            return null;
        } else if (isNonVariantBlock(vc)) {
            if (dontIncludeNonVariantsBlocks) {
                return null;
            }
            final STRContext result = model.composeContext(ref, vc);
            if (vc.getGenotype(0).getGQ() < minHomRefGQ) {
                logger.debug("Discarded STR at non-variant block since hom/ref call is not trustworthy " + result);
                return null;
            }
            if (result == null) {
                return null;
            } else {
                if (truthVC != null && !truthVC.getGenotype(0).isHom()) { // ignore het sites.
                    logger.debug("Discarded due to het truth: " + vc);
                    return null;
                } else if (truthVC != null && truthGQ(truthVC) < minTruthGQ) {
                    logger.debug("Discarded due to lack of confidence in truth call:" + vc + " " + truthVC);
                    return null;
                } else if (truthVC != null && !truthVC.getGenotype(0).isHomRef()) {
                    logger.debug("Discarded due to lack of confidence in truth call hom-ref: " + truthVC);
                    return null;
                }

                result.setGCContent(calculateGCContent(ref, result));
                return result;
            }
        } else if (vc.getStart() != ref.getLocus().getStart()) {
            return null;
        } else {
            // might return null if no repeat unit could be determined for the site.
            final STRContext result = model.composeContext(ref, vc);
            if (result == null) {
                return null;
            } else {
                if (truthVC != null && !truthVC.getGenotype(0).isHom()) { // ignore het sites.
                    logger.debug("Discarded due to het truth: " + vc);
                    return null;
                } else if (truthVC != null && truthGQ(truthVC) < minTruthGQ) {
                    logger.debug("Discarded due to lack of confidence in truth call:" + vc + " " + truthVC);
                    return null;
                }
                final int truthAlleleIndex = truthVC == null ? (truth == null ? -1 : result.getAlleles().alleleIndex(result.getAlleles().getReference())) : result.getAlleles().alleleIndex(truthVC.getReference(), truthVC.getGenotype(0).getAllele(0));
                final STRAllele truthAllele = truthAlleleIndex == -1 ? null : result.getAlleles().get(truthAlleleIndex);
                if (truthAllele == null || result.getMajorAlleleRepeatCount() != truthAllele.repeatCount) {
                    logger.debug("Discarded due to mismatching truth: " + vc);
                    return null; // skip cases where the truth allele is not even the major allele.
                }

                final double majorAlleleBalance = result.getMajorAlleleBalance();
                if (majorAlleleBalance < minimumMajorAlleleBalance) {
                    logger.debug("Discarded due to low MAB: " + vc);
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

    private VariantContext findTruthVC(final List<VariantContext> truths, final VariantContext vc) {
        if (truths == null || truths.isEmpty()) {
            return null;
        }
        for (final VariantContext truth : truths) {
            if (vc.getStart() == truth.getStart()) {
                return truth;
            }
        }
        for (final VariantContext truth : truths) {
            if (truth.getAlternateAlleles().isEmpty() || (truth.getAlternateAlleles().size() == 1 && truth.getAlternateAlleles().get(0).equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE))) {
                return truth;
            }
        }
        return null;
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
            return vc.getAlternateAllele(0).equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
        } else {
            return false;
        }
    }

    private double calculateGCContent(final ReferenceContext ref, final STRContext context) {
        final int startOffset = -gcContentPadding;
        final int stopOffset = context.getAlleles().getRepeatUnitLength() * context.getMajorAlleleRepeatCount() + gcContentPadding - 1;
        final GenomeLoc window = ref.getWindow();
        final int startOffsetInWindow = Math.max(window.getStart(), ref.getLocus().getStart() + startOffset) - window.getStart();
        final int stopOffsetInWindow = Math.min(window.getStop(), ref.getLocus().getStart() + stopOffset) - window.getStart();
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

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(final STRContext value, final Integer count) {
        if (value != null) {
            outputContext(value);
        }
        return count + 1;
    }

    private void outputContext(STRContext value) {
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
    public void onTraversalDone(final Integer count) {
        super.onTraversalDone(count);
        logger.info("Number of STR context found: " + count);
        outputWriter.flush();
    }
}
