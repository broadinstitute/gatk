package org.broadinstitute.hellbender.tools.walkers.concordance;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import org.apache.commons.collections4.iterators.FilterIterator;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import static org.broadinstitute.hellbender.tools.walkers.concordance.VariantStatusRecord.Status.*;

/**
 *
 * This tool evaluates a vcf against a validated truth vcf. We assume that the truth vcf only contains PASS variants.
 * The summary statistics (# true positives, # false positives, # false negatives, sensitivity, precision)
 * are reported as a summary tsv (--summary). Also reported as a tsv (--table) are evaluated variants of interest
 * (i.e. true positives, false positives, false negatives) and their annotations (ref allele, alt allele, minor allele fraction, etc.)
 * for further statistical analysis e.g. with an R script
 *
 * java -jar gatk.jar Concordance -V na12878-eval.vcf --truth na12878-truth.vcf --table table.tsv --summary summary.tsv --evalSampleName NA12878
 *
 * Created by Takuto Sato on 1/30/17.
 */

@CommandLineProgramProperties(
        summary = "Evaluate a vcf against a vcf of validated (true) variants",
        oneLineSummary = "Evaluate a vcf against a vcf of validated (true) variants",
        programGroup = VariantProgramGroup.class
)

public class Concordance extends VariantWalker {
    public static final String CONFIDENCE_REGION_LONG_NAME = "confidence";
    public static final String CONFIDENCE_REGION_SHORT_NAME = "C";
    public static final String SAMPLE_LONG_NAME = "sampleName";
    public static final String SAMPLE_SHORT_NAME = "sample";
    public static final String SUMMARY_LONG_NAME = "summary";
    public static final String SUMMARY_SHORT_NAME = "S";
    public static final String TRUTH_LONG_NAME = "truth";
    public static final String TRUTH_SHORT_NAME = "T";

    @Argument(doc = "truth vcf (tool assumes all sites in truth are PASS)",
            fullName= TRUTH_LONG_NAME,
            shortName = TRUTH_SHORT_NAME)
    protected File truth;

    @Argument(doc = "TO BE IMPLEMENTED",
            fullName= CONFIDENCE_REGION_LONG_NAME,
            shortName = CONFIDENCE_REGION_SHORT_NAME,
            optional = true)
    protected File highConfidenceRegion;

    @Argument(doc = "A table of evaluation variants of interest and their basic annotations",
            fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    protected File table;

    @Argument(doc = "A table of summary statistics (true positives, sensitivity, etc.)",
            fullName= SUMMARY_LONG_NAME,
            shortName = SUMMARY_SHORT_NAME)
    protected File summary;

    @Argument(doc = "sample name of the eval variants",
            fullName = SAMPLE_LONG_NAME,
            shortName = SAMPLE_SHORT_NAME)
    protected String evalSampleName;

    // TODO: output a vcf of false positives (and negative) sites?
    // TODO: take a low confidence region where a false positive there is not counted (masked in dream?)
    // TODO: what about other variant types besides SNP and INDEL?

    // we count evaluation status in these variables
    private long snpTruePositives = 0;
    private long snpFalsePositives = 0;
    private long snpFalseNegatives = 0;
    private long indelTruePositives = 0;
    private long indelFalsePositives = 0;
    private long indelFalseNegatives = 0;

    private VariantStatusRecord.Writer variantStatusWriter;


    private Iterator<VariantContext> truthIterator;
    private VariantContext currentTruthVariant;
    private long currentSumOfTPandFN = 0;
    private long numVariantsInTruth;
    private boolean exhaustedTruthVariants = false;



    private VariantContextComparator variantContextComparator;

    @Override
    public void onTraversalStart() {
        Utils.regularReadableUserFile(truth);
        Utils.validateArg(getHeaderForVariants().getSampleNamesInOrder().contains(evalSampleName),
                String.format("the eval sample %s does not exist in vcf", evalSampleName));

        variantStatusWriter = VariantStatusRecord.getWriter(table);

        SAMSequenceDictionary dictionary = getSequenceDictionaryForDrivingVariants();
        variantContextComparator = new VariantContextComparator(dictionary);

        final FeatureDataSource<VariantContext> truthSource = new FeatureDataSource<>(truth);
        numVariantsInTruth = truthSource.spliterator().getExactSizeIfKnown();
        // the iterator skips the SV and Symbolic variants
        truthIterator = new FilterIterator<>(truthSource.iterator(), vc -> ! vc.isSymbolicOrSV());
        if (truthIterator.hasNext()){
            currentTruthVariant = truthIterator.next();
        } else {
            throw new UserException(String.format("truth file %s does not have any variants", truth.toString()));
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        // TODO: When we give the tool an interval List (-L), variant walker subsets the range of driving variants
        // to these intervals for you behind the scene. But it doesn't do the same for the truth vcf
        // i.e. truth vcf starts at chr 1 when -L 21.

        // Case 0: We've reached the end of truth variants; eval variant is a FP unless filtered
        if ( exhaustedTruthVariants ){
            if (variant.isNotFiltered()){
                processVariant(variant, FALSE_POSITIVE);
            }
            return;
        }

        // case 1: Eval leapfrogged truth; we missed at least one variant.
        // Move the truth cursor forward until it catches up to eval.
        // Note we must check for this first, *then* check for the other two conditions;
        // otherwise we move the eval forward without giving it a diagnosis
        while (variantContextComparator.compare(variant, currentTruthVariant) > 0) {
            if (currentTruthVariant.isSNP()) {
                snpFalseNegatives++;
            } else {
                indelFalseNegatives++;
            }

            if (truthIterator.hasNext()) {
                checkInvariant();
                currentTruthVariant = truthIterator.next();
            } else {
                exhaustedTruthVariants = true;
                if (variant.isNotFiltered()){
                    processVariant(variant, FALSE_POSITIVE);
                }
                return;
            }
        }

        // Now the eval is either at the same position as the truth or ahead of it
        // Case 2: the position of the two variants match
        if (variantContextComparator.compare(variant, currentTruthVariant) == 0) {
            if (variant.isFiltered()) {
                processVariant(variant, FALSE_NEGATIVE);
            } else if (allelesMatch(variant, currentTruthVariant)) {
                processVariant(variant, TRUE_POSITIVE);
            } else {
                // the eval variant matches truth's position but their alleles don't match
                // e.g. genotype or alleles don't match
                // Increment the counts for both false positive and false negative
                // but the record will only say FALSE_POSITIVE (arbitrary design decision)
                processVariant(variant, FALSE_POSITIVE);

                if (variant.isSNP()) {
                    snpFalseNegatives++;
                } else {
                    indelFalseNegatives++;
                }
            }

            if (truthIterator.hasNext()) {
                checkInvariant();
                currentTruthVariant = truthIterator.next();
            } else {
                exhaustedTruthVariants = true;
            }

            return;
        }

        // case 3: truth got ahead of eval; the eval variant must be a false positive if not filtered
        if (variantContextComparator.compare(variant, currentTruthVariant) < 0) {
            if (variant.isNotFiltered()) {
                processVariant(variant, FALSE_POSITIVE);
            }

            // we don't care about true negatives - don't record them in the table
            return;
        }
    }

    /***
     *  When we went to increment the truth, check that we incremented TP + FN by exactly 1
     *  Call this right before we call truthIterator.next()
     */
    private void checkInvariant(){
        final long falseNegatives = snpFalseNegatives + indelFalseNegatives;
        final long truePositives = snpTruePositives + indelTruePositives;
        assert truePositives + falseNegatives == currentSumOfTPandFN + 1;
        currentSumOfTPandFN = truePositives + falseNegatives;
    }

    // Given a variant and its truth status, A) increment the counters and B) write the record
    private void processVariant(final VariantContext variant, final VariantStatusRecord.Status status){
         switch (status) {
            case TRUE_POSITIVE:
                if (variant.isSNP()) {
                    snpTruePositives++;
                } else {
                    indelTruePositives++;
                }
                break;
            case FALSE_POSITIVE:
                if (variant.isSNP()) {
                    snpFalsePositives++;
                } else {
                    indelFalsePositives++;
                }
                break;
            case FALSE_NEGATIVE:
                if (variant.isSNP()){
                    snpFalseNegatives++;
                } else {
                    indelFalseNegatives++;
                }
                break;
        }

        // Write the record to the output table
        // TODO: When the time comes we should handle the multi-allelic case more carefully
        final double[] alleleFractions = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(variant.getGenotype(evalSampleName),
                GATKVCFConstants.ALLELE_FRACTION_KEY, () -> new double[] {-1.0}, -1.0);
        final double maxAlleleFraction = MathUtils.arrayMax(alleleFractions);
        VariantStatusRecord record = new VariantStatusRecord(variant, status, maxAlleleFraction);
        try {
            variantStatusWriter.writeRecord(record);
        } catch (Exception e) {
            throw new UserException(String.format("Encountered an IO exception writing a record at chrom %d, pos ",
                    variant.getContig(), variant.getStart()), e);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (! exhaustedTruthVariants ) {
            // we land here off of case 3 above: that is, truth was ahead of eval and we ran out of evals
            // remaining variants in truth are false negatives
            while (true) {
                // first process the current truth variant, which we know to be non-empty
                if (currentTruthVariant.isSNP()) {
                    snpFalseNegatives++;
                } else {
                    indelFalseNegatives++;
                }

                // then move up the iterator if there are more elements
                if (truthIterator.hasNext()) {
                    checkInvariant();
                    currentTruthVariant = truthIterator.next();
                } else {
                    break;
                }
            }
        }

        try ( ConcordanceSummaryRecord.Writer concordanceSummaryWriter = ConcordanceSummaryRecord.getWriter(summary) ){
            concordanceSummaryWriter.writeRecord(new ConcordanceSummaryRecord(VariantContext.Type.SNP, snpTruePositives, snpFalsePositives, snpFalseNegatives));
            concordanceSummaryWriter.writeRecord(new ConcordanceSummaryRecord(VariantContext.Type.INDEL, indelTruePositives, indelFalsePositives, indelFalseNegatives));

            // We must close the variant status writer to flush the buffer.
            // Summary writer will be closed implicitly
            variantStatusWriter.close();
        } catch (IOException e){
            throw new UserException("Encountered an IO exception writing the concordance summary table", e);
        }

        return "SUCCESS";
    }

    // TODO: eventually this should be an abstract method to be overridden by a subclass
    // such that the user can adjust the stringency of allele comparison
    protected static boolean allelesMatch(final VariantContext eval, final VariantContext truth){
        final boolean sameContig = truth.getContig().equals(eval.getContig());
        final boolean sameStartPosition = truth.getStart() == eval.getStart();
        final boolean sameRefAllele = truth.getReference().equals(eval.getReference());
        // we assume that the truth has a single alternate allele
        final boolean containsAltAllele = eval.getAlternateAlleles().contains(truth.getAlternateAllele(0));

        return sameContig && sameStartPosition && sameRefAllele && containsAltAllele;
    }


}
