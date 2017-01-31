package org.broadinstitute.hellbender.tools.walkers.mutect;


import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

@CommandLineProgramProperties(
        summary = "Count ",
        oneLineSummary = "Count PASS variants",
        programGroup = VariantProgramGroup.class
)

/**
 * Created by tsato on 1/30/17.
 */
// TODO: maybe this tool should be a variant walker. Or a multiple variant walker.
public class Concordance extends GATKTool {
    @Argument(doc = "evaluation vcf", fullName= "eval", shortName = "E", optional = false)
    protected File eval;

    @Argument(doc = "truth vcf (tool assumes all sites in truth are PASS)", fullName= "truth", shortName = "T", optional = false)
    protected File truth;

    // TO BE IMPLEMENTED
    @Argument(doc = "???", fullName= "confidence", shortName = "C", optional = true)
    protected File confidence_region;

    @Argument(doc = "A table of variants. Prints out for each variant (row) its basic annotations, tumor alt allele fraction, and truth status", fullName="table", shortName = "O", optional = false)
    protected File table;

    @Argument(doc = "A table of summary statistics (true positives, sensitivity, etc.)", fullName="summary", shortName = "S", optional = false)
    protected File summary;

    @Argument(doc = "sample name of the tumor", fullName ="tumorSampleName", shortName = "tumor", optional = false)
    protected String tumorSampleName;

    // TODO: output a vcf of false positives (and negative) sites?
    @Override
    public void traverse() {
        // TODO: take a low confidence region where a false positive there is not counted (masked in dream?)
        // TODO: ideas. for stratifying, output a table of annotaitons (AF, DP, SNP/INDEL, and truth status)
        // @assumes that all variants in the truth vcf are REAL
        // TODO: make match function extendable, such that a user inherits the class to write his/her own evaluator (but get the iteration working first)
        Utils.regularReadableUserFile(eval);
        Utils.regularReadableUserFile(truth);

        try (VariantStatusTableWriter tableWriter = new VariantStatusTableWriter(table);
             SummaryTableWriter summaryWriter = new SummaryTableWriter(summary)) {

            long truePositives = 0;
            long falsePositives = 0;
            long falseNegatives = 0;

            // TODO: note evalSource is probably not a FeatureDataSource. Rather it's a VariantSource
            // in other words this tool should be a variant walker
            final FeatureDataSource<VariantContext> truthSource = new FeatureDataSource<>(truth);
            final FeatureDataSource<VariantContext> evalSource = new FeatureDataSource<>(eval);

            Iterator<VariantContext> truthIterator = truthSource.iterator();
            Iterator<VariantContext> evalIterator = evalSource.iterator();

            VariantContext truthvc = truthIterator.next();
            VariantContext evalvc = evalIterator.next();

            if (! evalvc.getSampleNames().contains(tumorSampleName)){
                // TODO: do I need to close files?
                throw new IllegalArgumentException(String.format("the tumor sample %s does not exit in vcf", tumorSampleName));
            }

            while (true) {
                // the position of two variants match; do genotypes match too?
                if (comparePositions(truthvc, evalvc) == 0) {
                    if (evalvc.isFiltered()) {
                        falseNegatives++;
                        // TODO: ensure that when we write a record, we advance the iterator once - no more, no less
                        // TODO: perhaps the increment and write record should be wrapped in a function
                        tableWriter.writeRecord(new EvalVariantRecord(evalvc, FALSE_NEGATIVE, tumorSampleName));
                    } else if (match(truthvc, evalvc)) {
                        truePositives++;
                        tableWriter.writeRecord(new EvalVariantRecord(evalvc, TRUE_POSITIVE, tumorSampleName));
                    } else {
                        // the eval variant matches truth's position but they don't match
                        // e.g. genotype or alleles doesn't match
                        falsePositives++;
                        falseNegatives++;
                        tableWriter.writeRecord(new EvalVariantRecord(evalvc, FALSE_POSITIVE_AND_FALSE_NEGATIVE, tumorSampleName));
                    }

                    // move on
                    if (!truthIterator.hasNext() || !evalIterator.hasNext()) {
                        break;
                    }

                    truthvc = truthIterator.next();
                    evalvc = evalIterator.next();
                    continue;
                }

                // truth is ahead of eval; eval variant is not in truth
                if (comparePositions(truthvc, evalvc) > 0) {
                    if (evalvc.isNotFiltered()) {
                        falsePositives++;
                        tableWriter.writeRecord(new EvalVariantRecord(evalvc, FALSE_POSITIVE, tumorSampleName));
                    }

                    // move on
                    if (!evalIterator.hasNext()) {
                        break;
                    }
                    evalvc = evalIterator.next();
                    continue;
                }

                // eval leapfrogged truth; we missed a variant
                if (comparePositions(truthvc, evalvc) < 0) {
                    falseNegatives++;

                    // move on
                    if (!truthIterator.hasNext()) {
                        break;
                    }

                    truthvc = truthIterator.next();
                    continue;
                }
            }

            // we exhausted variants in either vcf; process the rest of variants in the other
            if (! truthIterator.hasNext()){
                // exhausted the truth vcf first; the rest of variants in eval are false positives
                // TODO: should be able to use evalIterator.forEachRemaining

                while (evalIterator.hasNext()){
                    evalvc = evalIterator.next();
                    if (evalvc.isNotFiltered()){
                        falsePositives++;
                        tableWriter.writeRecord(new EvalVariantRecord(evalvc, FALSE_POSITIVE, tumorSampleName));
                    }
                }
            } else {
                // exhausted the eval vcf first; the rest of variants in truth are false negatives
                falseNegatives++;
                while (truthIterator.hasNext()){
                    truthvc = truthIterator.next();
                    if (truthvc.isNotFiltered()){
                        falseNegatives++;
                    }
                }
            }

            summaryWriter.writeRecord(new SummaryRecord(truePositives, falsePositives, falseNegatives));

        } catch (IOException exp) {
            throw new UserException.CouldNotCreateOutputFile(table, "Could not create the output table");
        }

    }

    private class SummaryTableWriter extends TableWriter<SummaryRecord> {
        private SummaryTableWriter(final File output) throws IOException {
            super(output, new TableColumnCollection(SummaryRecord.SUMMARY_TABLE_COLUMN_HEADER));
        }

        @Override
        protected void composeLine(final SummaryRecord record, final DataLine dataLine) {
            dataLine.set(SummaryRecord.TRUE_POSITIVE_COLUMN_NAME, record.getTruePositives())
                    .set(SummaryRecord.FALSE_POSITIVE_COLUMN_NAME, record.getFalsePositives())
                    .set(SummaryRecord.FALSE_NEGATIVE_COLUMN_NAME, record.getFalseNegatives())
                    .set(SummaryRecord.SENSITIVITY_COLUMN_NAME, record.getSensitivity())
                    .set(SummaryRecord.PRECISION_COLUMN_NAME, record.getPrecision());
        }

    }

    private class VariantStatusTableWriter extends TableWriter<EvalVariantRecord> {
        private VariantStatusTableWriter(final File output) throws IOException {
            super(output, new TableColumnCollection(EvalVariantRecord.VARIANT_TABLE_COLUMN_HEADERS));
        }

        @Override
        protected void composeLine(final EvalVariantRecord record, final DataLine dataLine) {
            dataLine.set(EvalVariantRecord.CHROMOSOME_COLUMN_NAME, record.getVariantContext().getContig())
                    .set(EvalVariantRecord.START_POSITION_COLUMN_NAME, record.getVariantContext().getStart())
                    .set(EvalVariantRecord.END_POSITION_COLUMN_NAME, record.getVariantContext().getEnd())
                    .set(EvalVariantRecord.REF_ALLELE_COLUMN_NAME, record.getVariantContext().getReference().toString())
                    .set(EvalVariantRecord.ALT_ALLELE_COLUMN_NAME, record.getVariantContext().getAlternateAlleles().toString())
                    .set(EvalVariantRecord.ALLELE_FRACTION_COLUMN_NAME, record.getTumorAlleleFraction()) // TODO: must be able to retrieve the allele fraction from tumor
                    .set(EvalVariantRecord.TRUTH_STATUS_COLUMN_NAME, record.getTruthStatus());
        }
    }



    // Possible statuses
    private static String FALSE_NEGATIVE = "FALSE_NEGATIVE";
    private static String TRUE_POSITIVE = "TRUE_POSITIVE";
    private static String FALSE_POSITIVE_AND_FALSE_NEGATIVE = "FALSE_POSITIVE_AND_FALSE_NEGATIVE";
    private static String FALSE_POSITIVE = "FALSE_POSITIVE";

    // TODO: make variantcontext comparable?
    // TODO: make protected to support subclassingu
    private static int comparePositions(final VariantContext truth, final VariantContext eval){

        if (! truth.getContig().equals(eval.getContig()) ) {
            // TODO: make sure x, y, and mt work as expected
            return truth.getContig().compareTo(eval.getContig());
        } else {
            // same chromosome, compare the start positions
            if (truth.getStart() > eval.getStart()) {
                return 1;
            } else if (truth.getStart() < eval.getStart()){
                return -1;
            } else {
                return 0;
            }
        }
    }

    // TODO: eventually this should be an abstract method to be overridden by a subclass
    // TODO: make protected to support subclassing
    private static boolean match(final VariantContext truth, final VariantContext eval){
        final boolean sameContig = truth.getContig().equals(eval.getContig());
        final boolean sameStartPosition = truth.getStart() == eval.getStart();
        final boolean sameRefAllele = truth.getReference().equals(eval.getReference());
        // TODO: assume single alt allele in truth?
        final boolean containsAltAllele = eval.getAlternateAlleles().contains(truth.getAlternateAllele(0));
        return sameContig && sameStartPosition && sameRefAllele && containsAltAllele;
    }


}
