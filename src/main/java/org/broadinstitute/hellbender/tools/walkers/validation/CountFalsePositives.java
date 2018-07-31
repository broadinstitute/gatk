package org.broadinstitute.hellbender.tools.walkers.validation;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.io.FilenameUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.List;
/**
 * Count variants which were not filtered in a VCF.
 * This is used to count false positives for Mutect2 normal-normal evaluation.
 *
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk --java-options "-Xmx4g" CountFalsePositives \
 *   -V filtered.vcf.gz \
 *   -R ref_fasta.fa \
 *   -O false-positives.txt
 * </pre>
 *
 * Created by Takuto Sato on 12/28/16.
 */

@CommandLineProgramProperties(
        summary = CountFalsePositives.USAGE_SUMMARY,
        oneLineSummary = CountFalsePositives.USAGE_ONE_LINE_SUMMARY,
        programGroup = VariantEvaluationProgramGroup.class
)

@DocumentedFeature
@BetaFeature
public class CountFalsePositives extends VariantWalker {
    static final String USAGE_ONE_LINE_SUMMARY = "Count PASS variants";
    static final String USAGE_SUMMARY = "Count PASS (false positive) variants in a vcf file for Mutect2 NA12878 normal-normal evaluation";

    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file"
    )
    private File outputFile = null;
    long indelFalsePositiveCount = 0;
    long snpFalsePositiveCount = 0;
    String id;

    // TODO: eventually use tumor and normal sample names instead of the file name. To do so we must extract them from the vcf, which I don't know how.
    @Override
    public void onTraversalStart() {
        // TODO: ideally we would identify the tumor sample name and normal sample name from the vcf header but
        // there doesn't seem to be a function that let us query a sample name by its ID (e.g. Tumor/Normal) in VCFHeader class.
        id = FilenameUtils.getBaseName(drivingVariantFile);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (variant.isFiltered()) {
            return;
        }

        if (variant.isIndel()) {
            indelFalsePositiveCount++;
        } else {
            snpFalsePositiveCount++;
        }
    }

    @Override
    public Object onTraversalSuccess() {
        final List<SimpleInterval> intervals =  intervalArgumentCollection.getIntervals(getReferenceDictionary());
        final long targetTerritory = intervals.stream().mapToLong(i -> i.size()).sum();

        try ( FalsePositiveTableWriter writer = new FalsePositiveTableWriter(outputFile) ) {
            FalsePositiveRecord falsePositiveRecord = new FalsePositiveRecord(id, snpFalsePositiveCount, indelFalsePositiveCount, targetTerritory);
            writer.writeRecord(falsePositiveRecord);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while opening %s", outputFile));
        }

        return "SUCCESS";
    }

    private class FalsePositiveTableWriter extends TableWriter<FalsePositiveRecord> {
        private FalsePositiveTableWriter(final File output) throws IOException {
            super(output, new TableColumnCollection(FalsePositiveRecord.ID_COLUMN_NAME, FalsePositiveRecord.SNP_COLUMN_NAME,
                    FalsePositiveRecord.INDEL_COLUMN_NAME, FalsePositiveRecord.SNP_FPR_COLUMN_NAME, FalsePositiveRecord.INDEL_FPR_COLUMN_NAME, FalsePositiveRecord.TARGET_TERRITORY_COLUMN_NAME));
        }

        @Override
        protected void composeLine(final FalsePositiveRecord falsePositiveRecord, final DataLine dataLine) {
            dataLine.set(FalsePositiveRecord.ID_COLUMN_NAME, falsePositiveRecord.getId())
                    .set(FalsePositiveRecord.SNP_COLUMN_NAME, falsePositiveRecord.getSnpFalsePositives())
                    .set(FalsePositiveRecord.INDEL_COLUMN_NAME, falsePositiveRecord.getIndelFalsePositives())
                    .set(FalsePositiveRecord.SNP_FPR_COLUMN_NAME, falsePositiveRecord.getSnpFalsePositiveRate())
                    .set(FalsePositiveRecord.INDEL_FPR_COLUMN_NAME, falsePositiveRecord.getIndelFalsePositiveRate())
                    .set(FalsePositiveRecord.TARGET_TERRITORY_COLUMN_NAME, falsePositiveRecord.getTargetTerritory());
        }
    }
}
