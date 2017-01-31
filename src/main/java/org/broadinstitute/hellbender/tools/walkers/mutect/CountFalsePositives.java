package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.mutable.MutableLong;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static org.broadinstitute.hellbender.tools.walkers.mutect.FalsePositiveRecord.*;

@CommandLineProgramProperties(
        summary = "Count PASS (false positive) variants in a vcf file for Mutect2 NA12878 normal-normal evaluation",
        oneLineSummary = "Count PASS variants",
        programGroup = VariantProgramGroup.class
)

/**
 * Created by tsato on 12/28/16.
 */
public class CountFalsePositives extends VariantWalker {
    @Argument(
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file"
    )
    private File outputFile = null;
    MutableLong indelFalsePositiveCount = new MutableLong();
    MutableLong snpFalsePositiveCount = new MutableLong();
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
            indelFalsePositiveCount.increment();
        } else {
            snpFalsePositiveCount.increment();
        }
    }

    @Override
    public Object onTraversalSuccess() {
        final List<SimpleInterval> intervals =  intervalArgumentCollection.getIntervals(getReferenceDictionary());
        final long targetTerritory = intervals.stream().mapToInt(i -> i.size()).sum();

        try ( FalsePositiveTableWriter writer = new FalsePositiveTableWriter(outputFile) ) {
            FalsePositiveRecord falsePositiveRecord = new FalsePositiveRecord(id, snpFalsePositiveCount.longValue(), indelFalsePositiveCount.longValue(), targetTerritory);
            writer.writeRecord(falsePositiveRecord);
        } catch (IOException e){
            throw new UserException(String.format("Encountered an IO exception while opening %s", outputFile));
        }

        return "SUCCESS";
    }

    private class FalsePositiveTableWriter extends TableWriter<FalsePositiveRecord> {
        private FalsePositiveTableWriter(final File output) throws IOException {
            super(output, new TableColumnCollection(ID_COLUMN_NAME, SNP_COLUMN_NAME,
                    INDEL_COLUMN_NAME, SNP_FPR_COLUMN_NAME, INDEL_FPR_COLUMN_NAME, TARGET_TERRITORY_COLUMN_NAME));
        }

        @Override
        protected void composeLine(final FalsePositiveRecord falsePositiveRecord, final DataLine dataLine) {
            dataLine.set(ID_COLUMN_NAME, falsePositiveRecord.getId())
                    .set(SNP_COLUMN_NAME, falsePositiveRecord.getSnpFalsePositives())
                    .set(INDEL_COLUMN_NAME, falsePositiveRecord.getIndelFalsePositives())
                    .set(SNP_FPR_COLUMN_NAME, falsePositiveRecord.getSnpFalsePositiveRate())
                    .set(INDEL_FPR_COLUMN_NAME, falsePositiveRecord.getIndelFalsePositiveRate())
                    .set(TARGET_TERRITORY_COLUMN_NAME, falsePositiveRecord.getTargetTerritory());
        }
    }
}
