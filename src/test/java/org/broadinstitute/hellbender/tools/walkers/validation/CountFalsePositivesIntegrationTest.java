package org.broadinstitute.hellbender.tools.walkers.validation;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.validation.FalsePositiveRecord;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.testng.Assert;
import org.testng.annotations.Test;
import java.io.File;
import java.io.IOException;

import static org.broadinstitute.hellbender.tools.walkers.validation.FalsePositiveRecord.*;

/**
 * Created by tsato on 12/28/16.
 */
public class CountFalsePositivesIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testSimple() throws Exception {
        final File dreamDir =  new File(publicTestDir, "org/broadinstitute/hellbender/tools/mutect/dream");
        final File output = createTempFile("output", ".txt");
        final String[] args = {
                "-V", dreamDir + "/vcfs/dream3-chr20.vcf",
                "-R", b37_reference_20_21,
                "-L", dreamDir + "/dream-chr20.interval_list",
                "-O", output.toString()
        };
        runCommandLine(args);

        FalsePositiveRecordReader reader = new FalsePositiveRecordReader(output);
        FalsePositiveRecord record = reader.readRecord();

        // check that the tool agrees with the false positive counts using AWK
        // SNP: grep PASS dream3-chr20.vcf | awk 'length($4) == length($5) { print $0 }' | wc -l
        // INDEL: grep PASS dream3-chr20.vcf | awk 'length($4) != length($5) { print $0 }' | wc -l
        Assert.assertEquals(record.getSnpFalsePositives(), 173);
        Assert.assertEquals(record.getIndelFalsePositives(), 276);
    }

    private class FalsePositiveRecordReader extends TableReader<FalsePositiveRecord> {
        private FalsePositiveRecordReader(final File falsePositiveTable) throws IOException {
            super(falsePositiveTable);
        }

        @Override
        protected FalsePositiveRecord createRecord(final DataLine dataLine) {
            final String id = dataLine.get(ID_COLUMN_NAME);
            final long snpFalsePositives = Long.parseLong(dataLine.get(SNP_COLUMN_NAME));
            final long indelFalsePositives = Long.parseLong(dataLine.get(INDEL_COLUMN_NAME));
            final long targetTerritory = Long.parseLong(dataLine.get(TARGET_TERRITORY_COLUMN_NAME));

            return new FalsePositiveRecord(id, snpFalsePositives, indelFalsePositives, targetTerritory);
        }
    }
}