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
    private final File dreamDir =  new File(toolsTestDir, "mutect/dream");

    // expected numbers are computed as follows:
    // SNP: grep PASS dream3-chr20.vcf | awk 'length($4) == length($5) { print $0 }' | wc -l
    // INDEL: grep PASS dream3-chr20.vcf | awk 'length($4) != length($5) { print $0 }' | wc -l
    private final int expectedSnpFPs = 173;
    private final int expectedIndelFPs = 276;

    @Test
    public void testSimple() throws Exception {
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

        Assert.assertEquals(record.getSnpFalsePositives(), expectedSnpFPs);
        Assert.assertEquals(record.getIndelFalsePositives(), expectedIndelFPs);
    }

    @Test
    public void testWGS() throws Exception {
        final File output = createTempFile("output", ".txt");
        final String[] args = {
                "-V", dreamDir + "/vcfs/dream3-chr20.vcf",
                "-R", b37_reference_20_21,
                "-L", wgsIntervalFile,
                "-O", output.toString()
        };
        runCommandLine(args);

        FalsePositiveRecordReader reader = new FalsePositiveRecordReader(output);
        FalsePositiveRecord record = reader.readRecord();

        final long wgsIntervalSize = 2832598739L;
        Assert.assertEquals(record.getTargetTerritory(), wgsIntervalSize);
        Assert.assertEquals(record.getSnpFalsePositiveRate(), (double) expectedSnpFPs/wgsIntervalSize * 1e6);
        Assert.assertEquals(record.getIndelFalsePositiveRate(), (double) expectedIndelFPs/wgsIntervalSize * 1e6);
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