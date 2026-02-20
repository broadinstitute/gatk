package org.broadinstitute.hellbender.tools.gvs.ingest;

import htsjdk.variant.variantcontext.*;
import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.MessageTypeParser;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.IngestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;

import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.*;

public class VetCreatorUnitTest {
    private static final long SAMPLE_ID = 100;
    private static final String SAMPLE_NAME = "NA1";
    private static final String PROJECT_ID = "test";
    private static final String DATASET_NAME = "test";
    private final File outputDirectory = new File("quickstart/output/");
    Path currentRelativePath = Paths.get("");
    private final CommonCode.OutputType outputType = CommonCode.OutputType.PARQUET;
    private static final String VET_FILETYPE_PREFIX = "vet_"; // should this live somewhere else--check out IngestConstants for instance--why is that a tsv?!?!
    String PREFIX_SEPARATOR = "_"; // should this live somewhere else?

    int sampleTableNumber = IngestUtils.getTableNumber(SAMPLE_ID, IngestConstants.partitionPerTable);
    String tableNumber = String.format("%03d", sampleTableNumber);
    String sampleIdentifierForOutputFileName = "parquet";
    public final MessageType PARQUET_SCHEMA = MessageTypeParser // do we want this in a utils file? or as part of a method?
            .parseMessageType("""
            message VariantRow {
            	required int64 sample_id;
            	required int64 location;
            	required binary ref (UTF8);
            	required binary alt (UTF8);
            	optional binary AS_RAW_MQ (UTF8);
            	optional binary AS_RAW_MQRankSum (UTF8);
            	optional binary AS_QUALapprox (UTF8);
            	optional binary AS_RAW_ReadPosRankSum (UTF8);
            	optional binary AS_SB_TABLE (UTF8);
            	optional binary AS_VarDP (UTF8);
            	required binary call_GT (UTF8);
            	optional binary call_AD (UTF8);
            	optional binary call_DP (UTF8);
            	required int64 call_GQ;
            	optional binary call_PGT (UTF8);
            	optional binary call_PID (UTF8);
            	optional binary call_PS (UTF8);
            	optional binary call_PL (UTF8);
            }
            """);



    @Test(enabled = false)
    public void testParquetOutputFile() throws IOException {
        String fullPath = String.join(currentRelativePath.toAbsolutePath().toString(), outputDirectory.toString());
        final File parquetOutputFile = new File(fullPath, VET_FILETYPE_PREFIX + tableNumber + PREFIX_SEPARATOR + sampleIdentifierForOutputFileName + ".parquet");

        String expected = String.join(fullPath, "vet_001_parquet.parquet");
        Assert.assertEquals(parquetOutputFile.getAbsoluteFile(), expected);
        Files.deleteIfExists(parquetOutputFile.toPath());
    }

}
