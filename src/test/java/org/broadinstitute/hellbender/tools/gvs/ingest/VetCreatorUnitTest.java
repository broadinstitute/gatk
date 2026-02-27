package org.broadinstitute.hellbender.tools.gvs.ingest;

import htsjdk.variant.variantcontext.*;
import org.apache.hadoop.fs.FileAlreadyExistsException;
import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.MessageTypeParser;
import org.broadinstitute.hellbender.exceptions.UserException;
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

public class VetCreatorUnitTest {
    private static final long SAMPLE_ID = 100;
    private static final String PROJECT_ID = "test";
    private static final String DATASET_NAME = "test";
    private final File outputDirectory = new File("quickstart/output/");
    Path currentRelativePath = Paths.get("");
    private final CommonCode.OutputType outputType = CommonCode.OutputType.PARQUET;

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



    @Test
    public void testParquetOutputFileNaming() throws IOException {
        String fullPath = Paths.get(currentRelativePath.toAbsolutePath().toString(), outputDirectory.toString()).toString();
        String filename = VetCreator.getOutputFileName(tableNumber, 27L, sampleIdentifierForOutputFileName, outputType);
        final File parquetOutputFile = new File(fullPath, filename);

        String expected = Paths.get(fullPath, "vet_001_27_parquet.parquet").toString();
        Assert.assertEquals(parquetOutputFile.getAbsolutePath(), expected);
        Files.deleteIfExists(parquetOutputFile.toPath());
    }

    @Test
    public void testErrorFile() throws IOException {
        // This test verifies that FileAlreadyExistsException is properly thrown when
        // attempting to create a VetCreator with a file that already exists

        // Ensure the output directory exists
        Files.createDirectories(outputDirectory.toPath());

        // Calculate the actual filename that will be generated
        String filename = VetCreator.getOutputFileName(tableNumber, SAMPLE_ID, sampleIdentifierForOutputFileName, outputType);
        final File parquetOutputFile = new File(outputDirectory, filename);

        // Clean up any existing file first
        Files.deleteIfExists(parquetOutputFile.toPath());

        // Create the file to simulate it already existing
        Files.createFile(parquetOutputFile.toPath());

        try {
            // This should throw a FileAlreadyExistsException
            new VetCreator(
                sampleIdentifierForOutputFileName,
                SAMPLE_ID,
                tableNumber,
                outputDirectory,
                outputType,
                PROJECT_ID,
                DATASET_NAME,
                true,
                false,
                PARQUET_SCHEMA
            );

            // If we get here, the test failed - no exception was thrown
            Assert.fail("Expected FileAlreadyExistsException to be thrown (and wrapped) when the file already exists");
        } catch (UserException e) {
            Assert.assertTrue(e.getCause() instanceof FileAlreadyExistsException, e.getCause().toString());
        } catch (Exception e) {
            Assert.fail(e.getMessage());
        } finally {
            // Clean up the test file
            Files.deleteIfExists(parquetOutputFile.toPath());
        }
    }
}

