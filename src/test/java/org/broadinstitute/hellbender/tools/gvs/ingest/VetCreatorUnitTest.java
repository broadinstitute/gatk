package org.broadinstitute.hellbender.tools.gvs.ingest;

import htsjdk.variant.variantcontext.*;
import org.apache.hadoop.fs.FileAlreadyExistsException;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static htsjdk.variant.vcf.VCFConstants.DEPTH_KEY;
import static org.broadinstitute.hellbender.utils.variant.GATKVCFConstants.*;

public class VetCreatorUnitTest {
    private static final long SAMPLE_ID = 100;
    private static final String SAMPLE_NAME = "NA1";
    private static final String PROJECT_ID = "test";
    private static final String DATASET_NAME = "test";
    private File outputDirectory = new File("quickstart/output/");
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



    @Test
    public void testParquetOutputFile() throws IOException {
        String fullPath = String.join(currentRelativePath.toAbsolutePath().toString(), outputDirectory.toString());
        final File parquetOutputFile = new File(fullPath, VET_FILETYPE_PREFIX + tableNumber + PREFIX_SEPARATOR + sampleIdentifierForOutputFileName + ".parquet");

        String expected = String.join(fullPath, "vet_001_parquet.parquet");
        Assert.assertEquals(parquetOutputFile.getAbsoluteFile(), expected);
        Files.deleteIfExists(parquetOutputFile.toPath());
    }

    //@Test(expected = FileAlreadyExistsException.class)
    @Test
    public void testErrorFile() throws IOException {
        VariantContextBuilder builderA =
                new VariantContextBuilder("a","1",10329,10329,
                        Arrays.asList(Allele.REF_C,Allele.ALT_A,Allele.NON_REF_ALLELE));


        Genotype g = new GenotypeBuilder(SAMPLE_NAME)
                .alleles(Arrays.asList(Allele.REF_C, Allele.ALT_A))
                .PL(new int[]{74,0,34,707,390,467})
                .DP(64)
                .GQ(36)
                .AD(new int[]{22,42,0})
                .attribute(STRAND_BIAS_BY_SAMPLE_KEY, "1,21,6,50")
                .make();

        builderA.attribute(AS_RAW_RMS_MAPPING_QUALITY_KEY,"29707.00|39366.00|2405.00")
                .attribute(AS_RAW_MAP_QUAL_RANK_SUM_KEY,"|-0.2,1|-2.5,1")
                .attribute(RAW_QUAL_APPROX_KEY,"74")
                .attribute(AS_RAW_QUAL_APPROX_KEY,"|74|0")
                .attribute(AS_RAW_READ_POS_RANK_SUM_KEY,"|2.4,1|1.5,1")
                .attribute(AS_SB_TABLE_KEY,"1,21|3,39|3,11")
                .attribute(AS_VARIANT_DEPTH_KEY,"22|42|0")
                .genotypes(Arrays.asList(g));

        VariantContext vc = builderA.make();

        final File parquetOutputFile = new File(outputDirectory, VET_FILETYPE_PREFIX + tableNumber + PREFIX_SEPARATOR + sampleIdentifierForOutputFileName + ".parquet");
        // Path tempFile = Files.createTempFile(parquetOutputFile.getAbsolutePath());
       // Files.createTempFile("vet_001_parquet", ".parquet");
        String sampleIdentifierForOutputFileName = "bleh";

        VetCreator vetCreator = new VetCreator(parquetOutputFile.getName(), SAMPLE_ID, tableNumber, outputDirectory, outputType, PROJECT_ID, DATASET_NAME, true, false, PARQUET_SCHEMA);
        List<String> row = vetCreator.createRow(10329,vc, SAMPLE_NAME);

        Assert.assertEquals("/ by zero", row.get(0));
        Files.deleteIfExists(parquetOutputFile.toPath());

       // Assert.assertEquals("/ by zero", );
    }

}
