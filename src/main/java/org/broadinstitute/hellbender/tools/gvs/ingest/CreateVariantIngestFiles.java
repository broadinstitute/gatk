package org.broadinstitute.hellbender.tools.gvs.ingest;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.parquet.schema.MessageTypeParser;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.DeprecatedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.common.*;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.IntervalUtils;

import static org.apache.parquet.schema.PrimitiveType.PrimitiveTypeName.BINARY;
import static org.apache.parquet.schema.PrimitiveType.PrimitiveTypeName.INT64;
import static org.apache.parquet.schema.Type.Repetition.OPTIONAL;
import static org.apache.parquet.schema.Type.Repetition.REQUIRED;

import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.PrimitiveType;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Ingest variant walker
 */
@CommandLineProgramProperties(
        summary = "Exome and Genome Ingest tool for the Genomic Variant Store",
        oneLineSummary = "Ingest tool for GVS",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = true
)
public final class CreateVariantIngestFiles extends VariantWalker {
    static final Logger logger = LogManager.getLogger(CreateVariantIngestFiles.class);

    private RefCreator refCreator;
    private VetCreator vetCreator;
    private VcfHeaderLineScratchCreator vcfHeaderLineScratchCreator;
    private SamplePloidyCreator samplePloidyCreator;
    private LoadStatus loadStatus;

    private final Map<String, Boolean> allLineHeaders = new HashMap<>();

    private GenomeLocSortedSet intervalArgumentGenomeLocSortedSet;

    private Long sampleId;

    // Inside the parent directory, a directory for each chromosome will be created, with a vet directory in each one.
    // Each vet directory will hold all the vet TSVs for each sample.
    // A sample_info directory will be created, with a sample_info tsv for each sample.

    @Argument(fullName = "ref-block-gq-to-ignore",
            shortName = "IG",
            doc = "Ref Block GQ band to ignore, bands of 10 e.g 0-9 get combined to 0, 20-29 get combined to 20",
            optional = true)
    public GQStateEnum gqStateToIgnore = GQStateEnum.SIXTY;

    @DeprecatedFeature(detail="Argument of no foreseeable use")
    @Argument(fullName = "ignore-above-gq-threshold",
            shortName = "GTIG",
            doc = "in addition to dropping the gq block specified by ref-block-gq-to-ignore, also drop higher gq blocks",
            optional = true)
    public boolean dropAboveGqThreshold = false;

    @Argument(fullName = "enable-reference-ranges",
            shortName = "rr",
            doc = "write reference ranges data",
            optional = true)
    public boolean enableReferenceRanges = true;

    @Argument(fullName = "enable-vet",
            shortName = "ev",
            doc = "write vet data",
            optional = true)
    public boolean enableVet = true;

    @Argument(
            fullName = "enable-vcf-headers",
            doc = "write VCF header lines",
            optional = true
    )
    public boolean enableVCFHeaders = false;

    @Argument(fullName = "sample-name-mapping",
            shortName = "SNM",
            doc = "Sample name to sample id mapping. This must be provided if gvs-sample-id is not",
            optional = true)
    public File sampleMap;

    @Argument(fullName = "gvs-sample-id",
            shortName = "GVSID",
            doc = "GVS identifier for the sample. Can be looked up by external-sample-name in the mapping file if provided.",
            optional = true)
    public Long sampleIdParam;

    @Argument(fullName = "sample-name",
            shortName = "SN",
            doc = "The external sample name used for the sample. If this parameter is not provided, the sample name in the gvcf file will be used. If providing a sample-name-mapping file, this is the name that must be mapped to the id.",
            optional = true)
    public String sampleNameParam;

    @Argument(fullName = "load-status-table-name",
            doc = "Table to insert the sample_id when a sample has been successfully loaded",
            optional = true)
    public String loadStatusTableName = "sample_load_status";

    @Argument(fullName = "output-type",
            shortName = "ot",
            doc = "[Experimental] Output file format: TSV, ORC, PARQUET or BQ [default=TSV].",
            optional = true)
    public CommonCode.OutputType outputType = CommonCode.OutputType.TSV;

    @Argument(
            fullName = "ref-version",
            doc = "Remove this option!!!! only for ease of testing. Valid options are 37 or 38",
            optional = true)
    public String refVersion = "37";

    @Argument(
            fullName = "output-directory",
            doc = "directory for output tsv files",
            optional = true)
    public File outputDir = new File(".");

    @Argument(
            fullName = "project-id",
            doc = "ID of the Google Cloud project where the dataset for vet tables exists",
            optional = true
    )
    public String projectID = null;

    @Argument(
            fullName = "dataset-name",
            doc = "Name of the dataset to update vet tables",
            optional = true
    )
    public String datasetName = null;

    @Argument(
            fullName = "force-loading-from-non-allele-specific",
            doc = "Even if there are allele-specific (AS) annotations, use backwards compatibility mode",
            optional = true
    )
    public boolean forceLoadingFromNonAlleleSpecific = false;

    @Argument(
            fullName = "skip-loading-vqsr-fields",
            doc = "Do not load data for fields specific to VQSR",
            optional = true
    )
    public boolean skipLoadingVqsrFields = false;

    @Argument(
            fullName = "use-compressed-refs",
            doc = "Store the ref_ranges data in a compressed format. Saves about 40% on storage",
            optional = true
    )
    public boolean storeCompressedReferences = false;

    private boolean shouldWriteReferencesLoadedStatusRow = false;

    private boolean shouldWriteVariantsLoadedStatusRow = false;

    private boolean shouldWriteVCFHeadersLoadedStatusRow = false;

    private final Set<GQStateEnum> gqStatesToIgnore = new HashSet<>();

    public final MessageType variantRowSchema = MessageTypeParser
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

    /* Not yet outputting ref_ranges rows */
    public final MessageType refRangesRowSchema = MessageTypeParser
            .parseMessageType("""
            message VariantRow {
            	required int64 sample_id;
            	required int64 location;
            	required int64 length;
            	required binary state (UTF8);
            }
            """);


    // getGenotypes() returns list of lists for all samples at variant
    // assuming one sample per gvcf, getGenotype(0) retrieves GT for sample at index 0
    public static boolean isNoCall(VariantContext variant) {
        return variant.getGenotype(0).isNoCall();
    }

    @Override
    public boolean requiresIntervals() {
        return true; // TODO -- do I need to check the boolean flag on this?
    }

    private String getInputFileName() {
        // this returns the full file name including extensions
        String[] pathParts = drivingVariantFile.toString().split("/");
        return pathParts[pathParts.length - 1];
    }

    @Override
    public void onTraversalStart() {
        if (!enableVet && !enableReferenceRanges && !enableVCFHeaders) {
            throw new RuntimeIOException("Invalid invocation: variants, references, and VCF header writing all set to false");
        }

        // Set up output directory:
        if (!outputDir.exists() && !outputDir.mkdir()) {
            throw new RuntimeIOException("Unable to create directory: " + outputDir.getAbsolutePath());
        }

        // Set reference version -- TODO remove this in the future, also, can we get ref version from the header?
        ChromosomeEnum.setRefVersion(refVersion);

        String sampleName = sampleNameParam == null ? IngestUtils.getSampleName(getHeaderForVariants()) : sampleNameParam;
        if (sampleIdParam == null && sampleMap == null) {
            throw new IllegalArgumentException("One of sample-id or sample-name-mapping must be specified");
        }
        if (sampleIdParam != null) {
            sampleId = sampleIdParam;
        } else {
            sampleId = IngestUtils.getSampleId(sampleName, sampleMap);
        }

        // TODO when we pass in the full file path or gvs_id as an input arg, use path here instead or gvs_id in addition
        // use input gvcf file name instead of sample name in output filenames
        String sampleIdentifierForOutputFileName = getInputFileName();

        // Mod the sample directories
        int sampleTableNumber = IngestUtils.getTableNumber(sampleId, IngestConstants.partitionPerTable);
        String tableNumber = String.format("%03d", sampleTableNumber);

        boolean refRangesRowsExist = false;
        boolean vetRowsExist = false;
        boolean vcfHeaderRowsExist = false;

//        List<ColumnDescriptor> cols = variantRowSchema.getColumns();
//        for (int i = 0; i < cols.size(); ++i) {
//            ColumnDescriptor col = cols.get(i);
//            logger.info("col.getPath()");
//            logger.info(col.getPath());
//            logger.info("col.getPrimitiveType()");
//            logger.info(col.getPrimitiveType());
//            logger.info("col.getPrimitiveType().getName()");
//            logger.info(col.getPrimitiveType().getName());
//
//        }


        /* TEST CODE BELOW WORKS
        java.nio.file.Path currentRelativePath = Paths.get("");
        String s = currentRelativePath.toAbsolutePath().toString();
        System.out.println("Current absolute path is: " + s);
        File testParquetOutputFile = new File(s+ "/testFile.parquet");
        System.out.println("Test file will be written to: " + testParquetOutputFile.getAbsolutePath());
        try {
            GvsVariantParquetFileWriter testWriter = new GvsVariantParquetFileWriter(new Path(testParquetOutputFile.toURI()), variantRowSchema, false, CompressionCodecName.SNAPPY);
            // feed the file some test data
            JSONArray testData = new JSONArray("""
                    [{"AS_QUALapprox":"25","AS_RAW_MQ":"0|0","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"25","alt":"C","call_AD":null,"call_GQ":"3","call_GT":"1/1","call_PGT":"0|1","call_PID":"26885022_G_C","call_PL":"25,3,0,25,3,25","location":"17000026885022","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"45","AS_RAW_MQ":"0|0","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"45","alt":"T","call_AD":null,"call_GQ":"3","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"45,3,0,45,3,45","location":"15000101914931","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"0","AS_RAW_MQ":"0|43200","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"1","alt":"C","call_AD":"0,0","call_GQ":"0","call_GT":"0/1","call_PGT":"0|1","call_PID":"73839326_CAAA_C","call_PL":"1,0,0,1,0,1","location":"16000073839326","ref":"CAAA","sample_id":"7"},
                    {"AS_QUALapprox":"0","AS_RAW_MQ":"0|6894","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"1","alt":"A","call_AD":"0,0","call_GQ":"0","call_GT":"0/1","call_PGT":"0|1","call_PID":"66798833_ATGG_A","call_PL":"1,0,0,1,0,1","location":"17000066798833","ref":"ATGG","sample_id":"7"},
                    {"AS_QUALapprox":"0","AS_RAW_MQ":"0|18000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"1","alt":"C","call_AD":"0,0","call_GQ":"0","call_GT":"0/1","call_PGT":null,"call_PID":null,"call_PL":"1,0,0,1,0,1","location":"16000080965202","ref":"CATATAT","sample_id":"7"},
                    {"AS_QUALapprox":"0","AS_RAW_MQ":"0|43200","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"1","alt":"C","call_AD":"0,0","call_GQ":"0","call_GT":"0/1","call_PGT":null,"call_PID":null,"call_PL":"1,0,0,1,0,1","location":"16000003489112","ref":"CTTTT","sample_id":"7"},
                    {"AS_QUALapprox":"0","AS_RAW_MQ":"0|45450","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"1","alt":"G","call_AD":"0,0","call_GQ":"0","call_GT":"0/1","call_PGT":null,"call_PID":null,"call_PL":"1,0,0,1,0,1","location":"16000046773931","ref":"GC","sample_id":"7"},
                    {"AS_QUALapprox":"0","AS_RAW_MQ":"0|28800","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"1","alt":"C","call_AD":"0,0","call_GQ":"0","call_GT":"0/1","call_PGT":null,"call_PID":null,"call_PL":"1,0,0,1,0,1","location":"15000070796035","ref":"CAAAAA","sample_id":"7"},
                    {"AS_QUALapprox":"0","AS_RAW_MQ":"0|36110","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"1","alt":"A","call_AD":"0,0","call_GQ":"0","call_GT":"0/1","call_PGT":"0|1","call_PID":"72146971_AG_A","call_PL":"1,0,0,1,1,1","location":"17000072146974","ref":"AGGAGG","sample_id":"7"},
                    {"AS_QUALapprox":"0","AS_RAW_MQ":"0|36110","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"1","alt":"A","call_AD":"0,0","call_GQ":"0","call_GT":"0/1","call_PGT":"0|1","call_PID":"72146971_AG_A","call_PL":"1,0,0,1,1,1","location":"17000072146971","ref":"AG","sample_id":"7"},
                    {"AS_QUALapprox":"0","AS_RAW_MQ":"0|54000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"1","alt":"CGTCTCTCTG","call_AD":"0,0","call_GQ":"0","call_GT":"0/1","call_PGT":"0|1","call_PID":"87178619_CG_C","call_PL":"1,0,0,2,2,4","location":"16000087178644","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"0","AS_RAW_MQ":"0|55820","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"1","alt":"T","call_AD":"0,0","call_GQ":"0","call_GT":"0/1","call_PGT":"0|1","call_PID":"13467380_TGGG_T","call_PL":"1,0,0,3,2,5","location":"17000013467408","ref":"A","sample_id":"7"},
                    {"AS_QUALapprox":"1","AS_RAW_MQ":"0|7200","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"1","alt":"T","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"1,1,0,1,1,1","location":"16000059226361","ref":"TATATATAATATAC","sample_id":"7"},
                    {"AS_QUALapprox":"2","AS_RAW_MQ":"0|142210","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,3","AS_VarDP":"0|0","QUALapprox":"2","alt":"CAA","call_AD":"0,0","call_GQ":"0","call_GT":"0/1","call_PGT":null,"call_PID":null,"call_PL":"2,0,0,1,1,1","location":"16000024911166","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"2","AS_RAW_MQ":"0|28289","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"2","alt":"T","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"2,1,0,1,1,1","location":"15000100466812","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"2","AS_RAW_MQ":"0|44336","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"2","alt":"T","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"2,1,0,2,1,1","location":"16000062126446","ref":"TG","sample_id":"7"},
                    {"AS_QUALapprox":"2","AS_RAW_MQ":"0|43825","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"2","alt":"C","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":"0|1","call_PID":"17066723_CTTTTTTT_C","call_PL":"2,1,0,2,1,2","location":"17000017066723","ref":"CTTTTTTT","sample_id":"7"},
                    {"AS_QUALapprox":"2","AS_RAW_MQ":"0|47529","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"2","alt":"G","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":"0|1","call_PID":"74109516_AAT_A","call_PL":"2,1,0,2,2,3","location":"16000074109606","ref":"T","sample_id":"7"},
                    {"AS_QUALapprox":"2","AS_RAW_MQ":"0|54000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"2","alt":"C","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":"0|1","call_PID":"87178619_CG_C","call_PL":"2,1,0,3,2,4","location":"16000087178640","ref":"CT","sample_id":"7"},
                    {"AS_QUALapprox":"3","AS_RAW_MQ":"0|4988","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"3","alt":"A","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"3,1,0,3,1,3","location":"18000012870365","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"3","AS_RAW_MQ":"0|25929","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"3","alt":"AAC","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":"0|1","call_PID":"74109516_AAT_A","call_PL":"3,1,0,3,1,3","location":"16000074109638","ref":"A","sample_id":"7"},
                    {"AS_QUALapprox":"3","AS_RAW_MQ":"0|29529","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"3","alt":"A","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":"0|1","call_PID":"74109516_AAT_A","call_PL":"3,1,0,3,1,3","location":"16000074109633","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"3","AS_RAW_MQ":"0|43929","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"3","alt":"TTA","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":"0|1","call_PID":"74109516_AAT_A","call_PL":"3,1,0,3,2,4","location":"16000074109622","ref":"T","sample_id":"7"},
                    {"AS_QUALapprox":"3","AS_RAW_MQ":"0|56334","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"3","alt":"A","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"3,1,0,4,3,5","location":"17000013467441","ref":"T","sample_id":"7"},
                    {"AS_QUALapprox":"3","AS_RAW_MQ":"0|25929","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"3","alt":"AAT","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":"0|1","call_PID":"74109516_AAT_A","call_PL":"3,2,0,3,2,3","location":"16000074109643","ref":"A","sample_id":"7"},
                    {"AS_QUALapprox":"3","AS_RAW_MQ":"0|57600","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"3","alt":"GCC","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":"0|1","call_PID":"87178619_CG_C","call_PL":"3,2,0,4,2,4","location":"16000087178628","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"4","AS_RAW_MQ":"0|72000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"4","alt":"CTTTCTCCCTCCTT","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"4,2,0,5,2,5","location":"16000087178618","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"4","AS_RAW_MQ":"0|61200","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"4","alt":"C","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":"0|1","call_PID":"87178619_CG_C","call_PL":"4,2,0,5,3,5","location":"16000087178619","ref":"CG","sample_id":"7"},
                    {"AS_QUALapprox":"5","AS_RAW_MQ":"0|55521","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"5","alt":"T","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"5,2,0,5,2,5","location":"15000079059222","ref":"TCCTTC","sample_id":"7"},
                    {"AS_QUALapprox":"6","AS_RAW_MQ":"0|55820","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"6","alt":"TTTCCC","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"6,3,0,5,2,4","location":"17000013467376","ref":"T","sample_id":"7"},
                    {"AS_QUALapprox":"6","AS_RAW_MQ":"0|54000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"6","alt":"A","call_AD":"0,0","call_GQ":"3","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"6,3,0,6,3,6","location":"17000009630990","ref":"AAAAGAAAG","sample_id":"7"},
                    {"AS_QUALapprox":"6","AS_RAW_MQ":"0|56947","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"6","alt":"T","call_AD":"0,0","call_GQ":"3","call_GT":"1/1","call_PGT":"0|1","call_PID":"71253570_TATATAC_T","call_PL":"6,3,0,6,4,7","location":"16000071253584","ref":"TATATAC","sample_id":"7"},
                    {"AS_QUALapprox":"6","AS_RAW_MQ":"0|53347","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"6","alt":"T","call_AD":"0,0","call_GQ":"3","call_GT":"1/1","call_PGT":"0|1","call_PID":"71253570_TATATAC_T","call_PL":"6,3,0,7,4,7","location":"16000071253570","ref":"TATATAC","sample_id":"7"},
                    {"AS_QUALapprox":"9","AS_RAW_MQ":"0|25200","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"9","alt":"TA","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"9,2,0,9,3,9","location":"17000067472332","ref":"T","sample_id":"7"},
                    {"AS_QUALapprox":"16","AS_RAW_MQ":"0|33129","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"16","alt":"CTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"16,3,0,3,0,0","location":"17000030026961","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"17","AS_RAW_MQ":"0|32400","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"17","alt":"CTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"17,3,0,3,0,0","location":"16000074081756","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"17","AS_RAW_MQ":"0|28800","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"17","alt":"CT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"17,3,0,3,0,0","location":"16000077464450","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"17","AS_RAW_MQ":"0|18000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"17","alt":"CAAAAAA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"17,3,0,3,0,0","location":"16000003276234","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"17","AS_RAW_MQ":"0|39600","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"17","alt":"CTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"17,5,0,5,0,0","location":"16000030379640","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"18","AS_RAW_MQ":"0|35764","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"18","alt":"CA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"18,3,0,3,0,0","location":"17000027579878","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"20","AS_RAW_MQ":"0|32400","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"20","alt":"C","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"20,3,0,3,0,0","location":"15000061980183","ref":"CAA","sample_id":"7"},
                    {"AS_QUALapprox":"23","AS_RAW_MQ":"0|46800","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"23","alt":"CTTTTTTTTTTTTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"23,1,0,1,0,0","location":"16000024669941","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"26","AS_RAW_MQ":"0|103504","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,2","AS_VarDP":"0|0","QUALapprox":"26","alt":"CAA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"26,4,0,4,0,0","location":"17000039382092","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"27","AS_RAW_MQ":"0|43200","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"27","alt":"CTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"27,2,0,2,0,0","location":"16000017514797","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"31","AS_RAW_MQ":"0|36000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"31","alt":"T","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"31,2,0,3,1,1","location":"16000087178667","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"34","AS_RAW_MQ":"0|28800","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"34","alt":"CA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"34,3,0,3,0,0","location":"17000028822802","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"36","AS_RAW_MQ":"0|32400","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"36","alt":"CTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"36,3,0,3,0,0","location":"15000061236834","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"40","AS_RAW_MQ":"0|28800","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,2","AS_VarDP":"0|0","QUALapprox":"40","alt":"GA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"40,6,0,6,0,0","location":"15000082302715","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"42","AS_RAW_MQ":"0|32400","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"42","alt":"CTTTTTTTTTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"42,2,0,2,0,0","location":"16000030913975","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"44","AS_RAW_MQ":"0|18000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"44","alt":"CAAAAAAAAAAAAAA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"44,2,0,2,0,0","location":"17000013560289","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"44","AS_RAW_MQ":"0|18000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"44","alt":"CAAAAAAAAAAA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"44,3,0,3,0,0","location":"17000041879647","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"44","AS_RAW_MQ":"0|43200","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"44","alt":"CTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"44,3,0,3,0,0","location":"16000087731046","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"44","AS_RAW_MQ":"0|28800","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"44","alt":"CTTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"44,3,0,3,0,0","location":"16000057292424","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"45","AS_RAW_MQ":"0|4129","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"45","alt":"AATAT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"45,3,0,3,0,0","location":"17000011812026","ref":"A","sample_id":"7"},
                    {"AS_QUALapprox":"51","AS_RAW_MQ":"0|39600","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"51","alt":"GTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"51,5,0,7,0,1","location":"16000075799147","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"52","AS_RAW_MQ":"0|10800","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"52","alt":"GTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"52,8,0,8,0,0","location":"17000076388639","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"54","AS_RAW_MQ":"0|21600","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"54","alt":"CAAAAAAAAAAAA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"54,5,0,5,0,0","location":"16000074958592","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"62","AS_RAW_MQ":"0|32400","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"62","alt":"CAAAAAAAAA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"62,5,0,5,0,0","location":"17000027903859","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"63","AS_RAW_MQ":"0|10800","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"63","alt":"ATATATATTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"63,6,0,6,0,0","location":"17000035720573","ref":"A","sample_id":"7"},
                    {"AS_QUALapprox":"71","AS_RAW_MQ":"0|54000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"71","alt":"CTTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"71,5,0,5,0,0","location":"15000066658402","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"72","AS_RAW_MQ":"0|31609","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"72","alt":"CTTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"72,5,0,5,0,0","location":"16000073507485","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"81","AS_RAW_MQ":"0|25929","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"81","alt":"CAAAAAAAAAAAAAAAAAAA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"81,4,0,4,0,0","location":"17000007233087","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"81","AS_RAW_MQ":"0|32400","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"81","alt":"CTTTTTTTTTTTTTTTTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"81,4,0,4,0,0","location":"16000072688811","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"85","AS_RAW_MQ":"0|36000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"85","alt":"GTTTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"85,7,0,7,0,0","location":"17000030138909","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"89","AS_RAW_MQ":"0|25200","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"89","alt":"CAA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"89,5,0,6,0,0","location":"16000074086867","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"128","AS_RAW_MQ":"0|14400","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"128","alt":"CTTTTTTTTTTTTTTTTTTTTTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"128,7,0,8,0,0","location":"16000077805865","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"133","AS_RAW_MQ":"0|18000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"133","alt":"CTTTTTTTTTTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"133,8,0,8,0,0","location":"16000005651040","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"134","AS_RAW_MQ":"0|14400","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"134","alt":"CTTTTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"134,9,0,9,0,0","location":"16000082637428","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"18","AS_RAW_MQ":"0|21600","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"18","alt":"GTTTTTTTTT","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"18,2,0,11,2,8","location":"17000030878718","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"10","AS_RAW_MQ":"0|68400","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"10","alt":"T","call_AD":"0,0","call_GQ":"5","call_GT":"1/1","call_PGT":"0|1","call_PID":"10782225_A_T","call_PL":"10,5,0,10,5,10","location":"18000010782336","ref":"A","sample_id":"7"},
                    {"AS_QUALapprox":"11","AS_RAW_MQ":"0|40221","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"11","alt":"GTA","call_AD":"0,0","call_GQ":"5","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"11,5,0,11,5,11","location":"17000013467357","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"11","AS_RAW_MQ":"0|61200","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"11","alt":"A","call_AD":"0,0","call_GQ":"5","call_GT":"1/1","call_PGT":"0|1","call_PID":"5308090_ATATATATACCTATACATATATG_A","call_PL":"11,5,0,11,5,11","location":"18000005308090","ref":"ATATATATACCTATACATATATG","sample_id":"7"},
                    {"AS_QUALapprox":"13","AS_RAW_MQ":"0|60570","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"13","alt":"A","call_AD":"0,0","call_GQ":"6","call_GT":"1/1","call_PGT":"0|1","call_PID":"92047557_TAA_T","call_PL":"13,6,0,13,6,13","location":"15000092047583","ref":"AT","sample_id":"7"},
                    {"AS_QUALapprox":"13","AS_RAW_MQ":"0|60570","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"13","alt":"T","call_AD":"0,0","call_GQ":"6","call_GT":"1/1","call_PGT":"0|1","call_PID":"92047557_TAA_T","call_PL":"13,6,0,13,6,13","location":"15000092047588","ref":"TATATA","sample_id":"7"},
                    {"AS_QUALapprox":"16","AS_RAW_MQ":"0|54000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"16","alt":"CTTTTTTTTTTTTTTT","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":"0|1","call_PID":"67984630_C_CTTTTTTTTTTTTTTT","call_PL":"16,1,0,17,3,19","location":"17000067984630","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"17","AS_RAW_MQ":"0|43929","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"17","alt":"C","call_AD":"0,0","call_GQ":"7","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"17,7,0,17,7,17","location":"16000013990612","ref":"CAT","sample_id":"7"},
                    {"AS_QUALapprox":"170","AS_RAW_MQ":"0|55170","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"170","alt":"CCCTTCCCTAT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"170,9,0,10,0,0","location":"15000079059257","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"18","AS_RAW_MQ":"0|18000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"18","alt":"CAA","call_AD":"0,0","call_GQ":"3","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"18,3,0,18,3,18","location":"16000000713318","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"27","AS_RAW_MQ":"0|18000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"27","alt":"ATTTTTTTTTTTTTTT","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"27,2,0,28,3,29","location":"16000062098154","ref":"A","sample_id":"7"},
                    {"AS_QUALapprox":"27","AS_RAW_MQ":"0|22337","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"27","alt":"CAAA","call_AD":"0,0","call_GQ":"3","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"27,3,0,27,3,27","location":"16000032126382","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"28","AS_RAW_MQ":"0|28800","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"28","alt":"CAA","call_AD":"0,0","call_GQ":"3","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"28,3,0,21,3,18","location":"17000006238053","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"36","AS_RAW_MQ":"0|54000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"36","alt":"CAAAAAA","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"36,2,0,37,3,38","location":"17000005494594","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"39","AS_RAW_MQ":"0|40129","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"39","alt":"CAAAAAAAA","call_AD":"0,0","call_GQ":"1","call_GT":"0/1","call_PGT":null,"call_PID":null,"call_PL":"39,0,3,16,1,14","location":"17000020360227","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"41","AS_RAW_MQ":"0|21600","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"41","alt":"CAAAAAAAAAAAAAAA","call_AD":"0,0","call_GQ":"1","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"41,1,0,17,3,16","location":"16000071736098","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"43","AS_RAW_MQ":"0|18000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"43","alt":"GTTTTTTTTTTTTTTTT","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"43,2,0,44,3,45","location":"17000019526269","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"45","AS_RAW_MQ":"0|3600","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"45","alt":"A","call_AD":"0,0","call_GQ":"3","call_GT":"1/1","call_PGT":"0|1","call_PID":"89603914_GA_G","call_PL":"45,3,0,45,3,45","location":"16000089603975","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"77","AS_RAW_MQ":"0|18000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"77","alt":"CAAAAAAAAAAAAAAAAAAAAAA","call_AD":"0,0","call_GQ":"8","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"77,8,0,78,9,79","location":"17000012602773","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"8","AS_RAW_MQ":"0|63171","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"8","alt":"CGCTGAGCAGTGTGTGGGGTGGAGCGGTTGCTGACCGGGGCGCTGAGCAGTGTGTGGGGTGGAGCGGGTGTTGACTGGGGT","call_AD":"0,0","call_GQ":"3","call_GT":"1/1","call_PGT":"0|1","call_PID":"89898371_C_CGCTGAGCAGTGTGTGGGGTGGAGCGGTTGCTGACCGGGGCGCTGAGCAGTGTGTGGGGTGGAGCGGGTGTTGACTGGGGT","call_PL":"8,3,0,14,10,22","location":"16000089898371","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"81","AS_RAW_MQ":"0|45555","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"81","alt":"AATATAAAATAATTACATATTATATATTTATATATAAT","call_AD":"0,0","call_GQ":"4","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"81,4,0,84,7,87","location":"16000065128304","ref":"A","sample_id":"7"},
                    {"AS_QUALapprox":"88","AS_RAW_MQ":"0|16401","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"88","alt":"GAAAAAAAAAAAAAA","call_AD":"0,0","call_GQ":"5","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"88,5,0,89,6,90","location":"17000043305380","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"88","AS_RAW_MQ":"0|36784","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"88","alt":"GTTTTTTTTTTTTT","call_AD":"0,0","call_GQ":"5","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"88,5,0,89,6,90","location":"16000083080871","ref":"G","sample_id":"7"},
                    {"AS_QUALapprox":"98","AS_RAW_MQ":"0|56809","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"98","alt":"CAAAAAAAA","call_AD":"0,0","call_GQ":"2","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"98,8,0,32,2,24","location":"15000077053451","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"161","AS_RAW_MQ":"0|46800","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"161","alt":"CTTTTTTTTTTTTTTTTTT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"161,11,0,11,0,0","location":"17000056368276","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"179","AS_RAW_MQ":"0|18000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"179","alt":"TC","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"179,15,0,15,0,0","location":"16000023902800","ref":"T","sample_id":"7"},
                    {"AS_QUALapprox":"199","AS_RAW_MQ":"0|43200","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"199","alt":"CAAAAAAAAAAAAAAAAA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"199,13,0,14,0,1","location":"17000028353520","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"223","AS_RAW_MQ":"0|39600","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"223","alt":"CATATATAT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"223,14,0,14,0,0","location":"17000074365044","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"24","AS_RAW_MQ":"0|18000","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"24","alt":"CAAAAAAAAAAAAA","call_AD":"0,0","call_GQ":"1","call_GT":"0/1","call_PGT":null,"call_PID":null,"call_PL":"24,0,40,24,1,25","location":"18000005328386","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"250","AS_RAW_MQ":"0|34301","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"250","alt":"CCAAA","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"250,27,0,27,0,0","location":"16000012606864","ref":"C","sample_id":"7"},
                    {"AS_QUALapprox":"266","AS_RAW_MQ":"0|43200","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"266","alt":"G","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"266,29,0,30,0,1","location":"15000061980212","ref":"A","sample_id":"7"},
                    {"AS_QUALapprox":"276","AS_RAW_MQ":"0|33625","AS_RAW_MQRankSum":null,"AS_RAW_ReadPosRankSum":null,"AS_SB_TABLE":"0,0|0,0","AS_VarDP":"0|0","QUALapprox":"276","alt":"GT","call_AD":"0,0","call_GQ":"0","call_GT":"1/1","call_PGT":null,"call_PID":null,"call_PL":"276,24,0,25,0,1","location":"15000062888285","ref":"G","sample_id":"7"}
                    ]
                    """);
            for (int i = 0; i < testData.length(); i++)
            {
                JSONObject jsonObj = testData.getJSONObject(i);
                testWriter.write(jsonObj);
            }
            testWriter.close();
        } catch (IOException exception) {
            exception.printStackTrace();
            System.exit(1);
        }



        logger.info("Exiting early, trying to set up Parquet writing");
        System.exit(0);
         */

        // If BQ, check the load status table to see if this sample has already been loaded.
        if (outputType == CommonCode.OutputType.BQ) {
            loadStatus = new LoadStatus(projectID, datasetName, loadStatusTableName);

            LoadStatus.LoadState state = loadStatus.getSampleLoadState(sampleId);

            // Legacy "FINISHED" state indicates variants and references completely loaded.
            if (state.isFinished()) {
                logger.info("Sample id " + sampleId + " was detected as already loaded, exiting successfully.");
                System.exit(0);
            }

            if (enableReferenceRanges) {
                if (state.areReferencesLoaded()) {
                    logger.info("Sample ID {}: Reference ranges writing enabled but REFERENCES_LOADED status row found, skipping.", sampleId);
                } else {
                    refRangesRowsExist = RefCreator.doRowsExistFor(outputType, projectID, datasetName, tableNumber, sampleId);
                    if (refRangesRowsExist) {
                        logger.warn("Reference ranges enabled for sample id = {}, name = {} but preexisting ref_ranges rows found, skipping ref_ranges writes.",
                                sampleId, sampleName);
                    }
                    shouldWriteReferencesLoadedStatusRow = true;
                }
            }

            if (enableVet) {
                if (state.areVariantsLoaded()) {
                    logger.info("Sample ID {}: Variant writing enabled but VARIANTS_LOADED status row found, skipping.", sampleId);
                } else {
                    vetRowsExist = VetCreator.doRowsExistFor(outputType, projectID, datasetName, tableNumber, sampleId);
                    if (vetRowsExist) {
                        logger.warn("Vet enabled for sample id = {}, name = {} but preexisting vet rows found, skipping vet writes.",
                                sampleId, sampleName);
                    }
                    shouldWriteVariantsLoadedStatusRow = true;
                }
            }

            if (enableVCFHeaders) {
                if (state.areHeadersLoaded()) {
                    logger.info("Sample ID {}: VCF Header writing enabled but HEADERS_LOADED status row found, skipping.", sampleId);
                } else {
                    vcfHeaderRowsExist = VcfHeaderLineScratchCreator.doScratchRowsExistFor(projectID, datasetName, sampleId);
                    if (vcfHeaderRowsExist) {
                        logger.warn("VCF header writing enabled for sample id = {}, name = {} but preexisting VCF header scratch rows found, skipping VCF header writes.",
                                sampleId, sampleName);
                    } else {
                        vcfHeaderRowsExist = VcfHeaderLineScratchCreator.doNonScratchRowsExistFor(projectID, datasetName, sampleId);
                        if (vcfHeaderRowsExist) {
                            logger.warn("VCF header writing enabled for sample id = {}, name = {} but preexisting VCF header non-scratch rows found, skipping VCF header writes.",
                                    sampleId, sampleName);
                        }
                    }
                    shouldWriteVCFHeadersLoadedStatusRow = true;
                }
            }
        }

        // This needs to be called *outside* the "if outputType == BQ" because it side-effects the initialization of
        // some class members that the CreateVariantIngestFilesTest expects to be initialized.
        SAMSequenceDictionary seqDictionary = initializeGQConfigurationAndIntervals();

        if (enableReferenceRanges && !refRangesRowsExist) {
            refCreator = new RefCreator(sampleIdentifierForOutputFileName, sampleId, tableNumber, seqDictionary, gqStatesToIgnore, outputDir, outputType, enableReferenceRanges, projectID, datasetName, storeCompressedReferences);

            // The ploidy table is really only needed for inferring reference ploidy, as variants store their genotypes
            // directly.  If we're not ingesting references, we can't compute and ingest ploidy either
            samplePloidyCreator = new SamplePloidyCreator(sampleId, projectID, datasetName, outputType);
        }

        if (enableVet && !vetRowsExist) {
            vetCreator = new VetCreator(sampleIdentifierForOutputFileName, sampleId, tableNumber, outputDir, outputType, projectID, datasetName, forceLoadingFromNonAlleleSpecific, skipLoadingVqsrFields, variantRowSchema);
        }
        if (enableVCFHeaders && !vcfHeaderRowsExist) {
            buildAllVcfLineHeaders();
            vcfHeaderLineScratchCreator = new VcfHeaderLineScratchCreator(sampleId, projectID, datasetName);
        }

        logger.info("enableReferenceRanges = {}, enableVet = {}, enableVCFHeaders = {}",
                enableReferenceRanges, enableVet, enableVCFHeaders);
        logger.info("shouldWriteReferencesLoadedStatus = {}, shouldWriteVariantsLoadedStatus = {}, shouldWriteVCFHeadersLoadedStatus = {}",
                shouldWriteReferencesLoadedStatusRow, shouldWriteVariantsLoadedStatusRow, shouldWriteVCFHeadersLoadedStatusRow);

        if (refCreator == null && vetCreator == null && vcfHeaderLineScratchCreator == null &&
                !shouldWriteReferencesLoadedStatusRow && !shouldWriteVariantsLoadedStatusRow && !shouldWriteVCFHeadersLoadedStatusRow) {
            logger.info("No data to be written, exiting successfully.");
            System.exit(0);
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (!enableReferenceRanges && !enableVet && enableVCFHeaders) {
            // Nothing to do here, exit ASAP.
            return;
        }

        // get the intervals this variant covers
        final GenomeLoc variantGenomeLoc = intervalArgumentGenomeLocSortedSet.getGenomeLocParser().createGenomeLoc(variant.getContig(), variant.getStart(), variant.getEnd());
        final List<GenomeLoc> intervalsToWrite = intervalArgumentGenomeLocSortedSet.getOverlapping(variantGenomeLoc);

        if (intervalsToWrite.isEmpty()) {
            throw new IllegalStateException("There are no intervals being covered by this variant, something went wrong with interval parsing");
        }

        // take the first interval(assuming this is returned in order) and make sure if it's a variant, that it starts at/after the interval start
        // we are going to ignore any deletions that start before an interval.
        if (!variant.isReferenceBlock() && intervalsToWrite.get(0).getStart() > variant.getStart()){
            return;
        }

        // if the only alt allele for a variant is `*`, we ignore it
        if (!variant.isReferenceBlock() &&  variant.getAlternateAlleles().size() == 2 && variant.hasAlternateAllele(Allele.SPAN_DEL)){
            return;
        }

        try {
            // write to VET if NOT reference block and NOT a no call
            if (!variant.isReferenceBlock() && !isNoCall(variant)) {
                if (vetCreator != null) {
                    vetCreator.apply(variant);
                }
            }
        } catch (IOException ioe) {
            throw new GATKException("Error writing VET", ioe);
        }

        try {
            if (refCreator != null) {
                refCreator.apply(variant, intervalsToWrite);
            }
        } catch (IOException ioe) {
            throw new GATKException("Error writing reference ranges", ioe);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (vcfHeaderLineScratchCreator != null) {
            try {
                vcfHeaderLineScratchCreator.apply(allLineHeaders);
            } catch (IOException ioe) {
                throw new GATKException("Error writing scratch header data", ioe);
            }
            // Wait until all data has been submitted and in pending state to commit
            vcfHeaderLineScratchCreator.commitData();
        }

        if (refCreator != null) {
            if ((gqStatesToIgnore.size() != 1) || (!gqStatesToIgnore.contains(GQStateEnum.ZERO))) {
                // We will write missing intervals as ZERO ('GQ0') unless that is the (ONLY???) GQ state that we are dropping.
                // If ZERO/GQ0 is the ONLY state that we are dropping then we do not write those intervals.
                try {
                    refCreator.writeMissingIntervals(intervalArgumentGenomeLocSortedSet);
                } catch (IOException ioe) {
                    throw new GATKException("Error writing missing intervals", ioe);
                }
            }
            // Wait until all data has been submitted and in pending state to commit
            refCreator.commitData();

            // this is likely an unnecessary check as it currently stands, but it can't hurt to have it in case we
            // later separate their creation, throwing the ploidy stuff explicity behind a flag
            if (samplePloidyCreator != null) {
                try {
                    samplePloidyCreator.apply(refCreator.getReferencePloidyData(), refCreator.getTotalRefEntries());
                } catch (IOException ioe) {
                    throw new GATKException("Error writing ploidy data", ioe);
                }

                samplePloidyCreator.commitData();
            }
        }

        if (vetCreator != null) {
            vetCreator.commitData();
        }

        if (outputType == CommonCode.OutputType.BQ) {
            if (shouldWriteVCFHeadersLoadedStatusRow) loadStatus.writeHeadersLoadedStatus(sampleId);
            if (shouldWriteVariantsLoadedStatusRow) loadStatus.writeVariantsLoadedStatus(sampleId);
            if (shouldWriteReferencesLoadedStatusRow) loadStatus.writeReferencesLoadedStatus(sampleId);
        }

        return 0;
    }

    @Override
    public void closeTool() {
        if (refCreator != null) {
            refCreator.closeTool();
        }
        if (vetCreator != null) {
            vetCreator.closeTool();
        }
        if (samplePloidyCreator != null) {
            samplePloidyCreator.closeTool();
        }
        if (vcfHeaderLineScratchCreator != null) {
            vcfHeaderLineScratchCreator.closeTool();
        }
    }

    private void buildAllVcfLineHeaders() {
        // get a set of header "lines" (command line INFO lines, chunks of non-command-line INFO lines)
        // to put into the header line scratch table
        Set<String> nonCommandLineHeaders = new HashSet<>();
        for (VCFHeaderLine line : getHeaderForVariants().getMetaDataInInputOrder()) {
            if (line.getKey().contains("CommandLine")) {
                allLineHeaders.put(line.toString(), true);
            }
            else {
                nonCommandLineHeaders.add(line.toString());
            }
        }
        allLineHeaders.put(StringUtils.join(nonCommandLineHeaders), false);
    }

    // TODO: there's almost certainly a better name for this chunk of logic...
    private SAMSequenceDictionary initializeGQConfigurationAndIntervals() {
        // To set up the missing positions
        SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();

        final GenomeLocParser genomeLocParser = new GenomeLocParser(seqDictionary);
        intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocParser, IntervalUtils.genomeLocsFromLocatables(genomeLocParser, intervalArgumentCollection.getIntervals(seqDictionary)));

        if (gqStateToIgnore != null) {
            gqStatesToIgnore.add(gqStateToIgnore);
            if (dropAboveGqThreshold) {
                gqStatesToIgnore.addAll(RefCreator.getGQStateEnumGreaterThan(gqStateToIgnore));
            }
        }
        return seqDictionary;
    }
}
