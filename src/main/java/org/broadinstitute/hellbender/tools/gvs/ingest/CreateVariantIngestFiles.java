package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.api.core.ApiFuture;
import com.google.cloud.bigquery.BigQuery;
import com.google.cloud.bigquery.BigQueryOptions;
import com.google.cloud.bigquery.Field;
import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.Schema;
import com.google.cloud.bigquery.StandardTableDefinition;
import com.google.cloud.bigquery.Table;
import com.google.cloud.bigquery.TableId;
import com.google.cloud.bigquery.TableResult;
import com.google.cloud.bigquery.storage.v1beta2.AppendRowsResponse;
import com.google.cloud.bigquery.storage.v1beta2.JsonStreamWriter;
import com.google.cloud.bigquery.storage.v1beta2.TableFieldSchema;
import com.google.cloud.bigquery.storage.v1beta2.TableName;
import com.google.cloud.bigquery.storage.v1beta2.TableSchema;
import com.google.protobuf.Descriptors;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.gvs.common.ChromosomeEnum;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.GQStateEnum;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.IngestUtils;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;

/**
 * Ingest variant walker
 */
@CommandLineProgramProperties(
        summary = "Exome and Genome Ingest tool for the Joint Genotyping in Big Query project",
        oneLineSummary = "Ingest tool for BQJG",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = true
)
public final class CreateVariantIngestFiles extends VariantWalker {
    static final Logger logger = LogManager.getLogger(CreateVariantIngestFiles.class);

    private RefCreator refCreator;
    private VetCreator vetCreator;
    private enum LoadStatus { STARTED, FINISHED };
    private TableName loadStatusTable;
    private TableSchema loadStatusTableSchema;

    private GenomeLocSortedSet intervalArgumentGenomeLocSortedSet;

    private String sampleName;
    private String sampleId;
    private List<SimpleInterval> userIntervals;

    // Inside the parent directory, a directory for each chromosome will be created, with a pet directory and vet directory in each one.
    // Each pet and vet directory will hold all of the pet and vet tsvs for each sample
    // A sample_info directory will be created, with a sample_info tsv for each sample

//    @Argument(fullName = "output-path",
//            shortName = "VPO",
//            doc = "Path to the directory where the variants TSVs and positions expanded TSVs should be written")
//    public GATKPathSpecifier parentOutputDirectory = null;
//    public Path parentDirectory = null;

    @Argument(fullName = "ref-block-gq-to-ignore",
            shortName = "IG",
            doc = "Ref Block GQ band to ignore, bands of 10 e.g 0-9 get combined to 0, 20-29 get combined to 20",
            optional = true)
    public GQStateEnum gqStateToIgnore = GQStateEnum.SIXTY;

    @Argument(fullName = "ignore-above-gq-threshold",
    shortName = "GTIG",
    doc = "in addition to dropping the gq block specified by ref-block-gq-to-ignore, also drop higher gq blocks",
    optional = true)
    public boolean dropAboveGqThreshold = false;

    @Argument(fullName = "enable-reference-ranges",
            shortName = "rr",
            doc = "write reference ranges data",
            optional = true)
    public boolean enableReferenceRanges = false;

    @Argument(fullName = "enable-vet",
            shortName = "ev",
            doc = "write vet data",
            optional = true)
    public boolean enableVet = true;

    @Argument(fullName = "enable-pet",
            shortName = "ep",
            doc = "write pet data",
            optional = true)
    public boolean enablePet = true;

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
    private String refVersion = "37";

    @Argument(
            fullName = "output-directory",
            doc = "directory for output tsv files",
            optional = true)
    private File outputDir = new File(".");

    @Argument(
            fullName = "project-id",
            doc = "ID of the Google Cloud project where the dataset for pet and vet tables exist",
            optional = true
    )
    protected String projectID = null;

    @Argument(
            fullName = "dataset-name",
            doc = "Name of the dataset to update pet and vet tables",
            optional = true
    )
    protected String datasetName = null;


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

    private void writeLoadStatus(LoadStatus status) {

        // This uses the _default stream since it (a) commits immediately and (b) doesn't count
        // towards the CreateStreamWriter quota
        try (JsonStreamWriter writer =
                     JsonStreamWriter.newBuilder(loadStatusTable.toString(), loadStatusTableSchema).build()) {

            // Create a JSON object that is compatible with the table schema.
            JSONArray jsonArr = new JSONArray();
            JSONObject jsonObject = new JSONObject();
            jsonObject.put("sample_id", Long.parseLong(sampleId));
            jsonObject.put("status", status.toString());
            jsonObject.put("event_timestamp", System.currentTimeMillis());
            jsonArr.put(jsonObject);

            ApiFuture<AppendRowsResponse> future = writer.append(jsonArr);
            AppendRowsResponse response = future.get();

            logger.info("Load status " + status + " appended successfully");
        } catch (ExecutionException | InterruptedException | Descriptors.DescriptorValidationException | IOException e) {
            // If the wrapped exception is a StatusRuntimeException, check the state of the operation.
            // If the state is INTERNAL, CANCELLED, or ABORTED, you can retry. For more information, see:
            // https://grpc.github.io/grpc-java/javadoc/io/grpc/StatusRuntimeException.html
            throw new GATKException(e.getMessage());
        }
    }

    @Override
    public void onTraversalStart() {
        //set up output directory
        if (!outputDir.exists() && !outputDir.mkdir()) {
            throw new RuntimeIOException("Unable to create directory: " + outputDir.getAbsolutePath());
        }

        // Set reference version -- TODO remove this in the future, also, can we get ref version from the header?
        ChromosomeEnum.setRefVersion(refVersion);

        // TODO should we reuse the SampleList class or move these methods there?
        // TODO if you change here, also change in CreateArrayIngestFiles
        // Get sample name
        final VCFHeader inputVCFHeader = getHeaderForVariants();
        sampleName = sampleNameParam == null ? IngestUtils.getSampleName(inputVCFHeader) : sampleNameParam;
        if (sampleIdParam == null && sampleMap == null) {
            throw new IllegalArgumentException("One of sample-id or sample-name-mapping must be specified");
        }
        if (sampleIdParam != null) {
            sampleId = String.valueOf(sampleIdParam);
        } else {
            sampleId = IngestUtils.getSampleId(sampleName, sampleMap);
        }

        // TODO when we pass in the full file path or gvs_id as an input arg, use path here instead or gvs_id in addition
        // use input gvcf file name instead of sample name in output filenames
        String sampleIdentifierForOutputFileName = getInputFileName();

        // Mod the sample directories
        int sampleTableNumber = IngestUtils.getTableNumber(sampleId, IngestConstants.partitionPerTable);
        String tableNumber = String.format("%03d", sampleTableNumber);

        // To set up the missing positions
        SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();
        userIntervals = intervalArgumentCollection.getIntervals(seqDictionary);

        final GenomeLocSortedSet genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));
        intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocSortedSet.getGenomeLocParser(), IntervalUtils.genomeLocsFromLocatables(genomeLocSortedSet.getGenomeLocParser(), intervalArgumentCollection.getIntervals(seqDictionary)));

        if (enablePet || enableReferenceRanges) {
            refCreator = new RefCreator(sampleIdentifierForOutputFileName, sampleId, tableNumber, seqDictionary, gqStateToIgnore, dropAboveGqThreshold, outputDir, outputType, enablePet, enableReferenceRanges, projectID, datasetName);
        }

        if (enableVet) {
            vetCreator = new VetCreator(sampleIdentifierForOutputFileName, sampleId, tableNumber, outputDir, outputType, projectID, datasetName);
        }

        // check the load status table to see if this sample has already been loaded...
        if (outputType == CommonCode.OutputType.BQ) {
            loadStatusTable = TableName.of(projectID, datasetName, loadStatusTableName);

            BigQuery bigquery = BigQueryUtils.getBigQueryEndPoint(projectID);
            Table table = bigquery.getTable(TableId.of(projectID, datasetName, loadStatusTableName));
            loadStatusTableSchema = getLoadStatusTableSchema();

            StandardTableDefinition tdd = table.getDefinition();
            if (BigQueryUtils.getEstimatedRowsInStreamingBuffer(projectID, datasetName, loadStatusTableName) > 0 ) {
                logger.info("Found estimated rows in streaming buffer!!! " + tdd.getStreamingBuffer().getEstimatedRows());
            }

            verifySampleIsNotLoaded();
        }

    }

    private TableSchema getLoadStatusTableSchema() {
        TableSchema.Builder builder = TableSchema.newBuilder();
        builder.addFields(
                TableFieldSchema.newBuilder().setName(SchemaUtils.SAMPLE_ID_FIELD_NAME).setType(TableFieldSchema.Type.NUMERIC).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName(SchemaUtils.LOAD_STATUS_FIELD_NAME).setType(TableFieldSchema.Type.STRING).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        builder.addFields(
                TableFieldSchema.newBuilder().setName(SchemaUtils.LOAD_STATUS_EVENT_TIMESTAMP_NAME).setType(TableFieldSchema.Type.TIMESTAMP).setMode(TableFieldSchema.Mode.REQUIRED).build()
        );
        return builder.build();
    }

    private void verifySampleIsNotLoaded() {
        String query = "SELECT " + SchemaUtils.LOAD_STATUS_FIELD_NAME +
                       " FROM `" + projectID + "." + datasetName + "." + loadStatusTableName + "` " +
                       " WHERE " + SchemaUtils.SAMPLE_ID_FIELD_NAME + " = " + sampleId;

        TableResult result = BigQueryUtils.executeQuery(projectID, query, true, null);

        int startedCount = 0;
        int finishedCount = 0;
        for ( final FieldValueList row : result.iterateAll() ) {
            final String status = row.get(0).getStringValue();
            if (LoadStatus.STARTED.toString().equals(status)) {
                startedCount++;
            } else if (LoadStatus.FINISHED.toString().equals(status)) {
                finishedCount++;
            }

        }

        // if fully loaded, exit successfully!
        if (startedCount == 1 && finishedCount == 1) {
            logger.info("Sample id " + sampleId + " was detected as already loaded, exiting successfully.");
            System.exit(0);
        }

        // otherwise if there are any records, exit
        if (startedCount > 0 || finishedCount > 0) {
            throw new GATKException("Sample Id " + sampleId + " has already been partially loaded!");
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // get the intervals this variant covers
        final GenomeLoc variantGenomeLoc = intervalArgumentGenomeLocSortedSet.getGenomeLocParser().createGenomeLoc(variant.getContig(), variant.getStart(), variant.getEnd());
        final List<GenomeLoc> intervalsToWrite = intervalArgumentGenomeLocSortedSet.getOverlapping(variantGenomeLoc);

        if (intervalsToWrite.size() == 0){
            throw new IllegalStateException("There are no intervals being covered by this variant, something went wrong with interval parsing");
        }

        // take the first interval(assuming this is returned in order) and make sure if its a variant, that it starts at/after the interval start
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
                if (enableVet) vetCreator.apply(variant, readsContext, referenceContext, featureContext);
            }
        } catch (IOException ioe) {
            throw new GATKException("Error writing VET", ioe);
        }

        try {
            if (refCreator != null) {
                refCreator.apply(variant, intervalsToWrite);
            }
        } catch (IOException ioe) {
            throw new GATKException("Error writing PET", ioe);
        }

    }


    @Override
    public Object onTraversalSuccess() {
        if (outputType == CommonCode.OutputType.BQ) {
            writeLoadStatus(LoadStatus.STARTED);
        }

        if (refCreator != null) {
            try {
                refCreator.writeMissingIntervals(intervalArgumentGenomeLocSortedSet);
            } catch (IOException ioe) {
                throw new GATKException("Error writing missing intervals", ioe);
            }
            // Wait until all data has been submitted and in pending state to commit
            refCreator.commitData();
        }

        if (vetCreator != null && enableVet) {
            vetCreator.commitData();
        }

        // upload the load status table
        if (outputType == CommonCode.OutputType.BQ) {
            writeLoadStatus(LoadStatus.FINISHED);
        }

        return 0;
    }

    @Override
    public void closeTool() {
        if (refCreator != null) {
            refCreator.closeTool();
        }
        if (vetCreator != null) {
            vetCreator.closeTool();;
        }
    }

}
