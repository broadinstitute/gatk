package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.cloud.bigquery.storage.v1beta2.BatchCommitWriteStreamsRequest;
import com.google.cloud.bigquery.storage.v1beta2.BatchCommitWriteStreamsResponse;
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
import org.broadinstitute.hellbender.utils.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

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

    private PetCreator petTsvCreator;
    private VetCreator vetTsvCreator;
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

        petTsvCreator = new PetCreator(sampleIdentifierForOutputFileName, sampleId, tableNumber, seqDictionary, gqStateToIgnore, dropAboveGqThreshold, outputDir, outputType, enablePet, enableReferenceRanges, projectID, datasetName);

        if (enableVet) {
            vetTsvCreator = new VetCreator(sampleIdentifierForOutputFileName, sampleId, tableNumber, outputDir, outputType, projectID, datasetName);
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
                if (enableVet) vetTsvCreator.apply(variant, readsContext, referenceContext, featureContext);
            }
        } catch (IOException ioe) {
            throw new GATKException("Error writing VET", ioe);
        }

        try {
            petTsvCreator.apply(variant, intervalsToWrite);
        } catch (IOException ioe) {
            throw new GATKException("Error writing PET", ioe);
        }

    }


    @Override
    public Object onTraversalSuccess() {
        try {
            petTsvCreator.writeMissingIntervals(intervalArgumentGenomeLocSortedSet);
        } catch (IOException ioe) {
            throw new GATKException("Error writing missing intervals", ioe);
        }
        // Wait until all data has been submitted and in pending state to commit
        vetTsvCreator.completeCreation();
        petTsvCreator.completeCreation();
        return 0;
    }

    @Override
    public void closeTool() {
        if (petTsvCreator != null) {
            petTsvCreator.closeTool();
        }
        if (vetTsvCreator != null) {
            vetTsvCreator.closeTool();;
        }
    }

}
