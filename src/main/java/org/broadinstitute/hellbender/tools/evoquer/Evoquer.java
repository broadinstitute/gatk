package org.broadinstitute.hellbender.tools.evoquer;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_MappingQualityRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_QualByDepth;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StandardAnnotation;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * Evoquer ("Evoker"):
 * Extract Variants Out of big QUERy.
 *
 * Connects to a BigQuery database and creates a VCF file from the given interval list based on the data in the database.
 *
 * The BigQuery database is assumed to have a schema and data that supports creation of {@link htsjdk.variant.variantcontext.VariantContext}
 * objects, and all data in this database are assumed to be based on the HG38 reference.
 *
 * Evoquer requires a sequence dictionary that can be inferred from a reference FASTA file or passed in directly.
 *
 * Example invocations are:
 *
 *   ./gatk Evoquer --intervals chr2:10000-1741634 -R ~/references/Homo_sapiens_assembly38.fasta
 *   ./gatk Evoquer --intervals chr2:10000-1741634 -R ~/references/Homo_sapiens_assembly38.fasta -O evoquer.g.vcf
 *   ./gatk Evoquer --intervals chr2:10000-1741634 --sequence-dictionary ~/references/Homo_sapiens_assembly38.dict -O evoquer.vcf
 *
 * Created by jonn on 4/16/19.
 */
@CommandLineProgramProperties(
        summary = "(\"Evoker\") - Extract Variants Out of big QUERy.  Tool to extract VCF information from a BigQuery table containing genomic variant data (joing genotyping / GVCF-equivalent data).",
        oneLineSummary = "Tool to extract VCF information from a BigQuery table containing genomic variant data (joing genotyping / GVCF-equivalent data).",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class Evoquer extends GATKTool {
    private static final Logger logger = LogManager.getLogger(Evoquer.class);

    //==================================================================================================================
    // Public Static Members:

    public static final String DEFAULT_PROJECT_ID = "broad-dsp-spec-ops";

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    @Argument(
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            doc = "Output VCF file to which annotated variants should be written.",
            optional = false
    )
    private String outputVcfPathString = null;

    @Argument(
            fullName = "project-id",
            doc = "ID of the Google Cloud project containing the dataset and tables from which to pull variant data",
            optional = false
    )
    private String projectID = null;
    
    @Argument(
            fullName = "dataset-map",
            doc = "Path to a file containing a mapping from contig name to BigQuery dataset name. Each line should consist of a contig name, followed by a tab, followed by a dataset name.",
            optional = false
    )
    private String datasetMapFile = null;

    @Argument(
            fullName = "query-record-limit",
            doc = "Limits the maximum number of records returned from each query on BiqQuery (for profiling/debugging purposes only). Set to 0 for no limit.",
            optional = true)
    private int queryRecordLimit = 0;

    @Argument(
            fullName = "run-query-only",
            doc = "If true, just do the query against BigQuery and retrieve the resulting records, but don't write a VCF",
            optional = true)
    private boolean runQueryOnly = false;

    @Argument(
            fullName = "disable-gnarly-genotyper",
            doc = "If true, don't run the Gnarly Genotyper after combining variant records for each site",
            optional = true)
    private boolean disableGnarlyGenotyper = false;

    @Argument(
            fullName = "skip-malformed-as-annotations",
            doc = "If true, treat sample records with malformed allele-specific annotations as missing",
            optional = true)
    private boolean skipMalformedASAnnotations = false;

    @Argument(
            fullName = "print-debug-information",
            doc = "If true, print extra debugging output",
            optional = true)
    private boolean printDebugInformation = false;

    private VariantContextWriter vcfWriter = null;
    private EvoquerEngine evoquerEngine;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    /**
     * {@inheritDoc}
     *
     * {@link Evoquer} requires intervals.
     */
    @Override
    public boolean requiresIntervals() {
        return true;
    }

    /**
     * {@inheritDoc}
     *
     * {@link Evoquer} doesn't require a reference, but DOES require a sequence dictionary for interval processing.
     */
    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean useVariantAnnotations() { return true; }

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Arrays.asList(
                StandardAnnotation.class, AS_StandardAnnotation.class
        );
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        final Map<String, EvoquerDataset> datasetMap = loadDatasetMapFile(datasetMapFile);

        logger.info("Disabling index creation on the output VCF, since this tool writes an unsorted VCF");
        createOutputVariantIndex = false;
        
        vcfWriter = createVCFWriter(IOUtils.getPath(outputVcfPathString));

        final VariantAnnotatorEngine annotationEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), null, Collections.emptyList(), false, false);

        // Set up our EvoquerEngine:
        evoquerEngine = new EvoquerEngine(projectID,
                                          datasetMap,
                                          queryRecordLimit,
                                          vcfWriter,
                                          getDefaultToolVCFHeaderLines(),
                                          annotationEngine,
                                          reference,
                                          runQueryOnly,
                                          disableGnarlyGenotyper,
                                          skipMalformedASAnnotations,
                                          printDebugInformation,
                                          progressMeter);

        vcfWriter.writeHeader(evoquerEngine.getHeader());
    }

    @Override
    public void traverse() {
        progressMeter.setRecordsBetweenTimeChecks(100L);

        for (final SimpleInterval interval : getTraversalIntervals()) {
            evoquerEngine.evokeInterval(interval);
        }
    }

    @Override
    protected void onShutdown() {
        super.onShutdown();

        logger.info(String.format("***Processed %s total sites", evoquerEngine.getTotalNumberOfSites()));
        logger.info(String.format("***Skipped %d/%d total variant records due to malformed allele-specific annotations", evoquerEngine.getNumberOfSkippedVariants(), evoquerEngine.getTotalNumberOfVariants()));

        // Close up our writer if we have to:
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

    private Map<String, EvoquerDataset> loadDatasetMapFile(final String datasetMapFile) {
        final Path datasetMapPath = IOUtils.getPath(datasetMapFile);
        final Map<String, EvoquerDataset> datasetMap = new HashMap<>();
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        try {
            final List<String> lines = Files.readAllLines(datasetMapPath);
            if ( lines.isEmpty() ) {
                throw new UserException.BadInput("Dataset map file " + datasetMapFile + " is empty");
            }

            for ( final String line : lines ) {
                final String[] split = line.split("\\s+");
                if ( split.length != 4 || split[0].trim().isEmpty() || split[1].trim().isEmpty() || split[2].trim().isEmpty() || split[3].trim().isEmpty() ) {
                    throw new UserException.BadInput("Dataset map file " + datasetMapFile + " contains a malformed line: \"" + line + "\". " +
                            "Format should be: \"contig  dataset  petTableName  vetTableName\"");
                }

                final String contig = split[0];
                final String contigDataset = split[1];
                final String contigPetTableName = split[2];
                final String contigVetTableName = split[3];

                if ( sequenceDictionary.getSequenceIndex(contig) == -1 ) {
                    throw new UserException.BadInput("Dataset map file " + datasetMapFile + " contains a contig not in our dictionary: " + contig);
                }

                datasetMap.put(contig, new EvoquerDataset(contigDataset, contigPetTableName, contigVetTableName));
            }
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(datasetMapFile, e);
        }

        return datasetMap;
    }

    public static final class EvoquerDataset {
        private final String datasetName;
        private final String petTableName;
        private final String vetTableName;

        public EvoquerDataset( final String datasetName, final String petTableName, final String vetTableName ) {
            this.datasetName = datasetName;
            this.petTableName = petTableName;
            this.vetTableName = vetTableName;
        }

        public String getDatasetName() { return datasetName; }

        public String getPetTableName() { return petTableName; }

        public String getVetTableName() { return vetTableName; }
    }

}
