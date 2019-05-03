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
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
            optional = true
    )
    private String outputVcfPathString;

    @Argument(
            fullName = "project-id",
            doc = "ID of the Google Cloud project containing the dataset and tables from which to pull variant data",
            optional = true
    )
    private String projectID = DEFAULT_PROJECT_ID;
    
    @Argument(
            fullName = "dataset-map",
            doc = "Path to a file containing a mapping from contig name to BigQuery dataset name. Each line should consist of a contig name, followed by a tab, followed by a dataset name.",
            optional = false
    )
    private String datasetMapFile;

    @Argument(
            fullName = "query-record-limit",
            doc = "Limits the maximum number of records returned from each query on BiqQuery (for profiling/debugging purposes only). Set to 0 for no limit.",
            optional = true)
    private int queryRecordLimit = 0;

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
        return false;
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        // Validate that somehow we have a dictionary:
        if ( getBestAvailableSequenceDictionary() == null ) {
            throw new UserException("Error: You must provide either a reference file or a sequence dictionary using:\n" +
                    "`-R REFERENCE_FILE`\n" +
                    "OR\n" +
                    "--sequence-dictionary SEQUENCE_DICTIONARY");
        }

        final Map<String, String> datasetMap = loadDatasetMapFile(datasetMapFile);

        // Set up our EvoquerEngine:
        evoquerEngine = new EvoquerEngine(projectID,
                                        datasetMap,
                                        queryRecordLimit);
    }

    @Override
    protected void onShutdown() {
        super.onShutdown();

        // Close up our writer if we have to:
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

    @Override
    public void traverse() {
        // Get our variants from BigQuery:
        final List<VariantContext> variants = evoquerEngine.evokeIntervals(getTraversalIntervals());

        // Now we can write out the VariantContexts:
        if ( outputVcfPathString != null ) {
            vcfWriter = createVCFWriter(IOUtils.getPath(outputVcfPathString) );
            vcfWriter.writeHeader( evoquerEngine.generateVcfHeader(getDefaultToolVCFHeaderLines(), getBestAvailableSequenceDictionary()) );

            logger.info( "Created the following variants:" );
            for ( final VariantContext variantContext : variants ) {
                logger.info( variantContext.toStringWithoutGenotypes() );
                vcfWriter.add( variantContext );
            }
        }
    }

    private Map<String, String> loadDatasetMapFile(final String datasetMapFile) {
        final Path datasetMapPath = IOUtils.getPath(datasetMapFile);
        final Map<String, String> datasetMap = new HashMap<>();
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();

        try {
            final List<String> lines = Files.readAllLines(datasetMapPath);
            if ( lines.isEmpty() ) {
                throw new UserException.BadInput("Dataset map file " + datasetMapFile + " is empty");
            }

            for ( final String line : lines ) {
                final String[] split = line.split("\\t",-1);
                if ( split.length != 2 || split[0].trim().isEmpty() || split[1].trim().isEmpty() ) {
                    throw new UserException.BadInput("Dataset map file " + datasetMapFile + " contains a malformed line: " + line);
                }

                final String contig = split[0];
                final String contigDataset = split[1];
                if ( sequenceDictionary.getSequenceIndex(contig) == -1 ) {
                    throw new UserException.BadInput("Dataset map file " + datasetMapFile + " contains a contig not in our dictionary: " + contig);
                }

                datasetMap.put(contig, contigDataset);
            }
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(datasetMapFile, e);
        }

        return datasetMap;
    }

}
