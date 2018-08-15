package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.MetagenomicsProgramGroup;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import scala.Tuple2;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * Build an annotated taxonomy datafile for a given microbe reference. The output file from this tool is required to run the PathSeq pipeline.
 *
 * <p>The tool reads the list of sequence accessions from the given reference. For each accession, it looks up the NCBI
 * taxonomic ID of the corresponding organism and builds a taxonomic tree containing only organisms that are
 * represented in the reference. The reference should only contain sequences from NCBI RefSeq and/or Genbank databases.</p>
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>An indexed microbe reference in FASTA format (NCBI RefSeq/Genbank sequences)</li>
 *     <li>Downloaded NCBI RefSeq (and/or GenBank) catalog archive file(s)</li>
 *     <li>Downloaded NCBI taxonomy archive file</li>
 * </ul>
 *
 * <p>See argument documentation for information about where to download the archive files.</p>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A binary file containing reference taxonomy information</li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 * gatk PathSeqBuildReferenceTaxonomy \
 *   --reference microbe_reference.fasta \
 *   --output taxonomy.db \
 *   --refseq-catalog RefSeq-releaseXX.catalog.gz \
 *   --tax-dump taxdump.tar.gz \
 *   --min-non-virus-contig-length 2000
 * </pre>
 *
 * <h3>Notes</h3>
 *
 * <p>Often there are inconsistencies between the reference sequences, NCBI catalog, and taxonomy archive. To
 * minimize this issue, ensure that the input files are retrieved on the same date.</p>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@DocumentedFeature
@CommandLineProgramProperties(summary = "Build an annotated taxonomy datafile for a given microbe reference. The output file from this tool is required to run the PathSeq pipeline.",
        oneLineSummary = "Builds a taxonomy datafile of the microbe reference",
        programGroup = MetagenomicsProgramGroup.class)
public class PathSeqBuildReferenceTaxonomy extends CommandLineProgram {

    public static final String REFSEQ_CATALOG_LONG_NAME = "refseq-catalog";
    public static final String REFSEQ_CATALOG_SHORT_NAME = "RC";
    public static final String GENBANK_CATALOG_LONG_NAME = "genbank-catalog";
    public static final String GENBANK_CATALOG_SHORT_NAME = "GC";
    public static final String TAX_DUMP_LONG_NAME = "tax-dump";
    public static final String TAX_DUMP_SHORT_NAME = "TD";
    public static final String MIN_NON_VIRUS_CONTIG_LENGTH_LONG_NAME = "min-non-virus-contig-length";
    public static final String MIN_NON_VIRUS_CONTIG_LENGTH_SHORT_NAME = MIN_NON_VIRUS_CONTIG_LENGTH_LONG_NAME;

    @ArgumentCollection
    protected final ReferenceInputArgumentCollection referenceArguments = new RequiredReferenceInputArgumentCollection();

    @Argument(doc = "Local path for the output file. By convention, the extension should be \".db\"",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public String outputPath;

    @Argument(doc = "Local path to catalog file " +
            "(RefSeq-releaseXX.catalog.gz available at ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/)",
            fullName = REFSEQ_CATALOG_LONG_NAME,
            shortName = REFSEQ_CATALOG_SHORT_NAME,
            optional = true)
    public String refseqCatalogPath = null;

    /**
     * This may be supplied alone or in addition to the RefSeq catalog in the case that sequences from GenBank are
     * present in the reference.
     */
    @Argument(doc = "Local path to Genbank catalog file " +
            "(gbXXX.catalog.XXX.txt.gz at ftp://ftp.ncbi.nlm.nih.gov/genbank/catalog/)",
            fullName = GENBANK_CATALOG_LONG_NAME,
            shortName = GENBANK_CATALOG_SHORT_NAME,
            optional = true)
    public String genbankCatalogPath = null;

    @Argument(doc = "Local path to taxonomy dump tarball (taxdump.tar.gz available at ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)",
            fullName = TAX_DUMP_LONG_NAME,
            shortName = TAX_DUMP_SHORT_NAME)
    public String taxdumpPath;

    /**
     * Sequences from non-virus organisms less than this length will be filtered out such that any reads aligning to them will
     * be ignored. This is a quality control measure to remove shorter sequences from draft genomes that are likely to
     * contain sequencing artifacts such as cross-species contamination or sequencing adapters. Note this may
     * remove some bacteria plasmid sequences.
     */
    @Argument(doc = "Minimum reference contig length for non-viruses",
            fullName = MIN_NON_VIRUS_CONTIG_LENGTH_LONG_NAME,
            shortName = MIN_NON_VIRUS_CONTIG_LENGTH_SHORT_NAME,
            minValue = 0,
            minRecommendedValue = 500,
            maxRecommendedValue = 10000)
    public int minNonVirusContigLength = 0;

    @Override
    public Object doWork() {

        if (refseqCatalogPath == null && genbankCatalogPath == null) {
            throw new UserException.BadInput("At least one of --refseq-catalog or --genbank-catalog must be specified");
        }

        logger.info("Parsing reference and files... (this may take a few minutes)");
        final ReferenceDataSource reference = ReferenceDataSource.of(referenceArguments.getReferencePath());
        if (reference.getSequenceDictionary() == null) {
            throw new UserException.BadInput("Reference sequence dictionary not found. Please build one using CreateSequenceDictionary.");
        }
        final List<SAMSequenceRecord> referenceRecords = reference.getSequenceDictionary().getSequences();

        //Parse reference index, filling in data to taxIdToProperties and accessionToNameAndLength where possible
        final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties = new HashMap<>();
        final Map<String, Tuple2<String, Long>> accessionToNameAndLength = PSBuildReferenceTaxonomyUtils.parseReferenceRecords(referenceRecords, taxIdToProperties);

        //Parse RefSeq catalog to determine taxonomic ID's of accession keys in accessionToNameAndLength
        Set<String> accessionsNotFound = null;
        if (refseqCatalogPath != null) {
            try (final BufferedReader refseqCatalogStreamReader = PSBuildReferenceTaxonomyUtils.getBufferedReaderGz(refseqCatalogPath)) {
                accessionsNotFound = PSBuildReferenceTaxonomyUtils.parseCatalog(refseqCatalogStreamReader, accessionToNameAndLength, taxIdToProperties, false, null);
            } catch (IOException e) {
                throw new GATKException("Error reading RefSeq catalog", e);
            }
        }

        //Parse Genbank catalog
        if (genbankCatalogPath != null) {
            try (final BufferedReader genbankCatalogStreamReader = PSBuildReferenceTaxonomyUtils.getBufferedReaderGz(genbankCatalogPath)) {
                accessionsNotFound = PSBuildReferenceTaxonomyUtils.parseCatalog(genbankCatalogStreamReader, accessionToNameAndLength, taxIdToProperties, true, accessionsNotFound);
            } catch (IOException e) {
                throw new GATKException("Error reading GenBank catalog", e);
            }
        }
        if (accessionsNotFound != null && !accessionsNotFound.isEmpty()) {
            PSUtils.logItemizedWarning(logger, accessionsNotFound, "Did not find entries in the catalog for the following reference accessions");
        }

        //Get the scientific name of every taxonomic node (even the ones not in the reference)
        try (final BufferedReader namesStreamReader = PSBuildReferenceTaxonomyUtils.getBufferedReaderTarGz(taxdumpPath, "names.dmp")) {
            PSBuildReferenceTaxonomyUtils.parseNamesFile(namesStreamReader, taxIdToProperties);
        } catch (IOException e) {
            throw new GATKException("Error reading taxdump names files", e);
        }

        //Gets the taxonomic rank (e.g. family, order, genus, species, etc.) and parent of each node
        try (final BufferedReader nodesStreamReader = PSBuildReferenceTaxonomyUtils.getBufferedReaderTarGz(taxdumpPath, "nodes.dmp")) {
            final Collection<Integer> taxNotFound = PSBuildReferenceTaxonomyUtils.parseNodesFile(nodesStreamReader, taxIdToProperties);
            PSUtils.logItemizedWarning(logger, taxNotFound, "Did not find entry from reference sequence names or the names file for following some tax ID's. Setting name to tax_<tax ID>");
        } catch (IOException e) {
            throw new GATKException("Error reading taxdump names files", e);
        }

        //Build the taxonomic tree and a map from accession ID's to taxonomic ID's
        logger.info("Building taxonomic database...");
        final PSTree tree = PSBuildReferenceTaxonomyUtils.buildTaxonomicTree(taxIdToProperties);
        PSBuildReferenceTaxonomyUtils.removeUnusedTaxIds(taxIdToProperties, tree);
        final Map<String, Integer> accessionToTaxId = PSBuildReferenceTaxonomyUtils.buildAccessionToTaxIdMap(taxIdToProperties, tree, minNonVirusContigLength);

        //Write output
        PSBuildReferenceTaxonomyUtils.writeTaxonomyDatabase(outputPath, new PSTaxonomyDatabase(tree, accessionToTaxId));

        return null;
    }

}
