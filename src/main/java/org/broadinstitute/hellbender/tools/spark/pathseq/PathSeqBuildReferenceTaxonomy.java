package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.PathSeqProgramGroup;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import scala.Tuple2;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * Builds a minimal taxonomic database for a given pathogen reference. The database contains the taxonomic tree of all
 * organisms contained in the reference, and a mapping from each reference sequence to its corresponding NCBI
 * taxonomic id.
 * <p>
 * The database is built from a RefSeq and/or Genbank catalog file and the NCBI taxonomy dump.
 * <p>
 * The database is written to an output file, which is required by the ClassifyReads tool.
 */
@DocumentedFeature
@CommandLineProgramProperties(summary = "Builds a taxonomic database of the pathogen reference that " +
        "is required to run the scoring tool. User must supply a pathogen reference, NCBI catalog, and NCBI taxonomy " +
        "'taxdump' archive.",
        oneLineSummary = "Builds a taxonomy database of the pathogen reference",
        programGroup = PathSeqProgramGroup.class)
@BetaFeature
public class PathSeqBuildReferenceTaxonomy extends CommandLineProgram {

    @ArgumentCollection
    protected final ReferenceInputArgumentCollection referenceArguments = new RequiredReferenceInputArgumentCollection();
    @Argument(doc = "Local path for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public String outputPath;
    @Argument(doc = "Local path to catalog file " +
            "(RefSeq-releaseXX.catalog.gz available at ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/)",
            fullName = "refseqCatalogPath",
            optional = true)
    public String refseqCatalogPath = null;
    @Argument(doc = "Local path to Genbank catalog file " +
            "(gbXXX.catalog.XXX.txt.gz at ftp://ftp.ncbi.nlm.nih.gov/genbank/catalog/)",
            fullName = "genbankCatalogPath",
            optional = true)
    public String genbankCatalogPath = null;
    @Argument(doc = "Local path to taxonomy dump tarball (taxdump.tar.gz available at ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)",
            fullName = "taxdumpPath")
    public String taxdumpPath;
    @Argument(doc = "Minimum reference contig length for non-viruses",
            fullName = "minNonVirusContigLength")
    public int minNonVirusContigLength = 0;

    @Override
    public Object doWork() {

        if (refseqCatalogPath == null && genbankCatalogPath == null) {
            throw new UserException.BadInput("At least one of --refseqCatalogPath or --genbankCatalogPath must be specified");
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
