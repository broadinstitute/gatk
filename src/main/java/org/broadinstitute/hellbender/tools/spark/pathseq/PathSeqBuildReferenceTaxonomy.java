package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import scala.Tuple2;

import java.io.BufferedReader;
import java.util.*;

/**
 * Builds a minimal taxonomic database for a given pathogen reference. The database contains the taxonomic tree of all
 * organisms contained in the reference, and a mapping from each reference sequence to its corresponding NCBI
 * taxonomic id. This tool is designed to be run locally.
 * <p>
 * The database is built from a RefSeq and/or Genbank catalog file and the NCBI taxonomy dump, which are passed in as tool arguments.
 * <p>
 * The database is written to an output file, which is required by the ClassifyReads tool.
 */
@CommandLineProgramProperties(summary = "Builds PathSeq taxonomy database for a given pathogen reference",
        oneLineSummary = "PathSeq taxonomy database builder",
        programGroup = ReadProgramGroup.class)
public class PathSeqBuildReferenceTaxonomy extends GATKTool {

    @Argument(doc = "Local path for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    public String OUTPUT_PATH;
    @Argument(doc = "Local path to catalog file " +
            "(RefSeq-releaseXX.catalog.gz available at ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/)",
            fullName = "refseqCatalogPath",
        optional = true)
    public String REFSEQ_CATALOG_PATH = null;
    @Argument(doc = "Local path to Genbank catalog file " +
            "(gbXXX.catalog.XXX.txt.gz at ftp://ftp.ncbi.nlm.nih.gov/genbank/catalog/)",
            fullName = "genbankCatalogPath",
            optional = true)
    public String GENBANK_CATALOG_PATH = null;
    @Argument(doc = "Local path to taxonomy dump tarball (taxdump.tar.gz available at ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)",
            fullName = "taxdumpPath")
    public String TAXDUMP_PATH;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onStartup() {
        super.onStartup();

        logger.info("Parsing reference and files... (this may take a several minutes)");
        final List<SAMSequenceRecord> dictList = getReferenceDictionary().getSequences();

        final Map<String, PSTaxInfo> taxToInfo = new HashMap<>();
        final Map<String, Tuple2<String, Long>> accToRefInfo = PSTaxonomyBuilderUtils.parseReferenceRecords(dictList, taxToInfo);

        Set<String> accNotFound = null;
        if (REFSEQ_CATALOG_PATH != null) {
            final BufferedReader refseqCatalogStreamReader = PSTaxonomyBuilderUtils.getBufferedReaderGz(REFSEQ_CATALOG_PATH);
            accNotFound = PSTaxonomyBuilderUtils.parseCatalog(refseqCatalogStreamReader, accToRefInfo, taxToInfo, false, null);
            PSTaxonomyBuilderUtils.closeReader(refseqCatalogStreamReader);
        }

        if (GENBANK_CATALOG_PATH != null) {
            final BufferedReader genbankCatalogStreamReader = PSTaxonomyBuilderUtils.getBufferedReaderGz(GENBANK_CATALOG_PATH);
            accNotFound = PSTaxonomyBuilderUtils.parseCatalog(genbankCatalogStreamReader, accToRefInfo, taxToInfo, true, accNotFound);
            PSTaxonomyBuilderUtils.closeReader(genbankCatalogStreamReader);
        }
        if (accNotFound != null && !accNotFound.isEmpty()) {
            PSUtils.logItemizedWarning(logger, accNotFound, "Did not find entries in the catalog for the following reference accessions");
        }

        final BufferedReader namesStreamReader = PSTaxonomyBuilderUtils.getBufferedReaderTarGz(TAXDUMP_PATH, "names.dmp");
        PSTaxonomyBuilderUtils.parseNamesFile(namesStreamReader, taxToInfo);
        PSTaxonomyBuilderUtils.closeReader(namesStreamReader);

        final BufferedReader nodesStreamReader = PSTaxonomyBuilderUtils.getBufferedReaderTarGz(TAXDUMP_PATH, "nodes.dmp");
        final Collection<String> taxNotFound = PSTaxonomyBuilderUtils.parseNodesFile(nodesStreamReader, taxToInfo);
        PSUtils.logItemizedWarning(logger, taxNotFound, "Did not find entry from reference sequence names or the names file for following some tax ID's. Setting name to tax_<tax ID>");
        PSTaxonomyBuilderUtils.closeReader(nodesStreamReader);

        logger.info("Building taxonomic database...");
        final PSTree tree = PSTaxonomyBuilderUtils.buildTaxonomicTree(taxToInfo);
        final Map<String, String> refNameToTax = PSTaxonomyBuilderUtils.buildReferenceNameToTaxMap(taxToInfo);
        PSUtils.writeKryoTwo(OUTPUT_PATH, tree, refNameToTax);
    }

    @Override
    public void traverse() {
    }


}
