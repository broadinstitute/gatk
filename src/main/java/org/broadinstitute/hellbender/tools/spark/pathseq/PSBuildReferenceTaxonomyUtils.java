package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

public final class PSBuildReferenceTaxonomyUtils {

    protected final static Logger logger = LogManager.getLogger(PSBuildReferenceTaxonomyUtils.class);
    private final static String ROOT_ID = "1"; //NCBI root node taxonomic id
    final static String VIRUS_ID = "10239";

    /**
     * Build set of accessions contained in the reference.
     * Returns: a map from accession to the name and length of the record. If the sequence name contains the
     * taxonomic ID, it instead gets added to taxIdToProperties. Later we merge both results into taxIdToProperties.
     * Method: First, look for either "taxid|<taxid>|" or "ref|<accession>|" in the sequence name. If neither of
     * those are found, use the first word of the name as the accession.
     */
    protected static Map<String, Tuple2<String, Long>> parseReferenceRecords(final List<SAMSequenceRecord> dictionaryList,
                                                                             final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties) {

        final Map<String, Tuple2<String, Long>> accessionToNameAndLength = new HashMap<>();
        for (final SAMSequenceRecord record : dictionaryList) {
            final String recordName = record.getSequenceName();
            final long recordLength = record.getSequenceLength();
            final String[] tokens = recordName.split("[|]");
            String recordAccession = null;
            String recordTaxId = null;
            for (int i = 0; i < tokens.length - 1 && recordTaxId == null; i++) {
                if (tokens[i].equals("ref")) {
                    recordAccession = tokens[i + 1];
                } else if (tokens[i].equals("taxid")) {
                    recordTaxId = tokens[i + 1];
                }
            }
            if (recordTaxId == null) {
                if (recordAccession == null) {
                    final String[] tokens2 = tokens[0].split(" "); //Default accession to first word in the name
                    recordAccession = tokens2[0];
                }
                accessionToNameAndLength.put(recordAccession, new Tuple2<>(recordName, recordLength));
            } else {
                updateAccessionTaxonProperties(recordTaxId, recordName, recordLength, taxIdToProperties);
            }
        }
        return accessionToNameAndLength;
    }

    /**
     * Helper classes for defining RefSeq and GenBank catalog formats. Columns should be given as 0-based indices.
     */
    private interface AccessionCatalogFormat {
        int getTaxIdColumn();
        int getAccessionColumn();
    }

    private static final class RefSeqCatalogFormat implements AccessionCatalogFormat {
        private static final int TAX_ID_COLUMN = 0;
        private static final int ACCESSION_COLUMN = 2;
        public int getTaxIdColumn() {
            return TAX_ID_COLUMN;
        }
        public int getAccessionColumn() {
            return ACCESSION_COLUMN;
        }
    }

    private static final class GenBankCatalogFormat implements AccessionCatalogFormat {
        private static final int TAX_ID_COLUMN = 6;
        private static final int ACCESSION_COLUMN = 1;
        public int getTaxIdColumn() {
            return TAX_ID_COLUMN;
        }
        public int getAccessionColumn() {
            return ACCESSION_COLUMN;
        }
    }

    /**
     * Builds maps of reference contig accessions to their taxonomic ids and vice versa.
     * Input can be a RefSeq or Genbank catalog file. accNotFound is an initial list of
     * accessions from the reference that have not been successfully looked up; if null,
     * will be initialized to the accToRefInfo key set by default.
     * <p>
     * Returns a collection of reference accessions that could not be found, if any.
     */
    protected static Set<String> parseCatalog(final BufferedReader reader,
                                              final Map<String, Tuple2<String, Long>> accessionToNameAndLength,
                                              final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties,
                                              final boolean bGenBank,
                                              final Set<String> accessionsNotFoundIn) {

        final Set<String> accessionsNotFoundOut;
        try {
            String line;
            final AccessionCatalogFormat catalogFormat = bGenBank ? new GenBankCatalogFormat() : new RefSeqCatalogFormat();
            final int taxIdColumnIndex = catalogFormat.getTaxIdColumn();
            final int accessionColumnIndex = catalogFormat.getAccessionColumn();
            if (accessionsNotFoundIn == null) {
                //If accessionsNotFoundIn is null, this is the first call to parseCatalog, so initialize the set to all accessions
                accessionsNotFoundOut = new HashSet<>(accessionToNameAndLength.keySet());
            } else {
                //Otherwise this is a subsequent call and we continue to look for any remaining accessions
                accessionsNotFoundOut = new HashSet<>(accessionsNotFoundIn);
            }
            final int minColumns = Math.max(taxIdColumnIndex, accessionColumnIndex) + 1;
            long lineNumber = 1;
            while ((line = reader.readLine()) != null && !line.isEmpty()) {
                final String[] tokens = line.trim().split("\t", minColumns + 1);
                if (tokens.length >= minColumns) {
                    final String taxId = tokens[taxIdColumnIndex];
                    final String accession = tokens[accessionColumnIndex];
                    if (accessionToNameAndLength.containsKey(accession)) {
                        final Tuple2<String, Long> nameAndLength = accessionToNameAndLength.get(accession);
                        updateAccessionTaxonProperties(taxId, nameAndLength._1, nameAndLength._2, taxIdToProperties);
                        accessionsNotFoundOut.remove(accession);
                    }
                } else {
                    throw new UserException.BadInput("Expected at least " + minColumns + " tab-delimited columns in " +
                            "GenBank catalog file, but only found " + tokens.length + " on line " + lineNumber);
                }
                lineNumber++;
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile("Error reading from catalog file", e);
        }
        return accessionsNotFoundOut;
    }

    /**
     * Parses scientific name of each taxon and puts it in taxIdToProperties
     */
    protected static void parseNamesFile(final BufferedReader reader, final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties) {
        try {
            String line;
            while ((line = reader.readLine()) != null) {
                //Split into columns delimited by <TAB>|<TAB>. Note the "\\|" is required to match "|" otherwise Java
                // thinks it's an escape character and throws and error
                final String[] tokens = line.split("\t\\|\t");
                if (tokens.length != 4) {
                    throw new UserException.BadInput("Expected 4 columns in tax dump names file but found " + tokens.length);
                }
                final String nameType = tokens[3];
                if (nameType.equals("scientific name\t|")) {
                    final String taxId = tokens[0];
                    final String name = tokens[1];
                    final PSPathogenReferenceTaxonProperties taxonProperties = taxIdToProperties.containsKey(taxId) ? taxIdToProperties.get(taxId) : new PSPathogenReferenceTaxonProperties();
                    taxonProperties.name = name;
                    taxIdToProperties.put(taxId, taxonProperties);
                }
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile("Error reading from taxonomy dump names file", e);
        }
    }

    /**
     * Gets the rank and parent of each taxon.
     * Returns a Collection of tax ID's found in the nodes file that are not in taxIdToProperties (i.e. were not found in
     * a reference sequence name using the taxid|\<taxid\> tag or the catalog file).
     */
    protected static Collection<String> parseNodesFile(final BufferedReader reader, final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties) {
        try {
            final Collection<String> taxIdsNotFound = new ArrayList<>();
            String line;
            while ((line = reader.readLine()) != null) {
                final String[] tokens = line.split("\t\\|\t");
                if (tokens.length != 13) {
                    throw new UserException.BadInput("Expected 13 columns in tax dump nodes file but found " + tokens.length);
                }
                final String taxId = tokens[0];
                final String parent = tokens[1];
                final String rank = tokens[2];
                final PSPathogenReferenceTaxonProperties taxonProperties;
                if (taxIdToProperties.containsKey(taxId)) {
                    taxonProperties = taxIdToProperties.get(taxId);
                } else {
                    taxonProperties = new PSPathogenReferenceTaxonProperties();
                    taxonProperties.name = "tax_" + taxId;
                    taxIdsNotFound.add(taxId);
                }
                taxonProperties.rank = rank;
                if (!taxId.equals(ROOT_ID)) { //keep root's parent set to null
                    taxonProperties.parentTaxId = parent;
                }
                taxIdToProperties.put(taxId, taxonProperties);
            }
            return taxIdsNotFound;
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile("Error reading from taxonomy dump nodes file", e);
        }
    }

    /**
     * Helper function for building the map from tax id to reference contig accession
     */
    private static void updateAccessionTaxonProperties(final String accession, final String name, final long length, final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties) {
        final PSPathogenReferenceTaxonProperties taxonProperties;
        if (taxIdToProperties.containsKey(accession)) {
            taxonProperties = taxIdToProperties.get(accession);
        } else {
            taxonProperties = new PSPathogenReferenceTaxonProperties();
        }
        taxonProperties.addAccession(name, length);
        taxIdToProperties.put(accession, taxonProperties);
    }

    /**
     * Removes nodes not in the tree from the tax_id-to-properties map
     */
    static void removeUnusedTaxIds(final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties,
                                   final PSTree tree) {
        final Set<String> unusedNodes = new HashSet<>(taxIdToProperties.keySet());
        unusedNodes.removeAll(tree.getNodeIDs());
        taxIdToProperties.keySet().removeAll(unusedNodes);
    }

    /**
     * Create reference_name-to-taxid map (just an inversion on taxIdToProperties)
     */
    protected static Map<String, String> buildAccessionToTaxIdMap(final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties,
                                                                  final PSTree tree,
                                                                  final int minNonVirusContigLength) {
        final Map<String, String> accessionToTaxId = new HashMap<>();
        for (final String taxId : taxIdToProperties.keySet()) {
            final boolean isVirus = tree.getPathOf(taxId).contains(VIRUS_ID);
            final PSPathogenReferenceTaxonProperties taxonProperties = taxIdToProperties.get(taxId);
            for (String name : taxonProperties.getAccessions()) {
                if (isVirus || taxonProperties.getAccessionLength(name) >= minNonVirusContigLength) {
                    accessionToTaxId.put(name, taxId);
                }
            }
        }
        return accessionToTaxId;
    }

    /**
     * Returns a PSTree representing a reduced taxonomic tree containing only taxa present in the reference
     */
    protected static PSTree buildTaxonomicTree(final Map<String, PSPathogenReferenceTaxonProperties> taxIdToProperties) {

        //Build tree of all taxa
        final PSTree tree = new PSTree(ROOT_ID);
        final Collection<String> invalidIds = new HashSet<>(taxIdToProperties.size());
        for (final String taxId : taxIdToProperties.keySet()) {
            if (!taxId.equals(ROOT_ID)) {
                final PSPathogenReferenceTaxonProperties taxonProperties = taxIdToProperties.get(taxId);
                if (taxonProperties.name != null && taxonProperties.parentTaxId != null && taxonProperties.rank != null) {
                    tree.addNode(taxId, taxonProperties.name, taxonProperties.parentTaxId, taxonProperties.getTotalLength(), taxonProperties.rank);
                } else {
                    invalidIds.add(taxId);
                }
            }
        }
        PSUtils.logItemizedWarning(logger, invalidIds, "The following taxonomic IDs did not have name/taxonomy information (this may happen when the catalog and taxdump files are inconsistent)");

        final Set<String> unreachableNodes = tree.removeUnreachableNodes();
        if (!unreachableNodes.isEmpty()) {
            PSUtils.logItemizedWarning(logger, unreachableNodes, "Removed " + unreachableNodes.size() + " unreachable tree nodes");
        }

        tree.checkStructure();

        //Trim tree down to nodes corresponding only to reference taxa (and their ancestors)
        final Set<String> relevantNodes = new HashSet<>();
        for (final String tax : taxIdToProperties.keySet()) {
            if (!taxIdToProperties.get(tax).getAccessions().isEmpty()) {
                if (tree.hasNode(tax)) {
                    relevantNodes.addAll(tree.getPathOf(tax));
                }
            }
        }
        if (relevantNodes.isEmpty()) {
            throw new UserException.BadInput("Did not find any taxa corresponding to reference sequence names.\n\n"
                    + "Check that reference names follow one of the required formats:\n\n"
                    + "\t...|ref|<accession.version>|...\n"
                    + "\t...|taxid|<taxonomy_id>|...\n"
                    + "\t<accession.version><mask>...");
        }
        tree.retainNodes(relevantNodes);

        return tree;
    }

    /**
     * Gets a buffered reader for a gzipped file
     * @param path   File path
     * @return  Reader for the file
     */
    public static BufferedReader getBufferedReaderGz(final String path) {
        try {
            return new BufferedReader(IOUtils.makeReaderMaybeGzipped(new File(path)));
        } catch (final IOException e) {
            throw new UserException.BadInput("Could not open file " + path, e);
        }
    }

    /**
     * Gets a Reader for a file in a gzipped tarball
     * @param tarPath   Path to the tarball
     * @param fileName   File within the tarball
     * @return   The file's reader
     */
    public static BufferedReader getBufferedReaderTarGz(final String tarPath, final String fileName) {
        try {
            InputStream result = null;
            final TarArchiveInputStream tarStream = new TarArchiveInputStream(new GZIPInputStream(new FileInputStream(tarPath)));
            TarArchiveEntry entry = tarStream.getNextTarEntry();
            while (entry != null) {
                if (entry.getName().equals(fileName)) {
                    result = tarStream;
                    break;
                }
                entry = tarStream.getNextTarEntry();
            }
            if (result == null) {
                throw new UserException.BadInput("Could not find file " + fileName + " in tarball " + tarPath);
            }
            return new BufferedReader(new InputStreamReader(result));
        } catch (final IOException e) {
            throw new UserException.BadInput("Could not open compressed tarball file " + fileName + " in " + tarPath, e);
        }
    }

    /**
     * Writes objects using Kryo to specified local file path.
     * NOTE: using setReferences(false), which must also be set when reading the file. Does not work with nested
     * objects that reference its parent.
     */
    public static void writeTaxonomyDatabase(final String filePath, final PSTaxonomyDatabase taxonomyDatabase) {
        try {
            final Kryo kryo = new Kryo();
            kryo.setReferences(false);
            Output output = new Output(new FileOutputStream(filePath));
            kryo.writeObject(output, taxonomyDatabase);
            output.close();
        } catch (final FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not serialize objects to file", e);
        }
    }
}
