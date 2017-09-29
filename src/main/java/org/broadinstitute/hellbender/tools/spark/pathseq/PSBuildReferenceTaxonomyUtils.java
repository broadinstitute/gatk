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

    protected static final Logger logger = LogManager.getLogger(PSBuildReferenceTaxonomyUtils.class);
    private static final String VERTICAL_BAR_DELIMITER_REGEX = "\\s*\\|\\s*";
    /**
     * Build set of accessions contained in the reference.
     * Returns: a map from accession to the name and length of the record. If the sequence name contains the
     * taxonomic ID, it instead gets added to taxIdToProperties. Later we merge both results into taxIdToProperties.
     * Method: First, look for either "taxid|<taxid>|" or "ref|<accession>|" in the sequence name. If neither of
     * those are found, use the first word of the name as the accession.
     */
    protected static Map<String, Tuple2<String, Long>> parseReferenceRecords(final List<SAMSequenceRecord> dictionaryList,
                                                                             final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties) {

        final Map<String, Tuple2<String, Long>> accessionToNameAndLength = new HashMap<>();
        for (final SAMSequenceRecord record : dictionaryList) {
            final String recordName = record.getSequenceName();
            final long recordLength = record.getSequenceLength();
            final String[] tokens = recordName.split(VERTICAL_BAR_DELIMITER_REGEX);
            String recordAccession = null;
            int recordTaxId = PSTree.NULL_NODE;
            for (int i = 0; i < tokens.length - 1 && recordTaxId == PSTree.NULL_NODE; i++) {
                if (tokens[i].equals("ref")) {
                    recordAccession = tokens[i + 1];
                } else if (tokens[i].equals("taxid")) {
                    recordTaxId = parseTaxonId(tokens[i + 1]);
                }
            }
            if (recordTaxId == PSTree.NULL_NODE) {
                if (recordAccession == null) {
                    final String[] tokens2 = tokens[0].split(" "); //Default accession to first word in the name
                    recordAccession = tokens2[0];
                }
                accessionToNameAndLength.put(recordAccession, new Tuple2<>(recordName, recordLength));
            } else {
                addReferenceAccessionToTaxon(recordTaxId, recordName, recordLength, taxIdToProperties);
            }
        }
        return accessionToNameAndLength;
    }

    private static int parseTaxonId(final String taxonId) {
        try {
            return Integer.valueOf(taxonId);
        } catch (final NumberFormatException e) {
            throw new UserException.BadInput("Expected taxonomy ID to be an integer but found \"" + taxonId + "\"", e);
        }
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
                                              final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties,
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
                    final int taxId = parseTaxonId(tokens[taxIdColumnIndex]);
                    final String accession = tokens[accessionColumnIndex];
                    if (accessionToNameAndLength.containsKey(accession)) {
                        final Tuple2<String, Long> nameAndLength = accessionToNameAndLength.get(accession);
                        addReferenceAccessionToTaxon(taxId, nameAndLength._1, nameAndLength._2, taxIdToProperties);
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
    protected static void parseNamesFile(final BufferedReader reader, final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties) {
        try {
            String line;
            while ((line = reader.readLine()) != null) {
                //Split into columns delimited by <TAB>|<TAB>
                final String[] tokens = line.split(VERTICAL_BAR_DELIMITER_REGEX);
                if (tokens.length < 4) {
                    throw new UserException.BadInput("Expected at least 4 columns in tax dump names file but found " + tokens.length);
                }
                final String nameType = tokens[3];
                if (nameType.equals("scientific name")) {
                    final int taxId = parseTaxonId(tokens[0]);
                    final String name = tokens[1];
                    if (taxIdToProperties.containsKey(taxId)) {
                        taxIdToProperties.get(taxId).setName(name);
                    } else {
                        taxIdToProperties.put(taxId, new PSPathogenReferenceTaxonProperties(name));
                    }
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
    protected static Collection<Integer> parseNodesFile(final BufferedReader reader, final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties) {
        try {
            final Collection<Integer> taxIdsNotFound = new ArrayList<>();
            String line;
            while ((line = reader.readLine()) != null) {
                final String[] tokens = line.split(VERTICAL_BAR_DELIMITER_REGEX);
                if (tokens.length < 3) {
                    throw new UserException.BadInput("Expected at least 3 columns in tax dump nodes file but found " + tokens.length);
                }
                final int taxId = parseTaxonId(tokens[0]);
                final int parent = parseTaxonId(tokens[1]);
                final String rank = tokens[2];
                final PSPathogenReferenceTaxonProperties taxonProperties;
                if (taxIdToProperties.containsKey(taxId)) {
                    taxonProperties = taxIdToProperties.get(taxId);
                } else {
                    taxonProperties = new PSPathogenReferenceTaxonProperties("tax_" + taxId);
                    taxIdsNotFound.add(taxId);
                }
                taxonProperties.setRank(rank);
                if (taxId != PSTaxonomyConstants.ROOT_ID) { //keep root's parent set to null
                    taxonProperties.setParent(parent);
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
    private static void addReferenceAccessionToTaxon(final int taxId, final String accession, final long length, final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties) {
        taxIdToProperties.putIfAbsent(taxId, new PSPathogenReferenceTaxonProperties());
        taxIdToProperties.get(taxId).addAccession(accession, length);
    }

    /**
     * Removes nodes not in the tree from the tax_id-to-properties map
     */
    static void removeUnusedTaxIds(final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties,
                                   final PSTree tree) {
        taxIdToProperties.keySet().retainAll(tree.getNodeIDs());
    }

    /**
     * Create reference_name-to-taxid map (just an inversion on taxIdToProperties)
     */
    protected static Map<String, Integer> buildAccessionToTaxIdMap(final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties,
                                                                  final PSTree tree,
                                                                  final int minNonVirusContigLength) {
        final Map<String, Integer> accessionToTaxId = new HashMap<>();
        for (final int taxId : taxIdToProperties.keySet()) {
            final boolean isVirus = tree.getPathOf(taxId).contains(PSTaxonomyConstants.VIRUS_ID);
            final PSPathogenReferenceTaxonProperties taxonProperties = taxIdToProperties.get(taxId);
            for (final String name : taxonProperties.getAccessions()) {
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
    protected static PSTree buildTaxonomicTree(final Map<Integer, PSPathogenReferenceTaxonProperties> taxIdToProperties) {

        //Build tree of all taxa
        final PSTree tree = new PSTree(PSTaxonomyConstants.ROOT_ID);
        final Collection<Integer> invalidIds = new HashSet<>(taxIdToProperties.size());
        for (final int taxId : taxIdToProperties.keySet()) {
            if (taxId != PSTaxonomyConstants.ROOT_ID) {
                final PSPathogenReferenceTaxonProperties taxonProperties = taxIdToProperties.get(taxId);
                if (taxonProperties.getName() != null && taxonProperties.getParent() != PSTree.NULL_NODE && taxonProperties.getRank() != null) {
                    tree.addNode(taxId, taxonProperties.getName(), taxonProperties.getParent(), taxonProperties.getTotalLength(), taxonProperties.getRank());
                } else {
                    invalidIds.add(taxId);
                }
            }
        }
        PSUtils.logItemizedWarning(logger, invalidIds, "The following taxonomic IDs did not have name/taxonomy information (this may happen when the catalog and taxdump files are inconsistent)");

        final Set<Integer> unreachableNodes = tree.removeUnreachableNodes();
        if (!unreachableNodes.isEmpty()) {
            PSUtils.logItemizedWarning(logger, unreachableNodes, "Removed " + unreachableNodes.size() + " unreachable tree nodes");
        }

        tree.checkStructure();

        //Trim tree down to nodes corresponding only to reference taxa (and their ancestors)
        final Set<Integer> relevantNodes = new HashSet<>();
        for (final int taxonId : taxIdToProperties.keySet()) {
            if (!taxIdToProperties.get(taxonId).getAccessions().isEmpty() && tree.hasNode(taxonId)) {
                relevantNodes.addAll(tree.getPathOf(taxonId));
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
