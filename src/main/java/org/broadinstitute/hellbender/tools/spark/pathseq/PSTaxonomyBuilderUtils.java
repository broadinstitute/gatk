package org.broadinstitute.hellbender.tools.spark.pathseq;

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

public final class PSTaxonomyBuilderUtils {

    protected final static Logger logger = LogManager.getLogger(org.broadinstitute.hellbender.tools.spark.pathseq.PSTaxonomyBuilderUtils.class);
    private final static String ROOT_ID = "1"; //NCBI root node taxonomic id

    /**
     * Build set of accessions contained in the reference.
     * Returns: a map from accession to the name and length of the record. If the sequence name contains the
     * taxonomic ID, it instead gets added to taxToInfo. Later we merge both results into taxToInfo.
     * Method: First, look for either "taxid|<taxid>|" or "ref|<accession>|" in the sequence name. If neither of
     * those are found, use the first word of the name as the accession.
     */
    protected static Map<String, Tuple2<String, Long>> parseReferenceRecords(final List<SAMSequenceRecord> dictList,
                                                                             final Map<String, PSTaxInfo> taxToInfo) {

        final Map<String, Tuple2<String, Long>> accToRefInfo = new HashMap<>();
        for (final SAMSequenceRecord entry : dictList) {
            final String name = entry.getSequenceName();
            final long length = entry.getSequenceLength();
            final String[] tokens = name.split("[|]");
            String acc = null;
            String tax = null;
            for (int i = 0; i < tokens.length - 1 && tax == null; i++) {
                if (tokens[i].equals("ref")) {
                    acc = tokens[i + 1];
                } else if (tokens[i].equals("taxid")) {
                    tax = tokens[i + 1];
                }
            }
            if (tax == null) {
                if (acc == null) {
                    final String[] tokens_2 = tokens[0].split(" "); //Default accession to first word in the name
                    acc = tokens_2[0];
                }
                accToRefInfo.put(acc, new Tuple2<>(name, length));
            } else {
                addRefNameToTaxInfo(tax, name, length, taxToInfo);
            }
        }
        return accToRefInfo;
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
                                              final Map<String, Tuple2<String, Long>> accToRefInfo,
                                              final Map<String, PSTaxInfo> taxToInfo,
                                              final boolean bGenBank, Set<String> accNotFound) {
        try {
            String line;
            int taxidIndex, accIndex;
            if (bGenBank) {
                taxidIndex = 6;
                accIndex = 1;
            } else {
                taxidIndex = 0;
                accIndex = 2;
            }
            if (accNotFound == null) {
                accNotFound = new HashSet<>(accToRefInfo.keySet());
            }
            final int minColumns = Math.max(taxidIndex, accIndex) + 1;
            long lineNo = 1;
            while ((line = reader.readLine()) != null && !line.isEmpty()) {
                final String[] tokens = line.trim().split("\t", minColumns + 1);
                if (tokens.length >= minColumns) {
                    final String taxid = tokens[taxidIndex];
                    final String acc = tokens[accIndex];
                    if (accToRefInfo.containsKey(acc)) {
                        final Tuple2<String, Long> acc_info = accToRefInfo.get(acc);
                        addRefNameToTaxInfo(taxid, acc_info._1, acc_info._2, taxToInfo);
                        accNotFound.remove(acc);
                    }
                } else {
                    throw new UserException.BadInput("Expected at least " + minColumns + " tab-delimited columns in " +
                            "GenBank catalog file, but only found " + tokens.length + " on line " + lineNo);
                }
                lineNo++;
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile("Error reading from catalog file", e);
        }
        return accNotFound;
    }

    /**
     * Parses scientific name of each taxon
     */
    protected static void parseNamesFile(final BufferedReader reader, final Map<String, PSTaxInfo> taxToInfo) {
        try {
            String line;
            while ((line = reader.readLine()) != null) {
                final String[] tokens = line.split("\t\\|\t");
                if (tokens.length != 4) {
                    throw new UserException.BadInput("Expected 4 columns in tax dump names file but found " + tokens.length);
                }
                final String name_type = tokens[3];
                if (name_type.equals("scientific name\t|")) {
                    final String tax = tokens[0];
                    final String name = tokens[1];
                    final PSTaxInfo info = taxToInfo.containsKey(tax) ? taxToInfo.get(tax) : new PSTaxInfo();
                    info.name = name;
                    taxToInfo.put(tax, info);
                }
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile("Error reading from taxonomy dump names file", e);
        }
    }

    /**
     * Gets the rank and parent of each taxon.
     * Returns a Collection of tax ID's found in the nodes file that are not in taxToInfo (i.e. were not found in
     * a reference sequence name using the taxid|\<taxid\> tag or the catalog file).
     */
    protected static Collection<String> parseNodesFile(final BufferedReader reader, final Map<String, PSTaxInfo> taxToInfo) {
        try {
            final Collection<String> taxNotFound = new ArrayList<>();
            String line;
            while ((line = reader.readLine()) != null) {
                final String[] tokens = line.split("\t\\|\t");
                if (tokens.length != 13) {
                    throw new UserException.BadInput("Expected 13 columns in tax dump nodes file but found " + tokens.length);
                }
                final String taxid = tokens[0];
                final String parent = tokens[1];
                final String rank = tokens[2];
                PSTaxInfo info;
                if (taxToInfo.containsKey(taxid)) {
                    info = taxToInfo.get(taxid);
                } else {
                    info = new PSTaxInfo();
                    info.name = "tax_" + taxid;
                    taxNotFound.add(taxid);
                }
                info.rank = rank;
                if (!taxid.equals(ROOT_ID)) { //keep root's parent set to null
                    info.parent_tax = parent;
                }
                taxToInfo.put(taxid, info);
            }
            return taxNotFound;
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile("Error reading from taxonomy dump nodes file", e);
        }
    }

    /**
     * Helper function for building the map from tax id to reference contig accession
     */
    private static void addRefNameToTaxInfo(final String key, final String name, final long length, final Map<String, PSTaxInfo> taxToInfo) {
        PSTaxInfo info;
        if (taxToInfo.containsKey(key)) {
            info = taxToInfo.get(key);
        } else {
            info = new PSTaxInfo();
        }
        info.ref_names.add(name);
        info.length += length;
        taxToInfo.put(key, info);
    }

    /**
     * Create reference_name-to-taxid map (just an inversion on taxToInfo)
     */
    protected static Map<String, String> buildReferenceNameToTaxMap(final Map<String, PSTaxInfo> taxToInfo) {
        final Map<String, String> refNameToTax = new HashMap<>();
        for (final String tax : taxToInfo.keySet()) {
            final PSTaxInfo info = taxToInfo.get(tax);
            for (String name : info.ref_names) {
                refNameToTax.put(name, tax);
            }
        }
        return refNameToTax;
    }

    /**
     * Returns a PSTree representing a reduced taxonomic tree containing only taxa present in the reference
     */
    protected static PSTree buildTaxonomicTree(final Map<String, PSTaxInfo> taxToInfo) {

        //Build tree of all taxa
        final PSTree tree = new PSTree(ROOT_ID);
        final Collection<String> invalidIDs = new ArrayList<>();
        for (final String taxID : taxToInfo.keySet()) {
            if (!taxID.equals(ROOT_ID)) {
                final PSTaxInfo info = taxToInfo.get(taxID);
                if (info.name != null && info.parent_tax != null && info.rank != null) {
                    tree.addNode(taxID, info.name, info.parent_tax, info.length, info.rank);
                } else {
                    invalidIDs.add(taxID);
                }
            }
        }
        PSUtils.logItemizedWarning(logger, invalidIDs, "The following taxonomic IDs did not have name/taxonomy information (this may happen when the catalog and taxdump files are inconsistent)");

        final Set<String> unreachable = tree.removeUnreachableNodes();
        if (!unreachable.isEmpty()) {
            PSUtils.logItemizedWarning(logger, unreachable, "Removed " + unreachable.size() + " unreachable tree nodes");
        }

        tree.checkStructure();

        //Trim tree down to nodes corresponding only to reference taxa (and their ancestors)
        final Set<String> relevantNodes = new HashSet<>();
        for (final String tax : taxToInfo.keySet()) {
            if (!taxToInfo.get(tax).ref_names.isEmpty()) {
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

    public static BufferedReader getBufferedReaderGz(final String path) {
        try {
            return new BufferedReader(IOUtils.makeReaderMaybeGzipped(new File(path)));
        } catch (final IOException e) {
            throw new UserException.BadInput("Could not open file " + path);
        }
    }

    public static BufferedReader getBufferedReaderTarGz(final String tarPath, final String filename) {
        try {
            InputStream result = null;
            final TarArchiveInputStream tarStream = new TarArchiveInputStream(new GZIPInputStream(new FileInputStream(tarPath)));
            TarArchiveEntry entry = tarStream.getNextTarEntry();
            while (entry != null) {
                if (entry.getName().equals(filename)) {
                    result = tarStream;
                    break;
                }
                entry = tarStream.getNextTarEntry();
            }
            if (result == null) {
                throw new UserException.BadInput("Could not find file " + filename + " in tarball " + tarPath);
            }
            return new BufferedReader(new InputStreamReader(result));
        } catch (final IOException e) {
            throw new UserException.BadInput("Stream error for compressed tarball file " + filename + " in " + tarPath);
        }
    }

    public static void closeReader(final BufferedReader r) {
        try {
            if (r != null) {
                r.close();
            }
        } catch (final IOException e) {
            throw new RuntimeException("Could not close file reader", e);
        }
    }
}
