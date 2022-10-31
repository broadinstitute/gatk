package org.broadinstitute.hellbender.tools.spark.pathseq;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.programgroups.MetagenomicsProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerShort;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Get PathSeq file information
 */
@DocumentedFeature
@CommandLineProgramProperties(summary = "Get PathSeq file information",
        oneLineSummary = "Get PathSeq file information",
        programGroup = MetagenomicsProgramGroup.class)
@BetaFeature
public final class PathSeqInfo extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger(PathSeqInfo.class);

    @Argument(doc = "Taxonomy file created with PathSeqBuildReferenceTaxonomy",
            fullName = "taxonomy-file",
            optional = true)
    public String taxonomyFile;

    @Argument(doc = "Species list output file (requires taxonomy-file)",
            fullName = "species-list",
            optional = true)
    public String speciesList;

    @Argument(doc = "K-mer file created with PathSeqBuildKmers",
            fullName = "kmer-file",
            optional = true)
    public String kmerFile;

    @Override
    protected Object doWork() {
        if (taxonomyFile == null && speciesList != null) {
            throw new UserException.BadInput("Argument \"taxonomy-file\" must be provided if \"species-list\" is specified");
        }
        if (taxonomyFile != null) {
            final PSTaxonomyDatabase taxonomyDatabase = PSScorer.readTaxonomyDatabase(taxonomyFile);
            writeTaxonomyTree(taxonomyDatabase, speciesList);
            printTaxonomyInfo(taxonomyDatabase);
        }
        if (kmerFile != null) {
            final PSKmerCollection kmerCollection = PSKmerUtils.readKmerFilter(kmerFile);
            printKmerInfo(kmerCollection);
        }
        return null;
    }

    private static void printKmerInfo(final PSKmerCollection kmerCollection) {
        Utils.nonNull(kmerCollection, "Cannot print info for null k-mer collection");
        logger.info("K-mer file properties:");
        logger.info("K-mer size: " + kmerCollection.kmerSize());
        logger.info("K-mer binary mask (2 bits per base): " + Long.toBinaryString(kmerCollection.getMask().getLong()));
        logger.info("Masked base positions: " + String.join(", ", getKmerMaskedBases(kmerCollection)));
        if (kmerCollection.getClass().isAssignableFrom(PSKmerSet.class)) {
            final PSKmerSet kmerSet = (PSKmerSet) kmerCollection;
            logger.info("Set type: Hopscotch set (hash table)");
            logger.info("Number of k-mers: " + kmerSet.setSize());
        } else if (kmerCollection.getClass().isAssignableFrom(PSKmerBloomFilter.class)) {
            final PSKmerBloomFilter kmerBloomFilter = (PSKmerBloomFilter) kmerCollection;
            logger.info("Set type: Bloom filter");
            logger.info("False positive probability: " + kmerBloomFilter.getFalsePositiveProbability());
        } else {
            logger.warn("Could not determine the class of the k-mer file");
        }
    }

    private static List<String> getKmerMaskedBases(final PSKmerCollection kmerCollection) {
        Utils.nonNull(kmerCollection, "Cannot get masked bases for null k-mer collection");
        final SVKmerShort kmerMask = kmerCollection.getMask();
        if (kmerMask == null) {
            throw new IllegalStateException(PSKmerCollection.class.getSimpleName() + " returned null mask");
        }
        long kmerMaskLong = kmerMask.getLong();
        final List<String> maskedBaseList = new ArrayList<>(32);
        for (int i = Long.SIZE - 2; i > 0; i -= 2) {
            if ((kmerMaskLong & 0x3L) == 0) {
                maskedBaseList.add(String.valueOf(i / 2));
            }
            kmerMaskLong >>= 2;
        }
        return maskedBaseList;
    }

    private static void writeTaxonomyTree(final PSTaxonomyDatabase taxonomyDatabase, final String speciesListPath) {
        Utils.nonNull(taxonomyDatabase, "Cannot print species for null taxonomy database");
        if (speciesListPath == null) return;
        final PSTree tree = taxonomyDatabase.tree;
        final Map<Integer, List<String>> taxIdToAccessionMap = new HashMap<>();
        for (final Map.Entry<String, Integer> entry : taxonomyDatabase.accessionToTaxId.entrySet()) {
            taxIdToAccessionMap.putIfAbsent(entry.getValue(), new ArrayList<>());
            taxIdToAccessionMap.get(entry.getValue()).add(entry.getKey());
        }
        try (final PrintStream printStream = new PrintStream(BucketUtils.createFile(speciesListPath))) {
            final List<Integer> list = new ArrayList<>(tree.getNodeIDs());
            Collections.sort(list);
            printStream.println("#taxonomy_id\tname\trank\tparent_id\tpath\tref_length\taccessions");
            for (final Integer nodeId : list) {
                final List<String> accessions = taxIdToAccessionMap.getOrDefault(nodeId, Collections.emptyList());
                final List<String> path = tree.getPathOf(nodeId).stream().map(String::valueOf).collect(Collectors.toList());
                printStream.println(nodeId + "\t" + tree.getNameOf(nodeId) + "\t" + tree.getRankOf(nodeId)
                        + "\t" + tree.getParentOf(nodeId) + "\t" + String.join(",", path)
                        + "\t" + tree.getLengthOf(nodeId) + "\t" + String.join(",", accessions));
            }
        }
    }

    private static void printTaxonomyInfo(final PSTaxonomyDatabase taxonomyDatabase) {
        Utils.nonNull(taxonomyDatabase, "Cannot print info for null taxonomy database");
        final PSTree tree = taxonomyDatabase.tree;
        logger.info("Taxonomy file properties:");
        logger.info("Number of contig accessions: " + taxonomyDatabase.accessionToTaxId.size());
        logger.info("Number of taxa: " + tree.getNodeIDs().size());
        logger.info("Total number of organisms: " + getNodesWithReference(tree).size());
        logger.info("Total length of references: " + getTaxonomyReferenceLength(tree));
        final List<Integer> kingdomNodeIds = getKingdomNodes(tree);
        for (final Integer nodeId : kingdomNodeIds) {
            final String kingdomName = tree.getNameOf(nodeId);
            final List<Integer> descendentsWithReferences = getDescendantNodesWithReference(tree, nodeId);
            final int numOrganisms = descendentsWithReferences.size();
            final long referenceLength = descendentsWithReferences.stream().map(tree::getLengthOf).mapToLong(Long::longValue).sum();
            logger.info("  " + kingdomName);
            logger.info("    Organisms: " + numOrganisms);
            logger.info("    Reference length: " + referenceLength + " bp");
        }
    }

    private static long getTaxonomyReferenceLength(final PSTree tree) {
        long totalLength = 0;
        for (final Integer nodeId : tree.getNodeIDs()) {
            final long nodeLength = tree.getLengthOf(nodeId);
            if (tree.getLengthOf(nodeId) > 0) {
                totalLength += nodeLength;
            }
        }
        return totalLength;
    }

    private static List<Integer> getKingdomNodes(final PSTree tree) {
        final List<Integer> list = new ArrayList<>();
        for (final Integer nodeId : tree.getNodeIDs()) {
            final String rank = tree.getRankOf(nodeId);
            if (rank.equals(PSTaxonomyConstants.KINGDOM_RANK_NAME) || rank.equals(PSTaxonomyConstants.SUPERKINGDOM_RANK_NAME)) {
                list.add(nodeId);
            }
        }
        return list;
    }

    private static List<Integer> getDescendantNodesWithReference(final PSTree tree, final Integer parentId) {
        final List<Integer> list = new ArrayList<>();
        final Queue<Integer> queue = new PriorityQueue<>();
        queue.add(parentId);
        while (!queue.isEmpty()) {
            final Integer nodeId = queue.poll();
            final Collection<Integer> children = tree.getChildrenOf(nodeId);
            queue.addAll(children);
            if (tree.getLengthOf(nodeId) > 0) {
                list.add(nodeId);
            }
        }
        return list;
    }

    private static List<Integer> getNodesWithReference(final PSTree tree) {
        final List<Integer> list = new ArrayList<>(tree.getNodeIDs().size());
        for (final Integer nodeId : tree.getNodeIDs()) {
            if (tree.getLengthOf(nodeId) > 0) {
                list.add(nodeId);
            }
        }
        return list;
    }
}
