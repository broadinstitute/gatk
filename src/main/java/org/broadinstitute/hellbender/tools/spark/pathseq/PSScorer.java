package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import htsjdk.samtools.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public final class PSScorer {

    public static final String HITS_TAG = "YP";
    //Convert scores to per million reference base pairs
    public static final double SCORE_GENOME_LENGTH_UNITS = 1e6;
    private static final Logger logger = LogManager.getLogger(PSScorer.class);

    private final PSScoreArgumentCollection scoreArgs;

    public PSScorer(final PSScoreArgumentCollection scoreArgs) {
        this.scoreArgs = scoreArgs;
    }

    public JavaRDD<GATKRead> scoreReads(final JavaSparkContext ctx,
                                        final JavaRDD<GATKRead> pairedReads,
                                        final JavaRDD<GATKRead> unpairedReads,
                                        final SAMFileHeader header) {

        //Group reads into pairs
        final JavaRDD<Iterable<GATKRead>> groupedReads = groupReadsIntoPairs(pairedReads,
                unpairedReads, scoreArgs.readsPerPartition);

        //Load taxonomy database, created by running PathSeqBuildReferenceTaxonomy with this reference
        final PSTaxonomyDatabase taxDB = readTaxonomyDatabase(scoreArgs.taxonomyDatabasePath);
        final Broadcast<PSTaxonomyDatabase> taxonomyDatabaseBroadcast = ctx.broadcast(taxDB);

        //Check header against database
        if (scoreArgs.headerWarningFile != null) {
            writeMissingReferenceAccessions(scoreArgs.headerWarningFile, header, taxDB, logger);
        }

        //Determine which alignments are valid hits and return their tax IDs in PSPathogenAlignmentHit
        //Also adds pathseq tags containing the hit IDs to the reads
        final JavaRDD<Tuple2<Iterable<GATKRead>, PSPathogenAlignmentHit>> readHits = mapGroupedReadsToTax(groupedReads,
                scoreArgs.minCoverage, scoreArgs.minIdentity, taxonomyDatabaseBroadcast);

        //Get the original reads, now with their pathseq hit tags set
        final JavaRDD<GATKRead> readsFinal = flattenIterableKeys(readHits);

        //Compute taxonomic scores from the alignment hits
        final JavaRDD<PSPathogenAlignmentHit> alignmentHits = readHits.map(Tuple2::_2);
        final JavaPairRDD<String, PSPathogenTaxonScore> taxScoresRdd = alignmentHits.mapPartitionsToPair(iter -> computeTaxScores(iter, taxonomyDatabaseBroadcast.value()));

        //Reduce scores by taxon and compute normalized scores
        Map<String, PSPathogenTaxonScore> taxScoresMap = new HashMap<>(taxScoresRdd.reduceByKey(PSPathogenTaxonScore::add).collectAsMap());
        taxScoresMap = computeNormalizedScores(taxScoresMap, taxDB.tree);

        //Write scores to file
        writeScoresFile(taxScoresMap, taxDB.tree, scoreArgs.scoresPath);

        return readsFinal;
    }

    /**
     * Collects the second elements of all tuples in an RDD
     */
    static <T, C> Iterable<C> collectValues(final JavaRDD<Tuple2<T, C>> tupleRdd) {
        return tupleRdd.map(tuple -> tuple._2).collect();
    }

    /**
     * Flattens RDD of tuples, in which the first elements are Iterable, to an RDD of the items in the Iterables
     */
    static <T, C> JavaRDD<T> flattenIterableKeys(final JavaRDD<Tuple2<Iterable<T>, C>> tupleRdd) {
        return tupleRdd.flatMap(tuple -> tuple._1.iterator());
    }

    /**
     * Moves reads from the same read template into an Iterable.
     * Paired reads must be queryname-sorted, and no pair of reads can be split across partitions.
     */
    static JavaRDD<Iterable<GATKRead>> groupReadsIntoPairs(final JavaRDD<GATKRead> pairedReads,
                                                           final JavaRDD<GATKRead> unpairedReads,
                                                           final int readsPerPartitionGuess) {
        JavaRDD<Iterable<GATKRead>> groupedReads;
        if (pairedReads != null) {
            groupedReads = pairedReads.mapPartitions(iter -> groupPairedReadsPartition(iter, readsPerPartitionGuess));
            if (unpairedReads != null) {
                groupedReads = groupedReads.union(unpairedReads.map(Collections::singletonList));
            }
        } else if (unpairedReads != null) {
            groupedReads = unpairedReads.map(Collections::singletonList);
        } else {
            throw new UserException.BadInput("No reads were loaded. Ensure --pairedInput and/or --unpairedInput are set and valid.");
        }
        return groupedReads;
    }

    /**
     * Helper for groupReadsIntoPairs()
     */
    private static Iterator<Iterable<GATKRead>> groupPairedReadsPartition(final Iterator<GATKRead> iter,
                                                                          final int readsPerPartitionGuess) {
        //Traverse name-sorted partition, pairing reads as we go
        final ArrayList<Iterable<GATKRead>> newPartitionList = new ArrayList<>(readsPerPartitionGuess / 2);
        while (iter.hasNext()) {
            final GATKRead read1 = iter.next();
            GATKRead read2 = null;
            if (iter.hasNext()) {
                read2 = iter.next();
            }
            if (read2 == null || !read1.getName().equals(read2.getName())) {
                throw new UserException.BadInput("Found an unpaired read but expected all reads to be paired: " + read1.getName());
            }
            final List<GATKRead> pair = new ArrayList<>(2);
            pair.add(read1);
            pair.add(read2);
            newPartitionList.add(pair);
        }

        //Minimize memory footprint (don't rely on readsPerPartitionGuess)
        newPartitionList.trimToSize();

        return newPartitionList.iterator();
    }

    /**
     * Writes accessions contained in a SAM header that do not exist in the taxonomy database
     */
    public static void writeMissingReferenceAccessions(final String path, final SAMFileHeader header, final PSTaxonomyDatabase taxDB,
                                                       final Logger logger) {
        if (header != null && header.getSequenceDictionary() != null && header.getSequenceDictionary().getSequences() != null) {
            final Set<String> unknownSequences = header.getSequenceDictionary().getSequences().stream()
                    .map(SAMSequenceRecord::getSequenceName)
                    .filter(name -> !taxDB.accessionToTaxId.containsKey(name))
                    .collect(Collectors.toSet());
            try (final PrintStream file = new PrintStream(BucketUtils.createFile(path))) {
                unknownSequences.stream().forEach(file::print);
                if (file.checkError()) {
                    logger.warn("Error writing to header warnings file");
                }
            }
        }

    }

    /**
     * Gets taxonomic IDs of contigs that aligned sufficiently well to the reads. If both mates are present, returns the
     * intersection of their hits. Also sets read tag HITS_TAG to a comma-separated list of the IDs.
     */
    static JavaRDD<Tuple2<Iterable<GATKRead>, PSPathogenAlignmentHit>> mapGroupedReadsToTax(final JavaRDD<Iterable<GATKRead>> pairs,
                                                                                            final double minCoverage,
                                                                                            final double minIdentity,
                                                                                            final Broadcast<PSTaxonomyDatabase> taxonomyDatabaseBroadcast) {
        return pairs.map(readIter -> {

            //Number of reads in the pair (1 for unpaired reads)
            final int numReads = (int) Utils.stream(readIter).count();

            //Get tax IDs of all alignments in all reads that meet the coverage/identity criteria.
            final Stream<String> taxIds = Utils.stream(readIter)
                    .flatMap(read -> getValidHits(read, taxonomyDatabaseBroadcast.value(), minCoverage, minIdentity).stream());

            //Get list of tax IDs that are hits in all reads
            final List<String> hitTaxIds;
            if (numReads > 1) {

                //Group the flattened stream by tax id, e.g. 3453 -> {3453, 3453}, 938 -> {938}, etc., so that the
                // length of the list is the number of reads with that tax ID. Then map the lists to list lengths.
                final Map<String, Long> taxIdCounts = taxIds.collect(Collectors.groupingBy(e -> e, Collectors.counting()));

                //Filter hits that didn't occur in all reads
                hitTaxIds = taxIdCounts.entrySet().stream().map(entry -> entry.getValue() == numReads ? entry.getKey() : null)
                        .filter(Objects::nonNull).collect(Collectors.toList());

            } else {
                //Unpaired reads
                hitTaxIds = taxIds.collect(Collectors.toList());
            }

            final PSPathogenAlignmentHit info = new PSPathogenAlignmentHit(hitTaxIds, numReads);

            //If there was at least one hit, append a tag to each read with the list of hits
            if (hitTaxIds.size() > 0) {
                final String hitString = String.join(",", hitTaxIds);
                Utils.stream(readIter).forEach(read -> read.setAttribute(HITS_TAG, hitString));
            }
            return new Tuple2<>(readIter, info);
        });
    }


    /**
     * Gets set of sufficiently well-mapped hits
     */
    private static Set<String> getValidHits(final GATKRead read,
                                            PSTaxonomyDatabase taxonomyDatabase,
                                            final double minCoverage,
                                            final double minIdentity) {

        //Short circuit if read isn't mapped
        if (read.isUnmapped()) return Collections.emptySet();

        //Check that NM tag is set
        if (!read.hasAttribute("NM")) {
            throw new UserException.BadInput("SAM flag indicates a read is mapped, but the NM tag is absent");
        }

        //Find and return all alignments meeting the coverage/identity criteria
        final Collection<String> hits = new ArrayList<>();
        if (isValidAlignment(read.getCigar(), read.getAttributeAsInteger("NM"), minCoverage, minIdentity)) {
            hits.add(SAMSequenceRecord.truncateSequenceName(read.getAssignedContig()));
        }
        hits.addAll(getValidAlternateHits(read, "XA", 0, 2, 3, minCoverage, minIdentity));
        hits.addAll(getValidAlternateHits(read, "SA", 0, 3, 5, minCoverage, minIdentity));

        //Throw out duplicates and accessions not in the taxonomic database so it returns a list of unique tax ID's
        // for each read in the pair
        return hits.stream().map(contig -> taxonomyDatabase.accessionToTaxId.containsKey(contig) ? taxonomyDatabase.accessionToTaxId.get(contig) : null)
                .filter(Objects::nonNull).collect(Collectors.toSet());
    }

    /**
     * Parses alternate alignments tag of form TAG:Z:val_A0,val_A1,val_A2,...;val_B0,val_B1,val_B2,...;...
     * Each alignment is delimited by a semicolon (";") and is represented as a sub-list of values containing the
     * accession, cigar, and number of mismatches, which are delimited by commas (",") and expected at the given indices.
     */
    private static Collection<String> getValidAlternateHits(final GATKRead read, final String tag, final int contigIndex,
                                                            final int cigarIndex, final int numMismatchesIndex, final double minCoverage,
                                                            final double minIdentity) {
        final Collection<String> alternateHits = new ArrayList<>();
        if (read.hasAttribute(tag)) {
            final int expectedTokens = Math.max(contigIndex, Math.max(cigarIndex, numMismatchesIndex)) + 1;
            final String tagValue = read.getAttributeAsString(tag);
            final String[] tagTokens = tagValue.split(";");
            for (final String tok : tagTokens) {
                final String[] subtokens = tok.split(",");
                if (subtokens.length < expectedTokens) {
                    throw new UserException.BadInput("Error parsing " + tag + " tag: expected at least " + expectedTokens + " values per alignment but found " + subtokens.length);
                }
                final String recordName = SAMSequenceRecord.truncateSequenceName(subtokens[contigIndex]);
                final Cigar cigar = TextCigarCodec.decode(subtokens[cigarIndex]);
                final int numMismatches = Integer.valueOf(subtokens[numMismatchesIndex]);
                if (isValidAlignment(cigar, numMismatches, minCoverage, minIdentity)) {
                    alternateHits.add(recordName);
                }
            }
        }
        return alternateHits;
    }


    /**
     * Returns true if the candidate alignment of the read meets coverage and identity criteria
     */
    public static boolean isValidAlignment(final Cigar cigar, final int numMismatches, final double minCoverage, final double minIdentity) {
        final int minCoveredBases = (int) Math.ceil(minCoverage * cigar.getReadLength());
        int numCoveredBases = 0;
        final List<CigarElement> cigarElements = cigar.getCigarElements();
        for (final CigarElement cigarElement : cigarElements) {
            if (cigarElement.getOperator().equals(CigarOperator.MATCH_OR_MISMATCH) || cigarElement.getOperator().equals(CigarOperator.INSERTION)) {
                numCoveredBases += cigarElement.getLength();
            }
        }
        final int minIdentityBases = (int) Math.ceil(minIdentity * numCoveredBases);
        return !(new HostAlignmentReadFilter(minCoveredBases, minIdentityBases)).test(cigar, numMismatches);
    }

    /**
     * Computes abundance scores and returns key-values of taxonomic id and scores
     */
    public static Iterator<Tuple2<String, PSPathogenTaxonScore>> computeTaxScores(final Iterator<PSPathogenAlignmentHit> taxonHits,
                                                                                  final PSTaxonomyDatabase taxonomyDatabase) {
        final PSTree tree = taxonomyDatabase.tree;
        final Map<String, PSPathogenTaxonScore> taxIdsToScores = new HashMap<>();
        final Set<String> invalidIds = new HashSet<>();
        while (taxonHits.hasNext()) {
            final PSPathogenAlignmentHit hit = taxonHits.next();
            final Collection<String> hitTaxIds = new ArrayList<>(hit.taxIDs);
            final Set<String> hitInvalidTaxIds = new HashSet<>(SVUtils.hashMapCapacity(hitTaxIds.size()));
            for (final String taxId : hitTaxIds) {
                if (!tree.hasNode(taxId) || tree.getLengthOf(taxId) == 0) hitInvalidTaxIds.add(taxId);
            }
            hitTaxIds.removeAll(hitInvalidTaxIds);
            invalidIds.addAll(hitInvalidTaxIds);

            //Number of genomes hit by this read and number of mates in the tuple (1 for single, 2 for pair)
            final int numHits = hitTaxIds.size();
            if (numHits == 0) continue;

            //Unambiguous read scores for the lowest common ancestor and its ancestors
            final String lowestCommonAncestor = tree.getLCA(hitTaxIds);
            final List<String> lcaPath = tree.getPathOf(lowestCommonAncestor);
            for (final String taxId : lcaPath) {
                getOrAddScoreInfo(taxId, taxIdsToScores, tree).unambiguousReads += hit.numMates;
            }

            //Scores normalized by genome length and degree of ambiguity (number of hits)
            final Set<String> hitPathNodes = new HashSet<>(); //Set of all unique hits and ancestors
            for (final String taxId : hitTaxIds) {
                if (!tree.hasNode(taxId) || tree.getLengthOf(taxId) == 0) {
                    invalidIds.add(taxId);
                    continue;
                }
                final Double score = SCORE_GENOME_LENGTH_UNITS * hit.numMates / (numHits * tree.getLengthOf(taxId));
                //Git list containing this node and its ancestors
                final List<String> path = tree.getPathOf(taxId);
                hitPathNodes.addAll(path);
                for (final String pathTaxId : path) {
                    final PSPathogenTaxonScore info = getOrAddScoreInfo(pathTaxId, taxIdsToScores, tree);
                    if (pathTaxId.equals(taxId)) {
                        info.selfScore += score;
                    } else {
                        info.descendentScore += score;
                    }
                    taxIdsToScores.put(pathTaxId, info);
                }
            }

            //"reads" score is the number of reads that COULD belong to each node i.e. an upper-bound
            for (final String taxId : hitPathNodes) {
                getOrAddScoreInfo(taxId, taxIdsToScores, tree).totalReads += hit.numMates;
            }
        }
        PSUtils.logItemizedWarning(logger, invalidIds, "The following taxonomic ID hits were ignored because " +
                "they either could not be found in the tree or had a reference length of 0 (this may happen when " +
                "the catalog file, taxdump file, and/or pathogen reference are inconsistent)");
        return taxIdsToScores.entrySet().stream().map(entry -> new Tuple2<>(entry.getKey(), entry.getValue())).iterator();
    }

    /**
     * Assigns scores normalized to 100%. For each taxon, its normalized score is own score divided by the sum
     * over all scores, plus the sum of its childrens' normalized scores.
     */
    final static Map<String, PSPathogenTaxonScore> computeNormalizedScores(final Map<String, PSPathogenTaxonScore> taxIdsToScores,
                                                                           final PSTree tree) {
        //Get sum of all scores that were assigned directly to each taxa (as opposed to being propagated up from descendents)
        double sum = 0.;
        for (final PSPathogenTaxonScore score : taxIdsToScores.values()) {
            sum += score.selfScore;
        }

        //Gets normalized selfScores and adds it to all ancestors
        for (final Map.Entry<String, PSPathogenTaxonScore> entry : taxIdsToScores.entrySet()) {
            final String taxId = entry.getKey();
            final double selfScore = entry.getValue().selfScore;
            final double normalizedScore;
            if (sum == 0) {
                normalizedScore = 0;
            } else {
                normalizedScore = 100.0 * selfScore / sum;
            }
            final List<String> path = tree.getPathOf(taxId);
            for (final String pathTaxId : path) {
                taxIdsToScores.get(pathTaxId).scoreNormalized += normalizedScore;
            }
        }
        return taxIdsToScores;
    }

    /**
     * Helper function for handling PSPathogenTaxonScore retrieval from the taxScores map
     */
    private static PSPathogenTaxonScore getOrAddScoreInfo(final String taxIds,
                                                          final Map<String, PSPathogenTaxonScore> taxScores,
                                                          final PSTree tree) {
        final PSPathogenTaxonScore score;
        if (taxScores.containsKey(taxIds)) {
            score = taxScores.get(taxIds);
        } else {
            score = new PSPathogenTaxonScore();
            score.referenceLength = tree.getLengthOf(taxIds);
            taxScores.put(taxIds, score);
        }
        return score;
    }

    /**
     * Reads taxonomy database that has been serialized to a file
     */
    @SuppressWarnings("unchecked")
    public static PSTaxonomyDatabase readTaxonomyDatabase(final String filePath) {
        final Kryo kryo = new Kryo();
        kryo.setReferences(false);
        final Input input = new Input(BucketUtils.openFile(filePath));
        final PSTaxonomyDatabase taxonomyDatabase = kryo.readObject(input, PSTaxonomyDatabase.class);
        input.close();
        return taxonomyDatabase;
    }

    /**
     * Output a tab-delimited table of taxonomic scores
     */
    public static void writeScoresFile(final Map<String, PSPathogenTaxonScore> scores,
                                       final PSTree tree, final String filePath) {
        final String header = "tax_id\ttaxonomy\ttype\tname\t" + PSPathogenTaxonScore.outputHeader;
        try (final PrintStream printStream = new PrintStream(BucketUtils.createFile(filePath))) {
            printStream.println(header);
            for (final String key : scores.keySet()) {
                final String name = tree.getNameOf(key);
                final String rank = tree.getRankOf(key);
                final List<String> path = tree.getPathOf(key).stream().map(tree::getNameOf).collect(Collectors.toList());
                Collections.reverse(path);
                final String taxonomy = String.join("|", path);
                final String line = key + "\t" + taxonomy + "\t" + rank + "\t" + name + "\t" + scores.get(key);
                printStream.println(line.replace(" ", "_"));
            }
            if (printStream.checkError()) {
                throw new UserException.CouldNotCreateOutputFile(filePath, new IOException());
            }
        }
    }

}
