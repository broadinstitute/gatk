package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import htsjdk.samtools.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.Charset;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public final class PSScoreUtils {

    public static final String HITS_TAG = "YP";
    public static final double SCORE_NORMALIZATION_FACTOR = 1e6;
    private static final Logger logger = LogManager.getLogger(PSScoreUtils.class);

    public static JavaRDD<GATKRead> scoreReads(final JavaSparkContext ctx,
                                               final JavaRDD<GATKRead> pairedReads,
                                               final JavaRDD<GATKRead> unpairedReads,
                                               final SAMFileHeader header,
                                               final PSScoreArgumentCollection scoreArgs) {

        //Group reads into pairs
        final JavaRDD<Iterable<GATKRead>> groupedReads = PSScoreUtils.groupReadsIntoPairs(pairedReads,
                unpairedReads, scoreArgs.readsPerPartition);

        //Load taxonomy database, created by running PathSeqBuildReferenceTaxonomy with this reference
        final PSTaxonomyDatabase taxDB = PSScoreUtils.readTaxonomyDatabase(scoreArgs.taxonomyDatabasePath);
        final Broadcast<Map<String, String>> accessionToTaxIdBroadcast = ctx.broadcast(taxDB.accessionToTaxId);

        //Check header against database
        if (scoreArgs.headerWarningFile != null) {
            PSScoreUtils.writeHeaderWarnings(scoreArgs.headerWarningFile, header, taxDB, logger);
        }

        //Determine which alignments are valid hits and return their tax IDs in PSPathogenAlignmentHit
        //Also adds pathseq tags containing the hit IDs to the reads
        final JavaRDD<Tuple2<Iterable<GATKRead>, PSPathogenAlignmentHit>> readHits = PSScoreUtils.mapGroupedReadsToTax(groupedReads,
                scoreArgs.minCoverage, scoreArgs.minIdentity, accessionToTaxIdBroadcast);

        //Collect PSPathogenAlignmentHit objects
        final Iterable<PSPathogenAlignmentHit> hitInfo = PSScoreUtils.collectTupleSecond(readHits);

        //Get the original reads, now with their pathseq hit tags set
        final JavaRDD<GATKRead> readsFinal = PSScoreUtils.flattenTupleFirstIterable(readHits);

        //Compute taxonomic scores
        final Map<String, PSPathogenTaxonScore> taxScores = PSScoreUtils.computeTaxScores(hitInfo, taxDB.tree);

        //Write scores to file
        PSScoreUtils.writeScoresFile(taxScores, taxDB.tree, scoreArgs.scoresPath);

        return readsFinal;
    }

    /**
     * Collects the second elements of all tuples in an RDD
     */
    public static <T, C> Iterable<C> collectTupleSecond(final JavaRDD<Tuple2<T, C>> tupleRdd) {
        return tupleRdd.map(tuple -> tuple._2).collect();
    }

    /**
     * Flattens RDD of tuples, in which the first elements are Iterable, to an RDD of the items in the Iterables
     */
    public static <T, C> JavaRDD<T> flattenTupleFirstIterable(final JavaRDD<Tuple2<Iterable<T>, C>> tupleRdd) {
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
        final List<Iterable<GATKRead>> newPartitionList = new ArrayList<>(readsPerPartitionGuess / 2);
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
        final List<Iterable<GATKRead>> newPartitionListResized = new ArrayList<>(newPartitionList.size());
        newPartitionListResized.addAll(newPartitionList);

        return newPartitionListResized.iterator();
    }

    /**
     * Combines SAM file header sequences and read groups. If both headers are not null, additional header entries
     * from the unpaired header will not be copied over.
     */
    static SAMFileHeader joinBamHeaders(final SAMFileHeader pairedHeader, final SAMFileHeader unpairedHeader) {
        SAMFileHeader header;
        if (pairedHeader != null) {
            header = pairedHeader;
            if (unpairedHeader != null && !header.equals(unpairedHeader)) {
                //Add sequences
                for (final SAMSequenceRecord rec : unpairedHeader.getSequenceDictionary().getSequences()) {
                    if (header.getSequenceDictionary().getSequence(rec.getSequenceName()) == null) {
                        header.addSequence(rec);
                    }
                }
                //Add read groups
                for (final SAMReadGroupRecord rec : unpairedHeader.getReadGroups()) {
                    if (header.getReadGroup(rec.getReadGroupId()) == null) {
                        header.addReadGroup(rec);
                    }
                }
            }
        } else if (unpairedHeader != null) {
            header = unpairedHeader;
        } else {
            throw new UserException.BadInput("No headers were loaded");
        }
        return header;
    }

    /**
     * Writes accessions contained in a SAM header that do not exist in the taxonomy database
     */
    public static void writeHeaderWarnings(final String path, final SAMFileHeader header, final PSTaxonomyDatabase taxDB,
                                           final Logger logger) {
        if (header != null && header.getSequenceDictionary() != null && header.getSequenceDictionary().getSequences() != null) {
            final Set<String> unknownSequences = header.getSequenceDictionary().getSequences().stream()
                    .map(SAMSequenceRecord::getSequenceName)
                    .filter(name -> !taxDB.accessionToTaxId.containsKey(name))
                    .collect(Collectors.toSet());
            try {
                final OutputStream file = BucketUtils.createFile(path);
                for (final String acc : unknownSequences) {
                    final String line = acc + "\n";
                    file.write(line.getBytes(Charset.defaultCharset()));
                }
                file.close();
            } catch (final IOException e) {
                logger.warn("Could not write header warnings to " + path, e);
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
                                                                                            final Broadcast<Map<String, String>> refNameToTaxBroadcast) {
        return pairs.map(readIter -> {

            final Collection<List<String>> taxIDLists = Utils.stream(readIter)
                    //Flat map the taxID's of all alignments meeting the coverage/identity criteria.
                    //Throw out duplicates so it returns a list of unique tax ID's for each read in the pair
                    .flatMap(read -> getValidHits(read, minCoverage, minIdentity).stream()
                            .map(contig -> refNameToTaxBroadcast.value().containsKey(contig) ? refNameToTaxBroadcast.value().get(contig) : null)
                            .filter(Objects::nonNull).distinct()
                    ).collect(Collectors.groupingBy(Function.identity())).values(); //Group the flattened stream by tax ID into lists

            final int numReads = (int) Utils.stream(readIter).count();

            //Count the length of each list - if it's equal to the number of reads, we have a hit
            final List<String> hitTaxIds = taxIDLists.stream()
                    .filter(ids -> ids.size() == numReads)
                    .map(ids -> ids.iterator().next())
                    .collect(Collectors.toList());

            final PSPathogenAlignmentHit info = new PSPathogenAlignmentHit(hitTaxIds, numReads);

            //If there was at least one hit, append a tag to each read with the list of hits
            if (hitTaxIds.size() > 0) {
                final StringBuilder stringBuilder = new StringBuilder();
                for (final String hit : hitTaxIds) {
                    stringBuilder.append(hit);
                    stringBuilder.append(",");
                }
                final String hitString = stringBuilder.substring(0, stringBuilder.length() - 1);
                Utils.stream(readIter).forEach(read -> read.setAttribute(HITS_TAG, hitString));
            } else {
                Utils.stream(readIter).forEach(GATKRead::setIsUnmapped);
            }
            return new Tuple2<>(readIter, info);
        });
    }


    /**
     * Gets set of sufficiently well-mapped hits
     */
    private static Collection<String> getValidHits(final GATKRead read, final double minCoverage, final double minIdentity) {

        //Short circuit if read isn't mapped
        if (read.isUnmapped()) return new ArrayList<>(0);

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
        return hits;
    }

    /**
     * Parses alternate alignments tag of form TAG:Z:val_A0,val_A1,val_A2,...;val_B0,val_B1,val_B2,...;...
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
     * Computes abundance scores and returns map from taxonomic id to scores
     */
    public static Map<String, PSPathogenTaxonScore> computeTaxScores(final Iterable<PSPathogenAlignmentHit> readTaxHits,
                                                                     final PSTree tree) {
        final Map<String, PSPathogenTaxonScore> taxIdsToScores = new HashMap<>();
        final Set<String> invalidIds = new HashSet<>();
        Double sum = 0.0;
        for (final PSPathogenAlignmentHit hit : readTaxHits) {
            final Collection<String> hitTaxIds = new ArrayList<>(hit.taxIDs);

            //Find and omit hits to tax ID's not in the database
            final Set<String> invalidHitIds = new HashSet<>(hitTaxIds.size());
            for (final String taxid : hitTaxIds) {
                if (!tree.hasNode(taxid) || tree.getLengthOf(taxid) == 0) {
                    invalidHitIds.add(taxid);
                }
            }
            hitTaxIds.removeAll(invalidHitIds);
            invalidIds.addAll(invalidHitIds);

            //Number of genomes hit by this read and number of mates in the tuple (1 for single, 2 for pair)
            final int numHits = hitTaxIds.size();
            if (numHits == 0) continue;

            //Unambiguous read scores for the lowest common ancestor and its ancestors
            final String lowestCommonAncestor = tree.getLCA(hitTaxIds);
            final List<String> lcaPath = tree.getPathOf(lowestCommonAncestor);
            for (final String taxId : lcaPath) {
                getOrAddScoreInfo(taxId, taxIdsToScores, tree).unambiguous += hit.numMates;
            }

            //Scores normalized by genome length and degree of ambiguity (number of hits)
            final Set<String> hitPathNodes = new HashSet<>(); //Set of all unique hits and ancestors
            for (final String taxId : hitTaxIds) {
                final Double score = SCORE_NORMALIZATION_FACTOR * hit.numMates / (numHits * tree.getLengthOf(taxId));
                sum += score;
                final List<String> path = tree.getPathOf(taxId);
                hitPathNodes.addAll(path);
                for (final String pathTaxID : path) {
                    PSPathogenTaxonScore info = getOrAddScoreInfo(pathTaxID, taxIdsToScores, tree);
                    info.score += score;
                    taxIdsToScores.put(pathTaxID, info);
                }
            }

            //"reads" score is the number of reads that COULD belong to each node i.e. an upper-bound
            for (final String taxId : hitPathNodes) {
                getOrAddScoreInfo(taxId, taxIdsToScores, tree).reads += hit.numMates;
            }
        }

        //Scores normalized to 100%
        for (final Map.Entry<String, PSPathogenTaxonScore> entry : taxIdsToScores.entrySet()) {
            final String readName = entry.getKey();
            final PSPathogenTaxonScore score = entry.getValue();
            if (sum == 0) {
                score.scoreNormalized = 0;
            } else {
                score.scoreNormalized = 100.0 * score.score / sum;
            }
            taxIdsToScores.replace(readName, score);
        }
        PSUtils.logItemizedWarning(logger, invalidIds, "The following taxonomic ID hits were ignored because " +
                "they either could not be found in the tree or had a reference length of 0 (this may happen when " +
                "the catalog file, taxdump file, and/or pathogen reference are inconsistent)");
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
            score.refLength = tree.getLengthOf(taxIds);
            taxScores.put(taxIds, score);
        }
        return score;
    }

    /**
     * Reads taxonomy database that has been written using PSUtils.writeKryoTwo()
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
        try {
            final OutputStream file = BucketUtils.createFile(filePath);
            final String header = "tax_id\ttaxonomy\ttype\tname\t" + PSPathogenTaxonScore.outputHeader + "\n";
            file.write(header.getBytes(Charset.defaultCharset()));
            for (final String key : scores.keySet()) {
                final String name = tree.getNameOf(key).replace(" ", "_");
                final String rank = tree.getRankOf(key).replace(" ", "_");
                final List<String> path = tree.getPathOf(key);
                StringBuilder stringBuilder = new StringBuilder();
                for (int i = path.size() - 1; i >= 0; i--) {
                    stringBuilder.append(tree.getNameOf(path.get(i)).replace(" ", "_"));
                    stringBuilder.append("|");
                }
                final String taxonomy = stringBuilder.toString();
                final String line = key + "\t" + taxonomy + "\t" + rank + "\t" + name + "\t" + scores.get(key) + "\n";
                file.write(line.getBytes(Charset.defaultCharset()));
            }
            file.close();
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(filePath, e);
        }
    }

}
