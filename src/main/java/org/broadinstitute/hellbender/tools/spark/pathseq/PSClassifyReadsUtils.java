package org.broadinstitute.hellbender.tools.spark.pathseq;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import htsjdk.samtools.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
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

public final class PSClassifyReadsUtils {

    public static final String HITS_TAG = "YP";
    public static final double SCORE_NORMALIZATION_FACTOR = 1e6;
    private static final Logger logger = LogManager.getLogger(PSClassifyReadsUtils.class);

    public static <T> Iterable<T> wrapInIterable(final T obj) {
        final List<T> list = new ArrayList<>(1);
        list.add(obj);
        return list;
    }

    public static <T, C> Iterable<C> collectTupleSecond(final JavaRDD<Tuple2<T, C>> tupleRDD) {
        return tupleRDD.map(tuple -> tuple._2).collect();
    }

    public static <T, C> JavaRDD<T> flattenTupleFirstIterable(final JavaRDD<Tuple2<Iterable<T>, C>> tupleRDD) {
        return tupleRDD.flatMap(tuple -> tuple._1.iterator());
    }

    public static void doHeaderWarnings(final String path, final SAMFileHeader header, final PSTaxonomyDatabase taxDB,
                                        final PipelineOptions options, final Logger logger) {
        if (header != null && header.getSequenceDictionary() != null && header.getSequenceDictionary().getSequences() != null) {
            final Set<String> unknownSequences = header.getSequenceDictionary().getSequences().stream()
                    .map(SAMSequenceRecord::getSequenceName)
                    .filter(name -> !taxDB.contigToTaxIDMap.containsKey(name))
                    .collect(Collectors.toSet());

            try {
                final OutputStream file = BucketUtils.createFile(path, options);
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
     * intersection of their hits.
     */
    public static JavaRDD<Tuple2<Iterable<GATKRead>, PSHitInfo>> mapGroupedReadsToTax(final JavaRDD<Iterable<GATKRead>> pairs,
                                                                                      final double min_cov,
                                                                                      final double min_ident,
                                                                                      final Broadcast<Map<String, String>> refNameToTaxBroadcast) {
        return pairs.map(readIter -> {

            final Collection<List<String>> taxIDLists = Utils.stream(readIter)
                    //Flat map the taxID's of all alignments meeting the coverage/identity criteria.
                    //Throw out duplicates so it returns a list of unique tax ID's for each read in the pair
                    .flatMap(read -> getValidHits(read, min_cov, min_ident).stream()
                            .map(contig -> refNameToTaxBroadcast.value().containsKey(contig) ? refNameToTaxBroadcast.value().get(contig) : null)
                            .filter(Objects::nonNull).distinct()
                    ).collect(Collectors.groupingBy(Function.identity())).values(); //Group the flattened stream by tax ID into lists

            final int numReads = (int) Utils.stream(readIter).count();

            //Count the length of each list - if it's equal to the number of reads, we have a hit
            final List<String> hits = taxIDLists.stream()
                    .filter(ids -> ids.size() == numReads)
                    .map(ids -> ids.iterator().next())
                    .collect(Collectors.toList());

            final PSHitInfo info = new PSHitInfo(hits, numReads);

            //If there was at least one hit, append a tag to each read with the list of hits
            if (hits.size() > 0) {
                String tmpStr = "";
                for (final String hit : hits) {
                    tmpStr += hit + ",";
                }
                final String hitStr = tmpStr.substring(0, tmpStr.length() - 1);
                Utils.stream(readIter).forEach(read -> read.setAttribute(HITS_TAG, hitStr));
            } else {
                Utils.stream(readIter).forEach(GATKRead::setIsUnmapped);
            }
            return new Tuple2<>(readIter, info);
        });
    }


    /**
     * Gets set of sufficiently well-mapped hits
     */
    private static Collection<String> getValidHits(final GATKRead read, final double min_cov, final double min_ident) {

        //Short circuit if read isn't mapped
        if (read.isUnmapped()) return new ArrayList<>(0);

        //Check that NM tag is set
        if (!read.hasAttribute("NM")) {
            throw new UserException.BadInput("SAM flag indicates a read is mapped, but the NM tag is absent");
        }

        //Find and return all alignments meeting the coverage/identity criteria
        final Collection<String> hits = new ArrayList<>();
        if (isValidAlignment(read.getCigar(), read.getAttributeAsInteger("NM"), min_cov, min_ident)) {
            hits.add(SAMSequenceRecord.truncateSequenceName(read.getAssignedContig()));
        }
        hits.addAll(getValidAlternateHits(read, "XA", 0, 2, 3, min_cov, min_ident));
        hits.addAll(getValidAlternateHits(read, "SA", 0, 3, 5, min_cov, min_ident));
        return hits;
    }

    /**
     * Parses alternate alignments tag of form TAG:Z:val_A0,val_A1,val_A2,...;val_B0,val_B1,val_B2,...;...
     */
    private static Collection<String> getValidAlternateHits(final GATKRead read, final String tag, final int contigIndex,
                                                            final int cigarIndex, final int NMIndex, final double min_cov,
                                                            final double min_ident) {
        final Collection<String> alternateHits = new ArrayList<>();
        if (read.hasAttribute(tag)) {
            final int expectedTokens = Math.max(contigIndex, Math.max(cigarIndex, NMIndex)) + 1;
            final String tagVal = read.getAttributeAsString(tag);
            final String[] tokens = tagVal.split(";");
            for (final String tok : tokens) {
                final String[] subtokens = tok.split(",");
                if (subtokens.length < expectedTokens) {
                    throw new UserException.BadInput("Error parsing " + tag + " tag: expected at least " + expectedTokens + " values per alignment but found " + subtokens.length);
                }
                final String contig = SAMSequenceRecord.truncateSequenceName(subtokens[contigIndex]);
                final Cigar c = TextCigarCodec.decode(subtokens[cigarIndex]);
                final int NM = Integer.valueOf(subtokens[NMIndex]);
                if (isValidAlignment(c, NM, min_cov, min_ident)) {
                    alternateHits.add(contig);
                }
            }
        }
        return alternateHits;
    }


    /**
     * Returns true if the candidate alignment of the read meets coverage and identity criteria
     */
    public static boolean isValidAlignment(final Cigar c, final int NM, final double min_cov, final double min_ident) {
        final int minCovBases = (int) Math.ceil(min_cov * c.getReadLength());
        int numCov = 0;
        final List<CigarElement> cigarElements = c.getCigarElements();
        for (final CigarElement e : cigarElements) {
            if (e.getOperator().equals(CigarOperator.MATCH_OR_MISMATCH) || e.getOperator().equals(CigarOperator.INSERTION)) {
                numCov += e.getLength();
            }
        }
        final int minIdentBases = (int) Math.ceil(min_ident * numCov);
        return !(new HostAlignmentReadFilter(minCovBases, minIdentBases)).test(c, NM);
    }

    /**
     * Computes abundance scores and returns map from taxonomic id to scores
     */
    public static Map<String, PSScoreInfo> computeTaxScores(final Iterable<PSHitInfo> readTaxHits,
                                                            final PSTree tree) {
        final Map<String, PSScoreInfo> taxScores = new HashMap<>();
        final Set<String> invalidIDs = new HashSet<>();
        Double sum = 0.0;
        for (final PSHitInfo hit : readTaxHits) {
            final Collection<String> taxIDs = new ArrayList<>(hit.taxIDs);

            //Find and omit hits to tax ID's not in the database
            final Set<String> invalidHitIDs = new HashSet<>(taxIDs.size());
            for (final String taxid : taxIDs) {
                if (!tree.hasNode(taxid) || tree.getLengthOf(taxid) == 0) {
                    invalidHitIDs.add(taxid);
                }
            }
            taxIDs.removeAll(invalidHitIDs);
            invalidIDs.addAll(invalidHitIDs);

            //Number of genomes hit by this read and number of mates in the tuple (1 for single, 2 for pair)
            final int numHits = taxIDs.size();
            if (numHits == 0) continue;

            //Unambiguous read scores for the lowest common ancestor and its ancestors
            final String lca = tree.getLCA(taxIDs);
            final List<String> lcaPath = tree.getPathOf(lca);
            for (final String taxid : lcaPath) {
                getOrAddScoreInfo(taxid, taxScores, tree).unambiguous += hit.numMates;
            }

            //Scores normalized by genome length and degree of ambiguity (number of hits)
            final Set<String> hitPathNodes = new HashSet<>(); //Set of all unique hits and ancestors
            for (final String taxid : taxIDs) {
                final Double score = SCORE_NORMALIZATION_FACTOR * hit.numMates / (numHits * tree.getLengthOf(taxid));
                sum += score;
                final List<String> path = tree.getPathOf(taxid);
                hitPathNodes.addAll(path);
                for (final String pathTaxID : path) {
                    PSScoreInfo info = getOrAddScoreInfo(pathTaxID, taxScores, tree);
                    info.score += score;
                    taxScores.put(pathTaxID, info);
                }
            }

            //"reads" score is the number of reads that COULD belong to each node i.e. an upper-bound
            for (final String taxid : hitPathNodes) {
                getOrAddScoreInfo(taxid, taxScores, tree).reads += hit.numMates;
            }
        }

        //Scores normalized to 100%
        for (final Map.Entry<String, PSScoreInfo> entry : taxScores.entrySet()) {
            final String readName = entry.getKey();
            final PSScoreInfo info = entry.getValue();
            if (sum == 0) {
                info.score_normalized = 0;
            } else {
                info.score_normalized = 100.0 * info.score / sum;
            }
            taxScores.replace(readName, info);
        }
        PSUtils.logItemizedWarning(logger, invalidIDs, "The following taxonomic ID hits were ignored because they either could not be found in the tree or had a reference length of 0 (this may happen when the catalog file, taxdump file, and/or pathogen reference are inconsistent)");
        return taxScores;
    }

    /**
     * Helper function for handling PSScoreInfo retrieval from the taxScores map
     */
    private static PSScoreInfo getOrAddScoreInfo(final String taxID, final Map<String, PSScoreInfo> taxScores, final PSTree tree) {
        PSScoreInfo info;
        if (taxScores.containsKey(taxID)) {
            info = taxScores.get(taxID);
        } else {
            info = new PSScoreInfo();
            info.ref_length = tree.getLengthOf(taxID);
            taxScores.put(taxID, info);
        }
        return info;
    }

    /**
     * Output a tab-delimited table of taxonomic scores
     */
    public static void writeScoresFile(final Map<String, PSScoreInfo> scores, final PSTree tree, final String filePath,
                                       final PipelineOptions options) {
        try {
            final OutputStream file = BucketUtils.createFile(filePath, options);
            final String header = "tax_id\ttaxonomy\ttype\tname\t" + PSScoreInfo.outputHeader + "\n";
            file.write(header.getBytes(Charset.defaultCharset()));
            for (final String key : scores.keySet()) {
                final String name = tree.getNameOf(key).replace(" ", "_");
                final String rank = tree.getRankOf(key).replace(" ", "_");
                final List<String> path = tree.getPathOf(key);
                String taxonomy = "";
                for (int i = path.size() - 1; i >= 0; i--) {
                    taxonomy += tree.getNameOf(path.get(i)).replace(" ", "_") + "|";
                }
                final String line = key + "\t" + taxonomy + "\t" + rank + "\t" + name + "\t" + scores.get(key) + "\n";
                file.write(line.getBytes(Charset.defaultCharset()));
            }
            file.close();
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(filePath, e);
        }
    }

}
