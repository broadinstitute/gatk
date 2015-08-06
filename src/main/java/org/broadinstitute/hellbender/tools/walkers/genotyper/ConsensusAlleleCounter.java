package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.AlignmentContextUtils;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Code for determining which indels are segregating among the samples.
 */
public class ConsensusAlleleCounter {
    final protected static Logger logger = Logger.getLogger(ConsensusAlleleCounter.class);
    private final int minIndelCountForGenotyping;
    private final boolean doMultiAllelicCalls;
    private final double minFractionInOneSample;

    public ConsensusAlleleCounter(final boolean doMultiAllelicCalls,
                                  final int minIndelCountForGenotyping,
                                  final double minFractionInOneSample) {
        this.minIndelCountForGenotyping = minIndelCountForGenotyping;
        this.doMultiAllelicCalls = doMultiAllelicCalls;
        this.minFractionInOneSample = minFractionInOneSample;
    }

    /**
     * Returns a list of Alleles at this locus that may be segregating
     * 
     * @param ref
     * @param contexts
     * @param contextType
     * @return
     */
    public List<Allele> computeConsensusAlleles(final ReferenceContext ref,
                                                final Map<String, AlignmentContext> contexts,
                                                final AlignmentContextUtils.ReadOrientation contextType) {
        final Map<String, Integer> consensusIndelStrings = countConsensusAlleles(ref, contexts, contextType);
        return consensusCountsToAlleles(ref, consensusIndelStrings);
    }

    private Map<String, Integer> countConsensusAlleles(final ReferenceContext ref,
                                                       final Map<String, AlignmentContext> contexts,
                                                       final AlignmentContextUtils.ReadOrientation contextType) {
        final Locatable loc = ref.getInterval();
        final HashMap<String, Integer> consensusIndelStrings = new HashMap<>();

        int insCount = 0, delCount = 0;
        // quick check of total number of indels in pileup
        for ( final Map.Entry<String, AlignmentContext> sample : contexts.entrySet() ) {
            final AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadPileup indelPileup = context.getBasePileup();
            insCount += indelPileup.getNumberOfElements(p->p.isBeforeInsertion());
            delCount += indelPileup.getNumberOfElements(p -> p.isBeforeDeletionStart());
        }

        if ( insCount < minIndelCountForGenotyping && delCount < minIndelCountForGenotyping ) {
            return Collections.emptyMap();
        }

        for (final Map.Entry<String, AlignmentContext> sample : contexts.entrySet()) {
            // todo -- warning, can be duplicating expensive partition here
            final AlignmentContext context = AlignmentContextUtils.stratify(sample.getValue(), contextType);

            final ReadPileup indelPileup = context.getBasePileup();

            final int nIndelReads = indelPileup.getNumberOfElements(p -> p.isBeforeInsertion()) + indelPileup.getNumberOfElements(p->p.isBeforeDeletionStart());
            final int nReadsOverall = indelPileup.size();

            if ( nIndelReads == 0 || (nIndelReads / (1.0 * nReadsOverall)) < minFractionInOneSample ) {
                continue;
            }

            for (final PileupElement p : indelPileup) {
                final GATKRead read = ReadClipper.hardClipAdaptorSequence(p.getRead());
                if (read == null) {
                    continue;
                }

                if ( p.isBeforeInsertion() ) {
                    final String insertionBases = p.getBasesOfImmediatelyFollowingInsertion();
                    // edge case: ignore a deletion immediately preceding an insertion as p.getBasesOfImmediatelyFollowingInsertion() returns null [EB]
                    if ( insertionBases == null ) {
                        continue;
                    }

                    boolean foundKey = false;
                    // copy of hashmap into temp arrayList
                    final ArrayList<Pair<String,Integer>> cList = new ArrayList<>();
                    for (final Map.Entry<String, Integer> s : consensusIndelStrings.entrySet()) {
                        cList.add(Pair.of(s.getKey(), s.getValue()));
                    }

                    if (read.getEnd() == loc.getStart()) {
                        // first corner condition: a read has an insertion at the end, and we're right at the insertion.
                        // In this case, the read could have any of the inserted bases and we need to build a consensus

                        for (int k=0; k < cList.size(); k++) {
                            final String s = cList.get(k).getLeft();
                            final int cnt = cList.get(k).getRight();
                            // case 1: current insertion is prefix of indel in hash map
                            if (s.startsWith(insertionBases)) {
                                cList.set(k, Pair.of(s,cnt+1));
                                foundKey = true;
                            }
                            else if (insertionBases.startsWith(s)) {
                                // case 2: indel stored in hash table is prefix of current insertion
                                // In this case, new bases are new key.
                                foundKey = true;
                                cList.set(k, Pair.of(insertionBases, cnt + 1));
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                        {
                            cList.add(Pair.of(insertionBases, 1));
                        }

                    }
                    else if (read.getStart() == loc.getStart()+1) {
                        // opposite corner condition: read will start at current locus with an insertion
                        for (int k=0; k < cList.size(); k++) {
                            final String s = cList.get(k).getLeft();
                            final int cnt = cList.get(k).getRight();
                            if (s.endsWith(insertionBases)) {
                                // case 1: current insertion (indelString) is suffix of indel in hash map (s)
                                cList.set(k,Pair.of(s,cnt+1));
                                foundKey = true;
                            }
                            else if (insertionBases.endsWith(s)) {
                                // case 2: indel stored in hash table is prefix of current insertion
                                // In this case, new bases are new key.
                                foundKey = true;
                                cList.set(k,Pair.of(insertionBases,cnt+1));
                            }
                        }
                        if (!foundKey)
                            // none of the above: event bases not supported by previous table, so add new key
                        {
                            cList.add(Pair.of(insertionBases, 1));
                        }


                    }
                    else {
                        // normal case: insertion somewhere in the middle of a read: add count to arrayList
                        final int cnt = consensusIndelStrings.containsKey(insertionBases)? consensusIndelStrings.get(insertionBases):0;
                        cList.add(Pair.of(insertionBases,cnt+1));
                    }

                    // copy back arrayList into hashMap
                    consensusIndelStrings.clear();
                    for (final Pair<String,Integer> pair : cList) {
                        consensusIndelStrings.put(pair.getLeft(),pair.getRight());
                    }

                }
                else if ( p.isBeforeDeletionStart() ) {
                    final String deletionString = String.format("D%d", p.getLengthOfImmediatelyFollowingIndel());
                    final int cnt = consensusIndelStrings.containsKey(deletionString)? consensusIndelStrings.get(deletionString):0;
                    consensusIndelStrings.put(deletionString,cnt+1);
                }
            }
        }

        return consensusIndelStrings;
    }

    private List<Allele> consensusCountsToAlleles(final ReferenceContext ref,
                                                  final Map<String, Integer> consensusIndelStrings) {
        final Locatable loc = ref.getInterval();
        final Collection<VariantContext> vcs = new ArrayList<>();
        int maxAlleleCnt = 0;
        Allele refAllele, altAllele;

        for (final Map.Entry<String, Integer> elt : consensusIndelStrings.entrySet()) {
            final String s = elt.getKey();
            final int curCnt = elt.getValue();
            int stop = 0;

            // if observed count if above minimum threshold, we will genotype this allele
            if (curCnt < minIndelCountForGenotyping) {
                continue;
            }

            if (s.startsWith("D")) {
                // get deletion length
                final int dLen = Integer.valueOf(s.substring(1));
                // get ref bases of accurate deletion
                final int startIdxInReference = 1 + loc.getStart() - ref.getWindow().getStart();
                stop = loc.getStart() + dLen;
                final byte[] refBases = Arrays.copyOfRange(ref.getBases(), startIdxInReference - 1, startIdxInReference + dLen);   // add reference padding

                if (Allele.acceptableAlleleBases(refBases, false)) {
                    refAllele = Allele.create(refBases, true);
                    altAllele = Allele.create(ref.getBase(), false);
                }
                else {
                    continue; // don't go on with this allele if refBases are non-standard
                }
            } else {
                // insertion case
                final String insertionBases = (char)ref.getBase() + s;  // add reference padding
                if (Allele.acceptableAlleleBases(insertionBases, false)) { // don't allow N's in insertions
                    refAllele = Allele.create(ref.getBase(), true);
                    altAllele = Allele.create(insertionBases, false);
                    stop = loc.getStart();
                }
                else {
                    continue; // go on to next allele if consensus insertion has any non-standard base.
                }
            }


            final VariantContextBuilder builder = new VariantContextBuilder().source("");
            builder.loc(loc.getContig(), loc.getStart(), stop);
            builder.alleles(Arrays.asList(refAllele, altAllele));
            builder.noGenotypes();
            if (doMultiAllelicCalls) {
                vcs.add(builder.make());
                if (vcs.size() >= GenotypeLikelihoods.MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED) {
                    break;
                }
            } else if (curCnt > maxAlleleCnt) {
                maxAlleleCnt = curCnt;
                vcs.clear();
                vcs.add(builder.make());
            }
        }

        if (vcs.isEmpty()) {
            return Collections.emptyList(); // nothing else to do, no alleles passed minimum count criterion
        }

        final VariantContext mergedVC = GATKVariantContextUtils.simpleMerge(vcs, null, GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, GATKVariantContextUtils.GenotypeMergeType.UNSORTED, false, false, null, false, false);
        return mergedVC.getAlleles();
    }
}
