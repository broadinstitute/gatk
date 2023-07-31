package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.walkers.ErrorCorrectHiFi.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ByteSequence;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.*;
import java.util.function.Function;

@CommandLineProgramProperties(
        summary = "Find and phase small variants near an SV breakpoint.",
        oneLineSummary = "Find and phase small variants near an SV breakpoint.",
        programGroup = VariantEvaluationProgramGroup.class
)
public class DetermineSVHaplotypes extends VariantWalker {
    public static final int PADDING = 150;
    public static final int WINDOW_SIZE = 2 * PADDING;
    public static final int MIN_QUAL = 10;
    public static final int MIN_SMALL_INDEL_COUNT = 3;
    public static final float MIN_SMALL_INDEL_FRACT = 0.1F;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public void apply( final VariantContext variant,
                       final ReadsContext readsContext,
                       final ReferenceContext refContext,
                       final FeatureContext featureContext ) {
        final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, null);
        if ( svType == null ) {
            return;
        }
        final boolean isDel = GATKSVVCFConstants.SYMB_ALT_STRING_DEL.equals(svType);

        // The VCF references the base before the actual variant for deletions and the base after for duplications.
        // At least it's true in the VCF I'm testing with.
        final int variantStart = variant.getStart() +
                (isDel ? 1 : GATKSVVCFConstants.SYMB_ALT_STRING_DUP.equals(svType) ? -1 : 0);
        final int paddedStart = Math.max(1, variantStart - PADDING);
        final SimpleInterval readQueryWindow =
                new SimpleInterval(variant.getContig(), paddedStart, paddedStart);

        // make an iterator for the calls in each read that overlap the variant region
        final List<CallIterator> callIterators = new ArrayList<>();
        readsContext.iterator(readQueryWindow)
                .forEachRemaining(read -> callIterators.add(new CallIterator(read, paddedStart)));

        int paddedEnd = paddedStart + WINDOW_SIZE;
        final int nReads = callIterators.size();
        Map<Call, List<Integer>> callMap = new HashMap<>(nReads);
        final List<HetSite> hetSites = new ArrayList<>();
        final List<Call> consensusCalls = new ArrayList<>(WINDOW_SIZE);

        // find the heterozygous sites in a window around the variant site
        for ( int refLoc = paddedStart; refLoc < paddedEnd && fillCallMap(callMap, callIterators); ++refLoc ) {
            if ( isVariant(callMap, consensusCalls) ) {
                hetSites.add(new HetSite(refLoc, callMap));
                paddedEnd = adjustWindowBoundary(paddedEnd, callMap.keySet());
                // make a new collection, the old one has been stashed in the HetSite
                callMap = new HashMap<>(nReads);
            } else {
                paddedEnd = adjustWindowBoundary(paddedEnd, callMap.keySet());
                // reuse the old collection after clearing it
                callMap.clear();
            }
        }

        // phase the variants -- save the HetSite at the variantStart locus, if present
        final HetSite ourHetSite = phase(hetSites, variantStart);

        final List<ByteSequence> haplotypes = getHaplotypes(paddedStart, consensusCalls, hetSites);
        final CigarMatcher cigarMatcher = new CigarMatcher(variant);
        if ( ourHetSite == null ) {
            String genotype = "0/0";
            if ( consensusCalls.size() > PADDING ) {
                final Call ourCall = consensusCalls.get(PADDING);
                if ( ourCall != null && cigarMatcher.matches(ourCall.getCigarElement()) ) {
                    genotype = "1/1";
                }
            }
            writeHaplotypes(variant, haplotypes, genotype);
            return;
        }

        final boolean alt1Matches = cigarMatcher.matches(ourHetSite.getAlt1().getCigarElement());
        final boolean alt2Matches = cigarMatcher.matches(ourHetSite.getAlt2().getCigarElement());
        if ( alt1Matches ) {
            if ( alt2Matches ) {
                writeHaplotypes(variant, haplotypes, "1/1");
            } else {
                final ByteSequence alt1Sequence = haplotypes.get(0);
                haplotypes.set(0, haplotypes.get(1));
                haplotypes.set(1, alt1Sequence);
                writeHaplotypes(variant, haplotypes, "0/1");
            }
        } else if ( alt2Matches ) {
            writeHaplotypes(variant, haplotypes, "0/1");
        } else {
            writeHaplotypes(variant, haplotypes, "0/0");
        }
    }

    private static int adjustWindowBoundary( final int currentEnd, final Set<Call> calls ) {
        int adjustment = 0;
        for ( final Call call : calls ) {
            final CigarElement cigarElement = call.getCigarElement();
            if ( cigarElement.getOperator() == CigarOperator.D ) {
                final int len = cigarElement.getLength();
                if ( len > adjustment ) {
                    adjustment = len;
                }
            }
        }
        return currentEnd + adjustment;
    }

    private static boolean fillCallMap( final Map<Call, List<Integer>> callMap,
                                        final List<CallIterator> callIterators ) {
        boolean exhausted = true;
        final int nReads = callIterators.size();
        for ( int readIdx = 0; readIdx != nReads; ++readIdx ) {
            final CallIterator callIterator = callIterators.get(readIdx);
            final Call call;
            if ( callIterator.hasNext() ) {
                call = callIterator.next();
                exhausted = false;
            } else {
                call = null;
            }
            if ( call != null && call.getMeanQual() >= MIN_QUAL ) {
                final int callIdx = readIdx;
                callMap.compute(call, ( k, v ) -> {
                    if ( v == null ) v = new ArrayList<>();
                    v.add(callIdx);
                    return v;
                });
            }
        }
        return !exhausted;
    }

    private static boolean isVariant( final Map<Call, List<Integer>> callMap,
                               final List<Call> consensusCalls ) {
        int maxIDs = 0;
        final int totalIds = callMap.values().stream().mapToInt(List::size).sum();
        Map.Entry<Call, List<Integer>> maxEntry = null;
        final Iterator<Map.Entry<Call, List<Integer>>> itr = callMap.entrySet().iterator();
        while ( itr.hasNext() ) {
            final Map.Entry<Call, List<Integer>> entry = itr.next();
            final int nIDs = entry.getValue().size();
            if ( nIDs == 1 ) { // assume that variation seen in a single read is bogus
                itr.remove();
                continue;
            }
            // remove small indels that aren't well attested -- that's our main sequencing error mode
            final Call call = entry.getKey();
            final int keyLen = call.getCalls().length();
            // if we have a single-base deletion, or a single-base insertion
            if ( (keyLen == 0 && call.getCigarElement().getLength() == 1) || keyLen == 2 ) {
                if ( nIDs <= MIN_SMALL_INDEL_COUNT || 1.F*nIDs/totalIds < MIN_SMALL_INDEL_FRACT ) {
                    itr.remove();
                    continue;
                }
            }
            if ( nIDs > maxIDs ) {
                maxEntry = entry;
                maxIDs = nIDs;
            }
        }

        consensusCalls.add(maxEntry == null ? null : maxEntry.getKey());
        return callMap.size() > 1;
    }

    private static HetSite phase( final List<HetSite> hetSites, final int variantStart ) {
        final int nHetSites = hetSites.size();
        for ( int idx = 0; idx != nHetSites; ++idx ) {
            final HetSite hetSite = hetSites.get(idx);
            final int cmp = Integer.compare(hetSite.getRefLoc(), variantStart);
            if ( cmp == 0 || nHetSites == 1 ) {
                final ArrayList<Map.Entry<Call, List<Integer>>> locs =
                        new ArrayList<>(hetSite.getCallMap().entrySet());
                locs.sort(Comparator.comparingInt(e -> -e.getValue().size()));
                hetSite.setPhase(locs.get(0).getKey(), locs.get(1).getKey());
                if ( nHetSites > 1 ) {
                    phaseExtend(hetSites, idx, idx);
                }
                return cmp == 0 ? hetSite : null;
            }
            if ( cmp > 0 ) {
                break;
            }
        }
        for ( int idx1 = 0; idx1 < nHetSites - 1; ++idx1 ) {
            for ( int idx2 = idx1 + 1; idx2 < nHetSites; ++idx2 ) {
                if ( phasePair(hetSites.get(idx1), hetSites.get(idx2)) ) {
                    phaseExtend(hetSites, idx1, idx2);
                    return null;
                }
            }
        }

        // if we have only two variants, and they can't be phased together
        if ( nHetSites == 2 ) {
            // remove the one that's farther from the variant of interest
            if ( Math.abs(hetSites.get(0).getRefLoc() - variantStart) <
                    Math.abs(hetSites.get(1).getRefLoc() - variantStart) ) {
                hetSites.remove(1);
            } else {
                hetSites.remove(0);
            }
            // arbitrarily set the phase of the one that remains
            final HetSite hetSite = hetSites.get(0);
            final ArrayList<Map.Entry<Call, List<Integer>>> locs =
                    new ArrayList<>(hetSite.getCallMap().entrySet());
            locs.sort(Comparator.comparingInt(e -> -e.getValue().size()));
            hetSite.setPhase(locs.get(0).getKey(), locs.get(1).getKey());
            return null;
        }
        // we tried everything we could think of, but couldn't phase
        hetSites.clear();
        return null;
    }

    private static boolean phasePair( final HetSite loc1, final HetSite loc2 ) {
        final List<CallPair> pairs = new ArrayList<>();
        for ( Map.Entry<Call, List<Integer>> entry1 : loc1.getCallMap().entrySet() ) {
            final List<Integer> list1 = entry1.getValue();
            CallPair pairForThisEntry = null;
            for ( Map.Entry<Call, List<Integer>> entry2 : loc2.getCallMap().entrySet() ) {
                final List<Integer> list2 = entry2.getValue();
                final int minCommon = Math.round(0.75F * Math.max(list1.size(), list2.size()));
                if ( minCommon > Math.min(list1.size(), list2.size()) ) {
                    continue;
                }
                if ( countCommonIDs(list1, list2) >= minCommon ) {
                    if ( pairForThisEntry == null ) {
                        pairForThisEntry = new CallPair(entry1.getKey(), entry2.getKey());
                    } else {
                        pairForThisEntry = null;
                        break;
                    }
                }
            }
            if ( pairForThisEntry != null ) {
                pairs.add(pairForThisEntry);
            }
        }
        if ( pairs.size() != 2 ) {
            return false;
        }
        loc1.setPhase(pairs.get(0).call1, pairs.get(1).call1);
        loc2.setPhase(pairs.get(0).call2, pairs.get(1).call2);
        return true;
    }

    private static void phaseExtend( final List<HetSite> hetSites, final int idx1, int idx2 ) {
        // phase or delete variants downstream of idx2
        HetSite phasedLoc = hetSites.get(idx2);
        int idx3 = idx2 + 1;
        while ( idx3 < hetSites.size() ) {
            final HetSite unphasedLoc = hetSites.get(idx3);
            if ( extendPhase(phasedLoc, unphasedLoc) ) {
                phasedLoc = unphasedLoc;
                idx3 += 1;
            } else {
                hetSites.remove(idx3);
            }
        }

        // phase or delete variants between idx1 and idx2
        phasedLoc = hetSites.get(idx1);
        int idx4 = idx1 + 1;
        while ( idx4 < idx2 ) {
            final HetSite unphasedLoc = hetSites.get(idx4);
            if ( extendPhase(phasedLoc, unphasedLoc) ) {
                phasedLoc = unphasedLoc;
                idx4 += 1;
            } else {
                hetSites.remove(idx4);
                idx2 -= 1;
            }
        }

        // phase or delete variants upstream of idx1
        phasedLoc = hetSites.get(idx1);
        int idx0 = idx1 - 1;
        while ( idx0 >= 0 ) {
            final HetSite unphasedLoc = hetSites.get(idx0);
            if ( extendPhase(phasedLoc, unphasedLoc) ) {
                phasedLoc = unphasedLoc;
            } else {
                hetSites.remove(idx0);
            }
            idx0 -= 1;
        }
    }

    private static boolean extendPhase( final HetSite phasedLoc, final HetSite unphasedLoc ) {
        final Call call1 = findEntry(phasedLoc.getListFor(phasedLoc.getAlt1()), unphasedLoc);
        if ( call1 == null ) {
            return false;
        }
        final Call call2 = findEntry(phasedLoc.getListFor(phasedLoc.getAlt2()), unphasedLoc);
        if ( call2 == null || call1.equals(call2) ) {
            return false;
        }
        unphasedLoc.setPhase(call1, call2);
        return true;
    }

    private static Call findEntry( final List<Integer> list1, final HetSite unphasedLoc ) {
        Call result = null;
        for ( Map.Entry<Call, List<Integer>> entry2 : unphasedLoc.getCallMap().entrySet() ) {
            final List<Integer> list2 = entry2.getValue();
            final int minCommon = Math.round(0.75F * Math.max(list1.size(), list2.size()));
            if ( minCommon > Math.min(list1.size(), list2.size()) ) {
                continue;
            }
            if ( countCommonIDs(list1, list2) >= minCommon ) {
                if ( result == null ) {
                    result = entry2.getKey();
                } else {
                    return null;
                }
            }
        }
        return result;
    }

    private static int countCommonIDs( final List<Integer> list1, final List<Integer> list2 ) {
        final Iterator<Integer> itr1 = list1.iterator();
        int val1 = itr1.next();
        final Iterator<Integer> itr2 = list2.iterator();
        int val2 = itr2.next();
        int commonCount = 0;
        while ( true ) {
            int cmp = Integer.compare(val1, val2);
            if ( cmp == 0 ) {
                commonCount += 1;
                if ( !itr1.hasNext() || !itr2.hasNext() ) return commonCount;
                val1 = itr1.next();
                val2 = itr2.next();
            } else if ( cmp > 0 ) {
                if ( !itr2.hasNext() ) return commonCount;
                val2 = itr2.next();
            } else {
                if ( !itr1.hasNext() ) return commonCount;
                val1 = itr1.next();
            }
        }
        // can't reach here
    }

    private static List<ByteSequence> getHaplotypes( final int windowStart,
                                                     final List<Call> consensusCalls,
                                                     final List<HetSite> hetSites ) {
        if ( hetSites.isEmpty() ) {
            return Collections.singletonList(getHaplotype(windowStart, consensusCalls, hetSites.iterator(), (h)->null));
        }
        final List<ByteSequence> result = new ArrayList<>(2);
        result.add(getHaplotype(windowStart, consensusCalls, hetSites.iterator(), HetSite::getAlt1));
        result.add(getHaplotype(windowStart, consensusCalls, hetSites.iterator(), HetSite::getAlt2));
        return result;
    }

    public static ByteSequence getHaplotype( final int windowStart,
                                             final List<Call> consensusCalls,
                                             final Iterator<HetSite> hetSiteIterator,
                                             final Function<HetSite, Call> hetSiteCallFunc ) {
        final List<ByteSequence> seqs = new ArrayList<>(WINDOW_SIZE);
        HetSite curHetSite = hetSiteIterator.hasNext() ? hetSiteIterator.next() : null;
        final int nCalls = consensusCalls.size();
        for ( int idx = 0; idx < nCalls; ++idx ) {
            final Call call;
            if ( curHetSite == null || windowStart + idx != curHetSite.getRefLoc() ) {
                call = consensusCalls.get(idx);
            } else {
                call = hetSiteCallFunc.apply(curHetSite);
                if ( hetSiteIterator.hasNext() ) {
                    curHetSite = hetSiteIterator.next();
                }
            }
            if ( call == null ) {
                continue;
            }
            final CigarElement cigarElement = call.getCigarElement();
            if ( cigarElement.getOperator() == CigarOperator.D ) {
                idx += cigarElement.getLength() - 1;
            } else {
                final ByteSequence seq = call.getCalls();
                seqs.add(seq);
            }
        }
        return ByteSequence.concat(seqs);
    }

    public void writeHaplotypes( final VariantContext variant,
                                 final List<ByteSequence> haplotypes,
                                 final String genotype ) {
        boolean vcfIsHet = variant.getGenotype(0).isHet();
        final StringBuilder sb = new StringBuilder();
        sb.append(variant.getContig())
                .append('\t')
                .append(variant.getStart())
                .append('\t')
                .append(variant.getID())
                .append('\t')
                .append(genotype)
                .append('\t')
                .append(vcfIsHet ? "0/1" : "1/1");
        for ( final ByteSequence haplotype : haplotypes ) {
            sb.append('\t');
            sb.append(haplotype);
        }
        System.out.println(sb);
    }

    public record CallPair( Call call1, Call call2 ) {}

    public static final class HetSite {
        private final int refLoc;
        private final Map<Call, List<Integer>> callMap;
        private Call alt1;
        private Call alt2;

        public HetSite( int refLoc, Map<Call, List<Integer>> callMap ) {
            this.refLoc = refLoc;
            this.callMap = callMap;
        }

        public int getRefLoc() { return refLoc; }
        public Map<Call, List<Integer>> getCallMap() { return callMap; }
        public Call getAlt1() { return alt1; }
        public Call getAlt2() { return alt2; }
        public List<Integer> getListFor( final Call call ) { return callMap.get(call); }

        public void setPhase( final Call alt1, final Call alt2 ) {
            this.alt1 = alt1;
            this.alt2 = alt2;
        }
    }

    public static final class CigarMatcher {
        final CigarOperator operator;
        final int length;

        public CigarMatcher( final VariantContext variant ) {
            final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, null);
            final int svLen = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
            if ( GATKSVVCFConstants.SYMB_ALT_STRING_INS.equals(svType) ||
                    GATKSVVCFConstants.SYMB_ALT_STRING_DUP.equals(svType) ) {
                operator = CigarOperator.I;
                length = svLen + 1;
            } else if ( GATKSVVCFConstants.SYMB_ALT_STRING_DEL.equals(svType) ) {
                operator = CigarOperator.D;
                length = -svLen;
            } else if ( GATKSVVCFConstants.BREAKEND_STR.equals(svType) ||
                    GATKSVVCFConstants.SYMB_ALT_STRING_INV.equals(svType) ) {
                operator = CigarOperator.M;
                length = -1;
            } else {
                throw new UnsupportedOperationException("Can't handle svType " + svType);
            }
        }

        public boolean matches( final CigarElement element ) {
            return element.getOperator() == operator &&
                    (length < 0 || Math.abs(element.getLength() - length) <= 2);
        }
    }
}
