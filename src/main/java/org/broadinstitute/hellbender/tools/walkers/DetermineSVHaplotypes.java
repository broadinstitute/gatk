package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.ErrorCorrectHiFi.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.*;

@CommandLineProgramProperties(
        summary = "Assemble long read fragments near SV breakpoint.",
        oneLineSummary = "Assemble long read fragments near SV breakpoint.",
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
        final int variantStart = variant.getStart() + 1; // vcf references base before the actual variant for indels
        final int paddedStart = Math.max(1, variantStart - PADDING);
        final SimpleInterval window = new SimpleInterval(variant.getContig(), paddedStart, paddedStart);
        final List<CallIterator> callIterators = new ArrayList<>();
        readsContext.iterator(window).forEachRemaining(read -> callIterators.add(new CallIterator(read, paddedStart)));
        final int paddedEnd = paddedStart + WINDOW_SIZE;
        final int nReads = callIterators.size();
        Map<Call, List<Integer>> callMap = new HashMap<>(nReads);
        final List<VariantLoc> variantLocs = new ArrayList<>();
        final List<Call> consensusCalls = new ArrayList<>(WINDOW_SIZE);
        for ( int refLoc = paddedStart; refLoc < paddedEnd; ++refLoc ) {
            for ( int readIdx = 0; readIdx != nReads; ++ readIdx ) {
                final CallIterator callIterator = callIterators.get(readIdx);
                final Call call = callIterator.hasNext() ? callIterator.next() : null;
                if ( call != null && call.getMeanQual() >= MIN_QUAL ) {
                    final int callIdx = readIdx;
                    callMap.compute(call, (k, v) -> {
                        if ( v == null ) v = new ArrayList<>();
                        v.add(callIdx);
                        return v;
                    });
                }
            }
            if ( isVariant(callMap, consensusCalls) ) {
                variantLocs.add(new VariantLoc(refLoc, callMap));

                // make a new collection, the old one has been stashed in the variantLoc
                callMap = new HashMap<>(nReads);
            } else {
                // reuse the old collection after clearing it
                callMap.clear();
            }
        }
        phase(variantLocs, variantStart);
        System.out.println(variant.getContig() + ":" + variant.getStart() + "\t" + variant.getID() + "\t" + variantLocs.size());
    }

    private boolean isVariant( final Map<Call, List<Integer>> callMap,
                               final List<Call> consensusCalls ) {
        int maxIDs = 0;
        final int totalIds = callMap.values().stream().mapToInt(List::size).sum();
        Map.Entry<Call, List<Integer>> maxEntry = null;
        final Iterator<Map.Entry<Call, List<Integer>>> itr = callMap.entrySet().iterator();
        while ( itr.hasNext() ) {
            final Map.Entry<Call, List<Integer>> entry = itr.next();
            final List<Integer> callIndices = entry.getValue();
            final int nIDs = callIndices.size();
            if ( nIDs == 1 ) { // assume that variation seen in a single read is bogus
                itr.remove();
                continue;
            }
            // remove small indels that aren't well attested -- that's our main sequencing error mode
            final Call call = entry.getKey();
            final int keyLen = call.getCalls().length();
            if ( keyLen == 2 || (keyLen == 0 && call.getCigarElement().getLength() <= 2) ) {
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

    private void phase( final List<VariantLoc> variantLocs, final int variantStart ) {
        final int nLocs = variantLocs.size();
        for ( int idx = 0; idx != nLocs; ++idx ) {
            final VariantLoc variantLoc = variantLocs.get(idx);
            if ( variantLoc.getRefLoc() == variantStart || nLocs == 1 ) {
                final ArrayList<Map.Entry<Call, List<Integer>>> locs =
                        new ArrayList<>(variantLoc.getCallMap().entrySet());
                locs.sort(Comparator.comparingInt(e -> -e.getValue().size()));
                variantLoc.setPhase(locs.get(0).getKey(), locs.get(1).getKey());
                phaseExtend(variantLocs, idx, idx);
                return;
            }
            if ( variantLoc.getRefLoc() > variantStart ) {
                break;
            }
        }
        for ( int idx1 = 0; idx1 < nLocs - 1; ++idx1 ) {
            for ( int idx2 = idx1 + 1; idx2 < nLocs; ++idx2 ) {
                if ( phasePair(variantLocs.get(idx1), variantLocs.get(idx2)) ) {
                    phaseExtend(variantLocs, idx1, idx2);
                    return;
                }
            }
        }

        // if we have only two variants, and they can't be phased together
        if ( nLocs == 2 ) {
            // remove the one that's farther from the variant of interest
            if ( Math.abs(variantLocs.get(0).getRefLoc() - variantStart) <
                    Math.abs(variantLocs.get(1).getRefLoc() - variantStart) ) {
                variantLocs.remove(1);
            } else {
                variantLocs.remove(0);
            }
            // arbitrarily set the phase of the one that remains
            final VariantLoc variantLoc = variantLocs.get(0);
            final ArrayList<Map.Entry<Call, List<Integer>>> locs =
                    new ArrayList<>(variantLoc.getCallMap().entrySet());
            locs.sort(Comparator.comparingInt(e -> -e.getValue().size()));
            variantLoc.setPhase(locs.get(0).getKey(), locs.get(1).getKey());
        }
        variantLocs.clear();
    }

    private boolean phasePair( final VariantLoc loc1, final VariantLoc loc2 ) {
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

    private void phaseExtend( final List<VariantLoc> variantLocs, final int idx1, int idx2 ) {
        // phase or delete variants downstream of idx2
        VariantLoc phasedLoc = variantLocs.get(idx2);
        int idx3 = idx2 + 1;
        while ( idx3 < variantLocs.size() ) {
            final VariantLoc unphasedLoc = variantLocs.get(idx3);
            if ( extendPhase(phasedLoc, unphasedLoc) ) {
                phasedLoc = unphasedLoc;
                idx3 += 1;
            } else {
                variantLocs.remove(idx3);
            }
        }

        // phase or delete variants between idx1 and idx2
        phasedLoc = variantLocs.get(idx1);
        int idx4 = idx1 + 1;
        while ( idx4 < idx2 ) {
            final VariantLoc unphasedLoc = variantLocs.get(idx4);
            if ( extendPhase(phasedLoc, unphasedLoc) ) {
                phasedLoc = unphasedLoc;
                idx4 += 1;
            } else {
                variantLocs.remove(idx4);
                idx2 -= 1;
            }
        }

        // phase or delete variants upstream of idx1
        phasedLoc = variantLocs.get(idx1);
        int idx0 = idx1 - 1;
        while ( idx0 >= 0 ) {
            final VariantLoc unphasedLoc = variantLocs.get(idx0);
            if ( extendPhase(phasedLoc, unphasedLoc) ) {
                phasedLoc = unphasedLoc;
            } else {
                variantLocs.remove(idx0);
            }
            idx0 -= 1;
        }
    }

    private boolean extendPhase( final VariantLoc phasedLoc, final VariantLoc unphasedLoc ) {
        final Call call1 = findEntry(phasedLoc.getAlt1(), unphasedLoc);
        if ( call1 == null ) {
            return false;
        }
        final Call call2 = findEntry(phasedLoc.getAlt2(), unphasedLoc);
        if ( call2 == null || call1.equals(call2) ) {
            return false;
        }
        unphasedLoc.setPhase(call1, call2);
        return true;
    }

    private Call findEntry( final List<Integer> list1, final VariantLoc unphasedLoc ) {
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

    private int countCommonIDs( final List<Integer> list1, final List<Integer> list2 ) {
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

    public record CallPair( Call call1, Call call2 ) {}

    public static final class VariantLoc {
        private final int refLoc;
        private final Map<Call, List<Integer>> callMap;
        private Call alt1;
        private Call alt2;

        public VariantLoc( int refLoc, Map<Call, List<Integer>> callMap ) {
            this.refLoc = refLoc;
            this.callMap = callMap;
        }

        public int getRefLoc() { return refLoc; }
        public Map<Call, List<Integer>> getCallMap() { return callMap; }
        public List<Integer> getAlt1() { return alt1 == null ? null : callMap.get(alt1); }
        public List<Integer> getAlt2() { return alt1 == null ? null : callMap.get(alt2); }

        public void setPhase( final Call alt1, final Call alt2 ) {
            this.alt1 = alt1;
            this.alt2 = alt2;
        }
    }
}
