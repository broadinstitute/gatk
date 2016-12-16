package org.broadinstitute.hellbender.transformers;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.*;

/**
 * Masks read bases and base qualities using the symmetric DUST algorithm
 */
public class DUSTReadTransformer implements ReadTransformer {
    public static final long serialVersionUID = 1L;
    private static final Logger logger = LogManager.getLogger(DUSTReadTransformer.class);

    private int MASK_PHRED = 15; //Phred score to mask low-complexity sequences with
    private int W = 64; //DUST window size
    private double T = 20.0f; //DUST score threshold

    private static final class PStruct {
        public int start = 0;
        public int finish = 0;
        public double score = 0;
    }

    private static final class IntTuple {
        public int a;
        public int b;
        public IntTuple(final int a, final int b) {this.a = a; this.b = b;}
    }


    public DUSTReadTransformer(final int dust_mask, final int dust_w, final double dust_t) {
        MASK_PHRED = dust_mask;
        W = dust_w;
        T = dust_t;
    }

    @Override
    public GATKRead apply(final GATKRead read) {

        List<IntTuple> res = runDUST(read,W,T);

        //Mask the intervals and base qualities
        byte[] q_new = read.getBases();
        for (IntTuple interval : res) {
            for (int i = interval.a; i <= interval.b; i++) {
                q_new[i] = 'N';
            }
        }
        read.setBases(q_new);

        if (read.getBaseQualityCount() == read.getLength()) {
            byte[] qual_new = read.getBaseQualities();
            for (IntTuple interval : res) {
                for (int i = interval.a; i <= interval.b; i++) {
                    qual_new[i] = (byte)MASK_PHRED;
                }
            }
            read.setBaseQualities(qual_new);
        }

        return read;
    }

    /**
     * Wrapper method for runDUST so that we don't have to expose the IntTuple class.
     */
    public static Collection<Tuple2<Integer,Integer>> sDUST(final GATKRead read, final int window, final double tscore) {
        Collection<IntTuple> res = runDUST(read,window,tscore);
        List<Tuple2<Integer,Integer>> res_tuple2 = new ArrayList<>(res.size());
        for (IntTuple tuple : res) {
            res_tuple2.add(new Tuple2<>(tuple.a,tuple.b));
        }
        return res_tuple2;
    }

    /**
     *   Implementation of the sDUST algorithm described in:
     //  Liebert et al.(2006). A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA Sequences.
     //  Journal of Computational Biology, 13(5), 1028 - 1040.
     */
    private static List<IntTuple> runDUST(final GATKRead read, final int window, final double tscore) {

        final byte[] q = read.getBases();

        //Find convert to 2-bit and assign random bases to N's
        Random rng = null;

        for (int i = 0; i < q.length; i++) {
            if (q[i] == 'N' || q[i] == 'n') {
                if (rng == null) {rng = new Random(); rng.setSeed(0);}
                q[i] = (byte)rng.nextInt(4);
            } else {
                if (q[i] == 'A' || q[i] == 'a') {
                    q[i] = 0x0;
                } else if (q[i] == 'T' || q[i] == 't') {
                    q[i] = 0x1;
                } else if (q[i] == 'C' || q[i] == 'c') {
                    q[i] = 0x2;
                } else if (q[i] == 'G' || q[i] == 'g') {
                    q[i] = 0x3;
                } else {
                    logger.warn("Invalid base " + (char) q[i] + " in read " + read.getName() + ", assigning to random base.");
                    if (rng == null) {rng = new Random(); rng.setSeed(0);}
                    q[i] = (byte)rng.nextInt(4);
                }
            }
        }

        //Main sDUST algorithm as written in the pseudocode at the end of the paper
        //  Subroutines ADD_TRIPLET_INFO and REMOVE_TRIPLET_INFO are converted to 1-line statements
        //  SAVE_MASKED_REGIONS is implemented as a separate function
        //  SHIFT_WINDOW and FIND_PERFECT are written in place

        //List of low-complexity base intervals ("res" in pseudocode)
        final ArrayList<IntTuple> final_intervals = new ArrayList<>();

        //List of perfect intervals ("P"). A perfect interval is a subsequence of the read at most length window and with
        //score t > tscore whose subsequences all have score t_sub <= t. The list is maintained in sorted order, first
        //ascending by start position, then descending by end position.
        final LinkedList<PStruct> perfect_intervals = new LinkedList<>();

        //Triplet counts in the largest suffix, v, of the current window, w ("cv[t]" in pseudocode)
        //cv[i] corresponds to the count of triplet i, where i is the decimal value of the concatenation of the 2-bit
        //represenations of the 3 bases (64 total possible, see the triplet() function)
        final short[] triplet_count_suffix = new short[64];

        //Triplet counts in w, analogous to cv ("cw[t]")
        final short[] triplet_count_window = new short[64];

        //Running score counts for w and v ("rw" and "rv")
        int running_score_suffix = 0;
        int running_score_window = 0;

        //List of triplets in w ("w")
        final LinkedList<Integer> window_triplets = new LinkedList<>();

        //Number of triplets in v occuring at most 2*tscore times (Proposition 1 in manuscript)
        int L = 0;

        //Window start position
        int wstart;

        //Loop over all windows in the read
        for (int wfinish = 2; wfinish < read.getLength(); wfinish++) {
            wstart = Math.max(wfinish - window + 1, 0);

            save_masked_regions(perfect_intervals, final_intervals, wstart);

            //Triplet at end of the window
            int t = triplet(q[wfinish-2],q[wfinish-1],q[wfinish]);

            //Updates triplet counts and running scores of w and v, window_triplets, and L
            //Begin: SHIFT_WINDOW
            int s;
            if (window_triplets.size() >= window - 2) {
                s = window_triplets.pop();
                running_score_window -= --triplet_count_window[s]; //REMOVE_TRIPLET_INFO
                if (L > window_triplets.size()) {
                    L--;
                    running_score_suffix -= --triplet_count_suffix[s]; //REMOVE_TRIPLET_INFO
                }
            }
            window_triplets.add(t);
            L++;
            running_score_window += triplet_count_window[t]++; //ADD_TRIPLET_INFO
            running_score_suffix += triplet_count_suffix[t]++; //ADD_TRIPLET_INFO
            if (triplet_count_suffix[t]*10 > 2*tscore) {
                do {
                    s = window_triplets.get(window_triplets.size()-L);
                    running_score_suffix -= --triplet_count_suffix[s]; //REMOVE_TRIPLET_INFO
                    L--;
                } while (s != t);
            }
            //End: SHIFT_WINDOW

            if (running_score_window*10 > L*tscore) {

                //Adds perfect intervals in w that are suffixes of w to perfect_intervals
                //Begin: FIND_PERFECT
                final short[] c = triplet_count_suffix.clone();
                int r = running_score_suffix;
                ListIterator<PStruct> iter = perfect_intervals.listIterator();
                PStruct perf = iter.hasNext() ? iter.next() : null;
                double max_score = 0;
                for (int i = window_triplets.size() - L - 1; i >= 0; i--) {
                    t = window_triplets.get(i);
                    r += c[t];
                    c[t]++;
                    final double new_score = ((double)r)/(window_triplets.size() - i - 1);
                    if (new_score*10 > tscore) {
                        while (perf != null && perf.start >= i + wstart) {
                            max_score = Math.max(max_score,perf.score);
                            perf = iter.hasNext() ? iter.next() : null;
                        }
                        if (new_score >= max_score) {
                            max_score = new_score;
                            PStruct new_perf = new PStruct();
                            new_perf.start = i + wstart;
                            new_perf.finish = window_triplets.size() + 1 + wstart;
                            new_perf.score = new_score;
                            if (iter.hasPrevious()) {
                                iter.previous();
                            } else {
                                iter = perfect_intervals.listIterator();
                            }
                            iter.add(new_perf);
                            if (iter.hasNext()) {iter.next();}
                        }
                    }
                }
                //End: FIND_PERFECT

            }
        } //End of main loop


        //Save masked regions from final window
        wstart = Math.max(0,read.getLength() - window + 1);
        while (!perfect_intervals.isEmpty()) {
            save_masked_regions(perfect_intervals, final_intervals, wstart);
            wstart++;
        }

        return final_intervals;
    }

    /**
     * Saves masked intervals found in the previous window that do not overlap the current window.
     * Note that the LinkedList and ArrayList are modified.
     */
    private static void save_masked_regions(final LinkedList<PStruct> perfect_intervals,
                                            final ArrayList<IntTuple> final_intervals,
                                            final int wstart) {
        if (!perfect_intervals.isEmpty()) {
            PStruct p1 = perfect_intervals.getLast();
            if (p1.start < wstart) {
                final int l = final_intervals.size();
                if (l > 0) {
                    IntTuple rt = final_intervals.get(l - 1);
                    if (p1.start < rt.b + 1) {
                        final IntTuple new_rt = new IntTuple(rt.a, Math.max(rt.b, p1.finish));
                        final_intervals.set(l - 1, new_rt);
                    } else {
                        final_intervals.add(new IntTuple(p1.start, p1.finish));
                    }
                } else {
                    final_intervals.add(new IntTuple(p1.start, p1.finish));
                }
                while (p1 != null && p1.start < wstart) {
                    perfect_intervals.removeLast();
                    p1 = perfect_intervals.isEmpty() ? null : perfect_intervals.getLast();
                }
            }
        }
    }

    /**
     * Returns integer 0-63 of a given triplet of 2-bit bases
     */
    private static int triplet(final byte b1, final byte b2, final byte b3) {
        return (b1 | (b2 << 2) | (b3 << 4));
    }

}
