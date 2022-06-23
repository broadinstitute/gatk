package org.broadinstitute.hellbender.tools.walkers.longreads;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import shaded.cloud_nio.com.google.common.math.StatsAccumulator;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * Compute simple metrics on long read data
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>An unaligned or aligned BAM file</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 *     <li>A table with metrics</li>
 * </ul>
 *
 * <h3>Usage Example</h3>
 * <h4>Restore annotations from unaligned BAM files that may be discarded during conversion by samtools fastq</h4>
 * <pre>
 *   gatk ComputeLongReadMetrics \
 *     -I my.bam \
 *     -O metrics
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Compute simple metrics on long read data",
        oneLineSummary = "Compute simple metrics on long read data",
        programGroup = ReadDataManipulationProgramGroup.class
)
public final class ComputeLongReadMetrics extends ReadWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write output to this file prefix")
    public String OUTPUT;

    @Argument(fullName = "roundLengthsToNearest", shortName="rn", doc="Round read lengths to nearest <value>")
    public Integer ROUND_TO_NEAREST = 100;

    private long num_reads = 0;
    private long yield = 0;
    private long yield_times_passes = 0;

    private Map<Integer, Long> prl_hist = new TreeMap<>();
    private Map<Integer, Long> prl_yield_hist = new TreeMap<>();
    private Map<Integer, Long> rl_hist = new TreeMap<>();
    private Map<Integer, Long> rl_yield_hist = new TreeMap<>();
    private Map<Long, Long> zmw_num_hist = new TreeMap<>();
    private Map<Integer, Long> np_hist = new TreeMap<>();
    private Map<Integer, Long> range_gap_hist = new TreeMap<>();

    private List<Integer> polymerase_read_lengths = new ArrayList<>();
    private List<Integer> read_lengths = new ArrayList<>();

    private long last_zmw = -1;
    private int prl = 0;
    private int last_range_stop = 0;

    @Override
    public void apply( GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        SAMRecord sr = read.convertToSAMRecord(this.getHeaderForReads());

        if (!sr.isSecondaryOrSupplementary() || sr.getReadUnmappedFlag()) {
            int np = sr.hasAttribute("np") ? sr.getIntegerAttribute("np") : 1;
            int rl = sr.getReadLength();
            long zmw = last_zmw + 1;

            int range_start = 0;
            int range_stop = 0;

            if (sr.getReadName().contains("/")) {
                String[] pieces = sr.getReadName().split("/");
                if (pieces.length == 3) {
                    zmw = Long.valueOf(pieces[1]);

                    if (pieces[2].contains("_")) {
                        String[] range = pieces[2].split("_");
                        range_start = Integer.valueOf(range[0]);
                        range_stop = Integer.valueOf(range[1]);
                    }
                }
            }

            num_reads++;
            this.yield += rl;
            yield_times_passes += rl * np;

            increment(rl_hist, roundToNearest(rl, ROUND_TO_NEAREST));
            increment(zmw_num_hist, zmw);
            increment(np_hist, np);

            add(rl_yield_hist, roundToNearest(rl, ROUND_TO_NEAREST), (long) rl);

            if (last_range_stop != 0 && last_zmw == zmw) {
                increment(range_gap_hist, range_start - last_range_stop);
            }

            read_lengths.add(rl);

            if (last_zmw != -1 && last_zmw != zmw) {
                increment(prl_hist, roundToNearest(prl, ROUND_TO_NEAREST));
                add(prl_yield_hist, roundToNearest(prl, ROUND_TO_NEAREST), (long) prl);

                polymerase_read_lengths.add(prl);

                prl = 0;
            }

            prl = (zmw == last_zmw) ? prl + rl : rl;

            last_zmw = zmw;
            last_range_stop = range_stop;
        }
    }

    private int roundToNearest(int num, int nearest) {
        return ((num + (nearest - 1)) / nearest) * nearest;
    }

    private void increment(Map<Long, Long> h, Long k) {
        if (!h.containsKey(k)) {
            h.put(k, 1L);
        } else {
            h.put(k, h.get(k) + 1L);
        }
    }

    private void increment(Map<Integer, Long> h, Integer k) {
        if (!h.containsKey(k)) {
            h.put(k, 1L);
        } else {
            h.put(k, h.get(k) + 1L);
        }
    }

    private void add(Map<Integer, Long> h, Integer k, Long v) {
        if (!h.containsKey(k)) {
            h.put(k, v);
        } else {
            h.put(k, h.get(k) + v);
        }
    }

    private int NX(List<Integer> l, int X) {
        double sum = 0.0;
        for (int i : l) { sum += i; }

        double h = (((double) X)/100.0)*sum;

        double cs = 0.0;
        int last = l.get(0);
        for (int i : l) {
            cs += i;
            if (cs > h) {
                break;
            }

            last = i;
        }

        return last;
    }

    @Override
    public void closeTool() {
        Collections.sort(read_lengths, Collections.reverseOrder());
        Collections.sort(polymerase_read_lengths, Collections.reverseOrder());

        StatsAccumulator rsa = new StatsAccumulator();
        rsa.addAll(read_lengths);

        StatsAccumulator psa = new StatsAccumulator();
        psa.addAll(polymerase_read_lengths);

        PrintStream prl_counts_out;
        PrintStream prl_hist_out;
        PrintStream prl_yield_hist_out;
        PrintStream prl_nx_out;

        PrintStream rl_counts_out;
        PrintStream rl_hist_out;
        PrintStream rl_yield_hist_out;
        PrintStream rl_nx_out;

        PrintStream zmw_hist_out;
        PrintStream np_hist_out;
        PrintStream range_gap_hist_out;

        // TODO: add data for number of passes vs. read length boxplot

        try {
            prl_counts_out = new PrintStream(new FileOutputStream(OUTPUT + ".prl_counts.txt"));
            prl_hist_out = new PrintStream(new FileOutputStream(OUTPUT + ".prl_hist.txt"));
            prl_yield_hist_out = new PrintStream(new FileOutputStream(OUTPUT + ".prl_yield_hist.txt"));
            prl_nx_out = new PrintStream(new FileOutputStream(OUTPUT + ".prl_nx.txt"));

            rl_counts_out = new PrintStream(new FileOutputStream(OUTPUT + ".rl_counts.txt"));
            rl_hist_out = new PrintStream(new FileOutputStream(OUTPUT + ".rl_hist.txt"));
            rl_yield_hist_out = new PrintStream(new FileOutputStream(OUTPUT + ".rl_yield_hist.txt"));
            rl_nx_out = new PrintStream(new FileOutputStream(OUTPUT + ".rl_nx.txt"));

            zmw_hist_out = new PrintStream(new FileOutputStream(OUTPUT + ".zmw_hist.txt"));
            np_hist_out = new PrintStream(new FileOutputStream(OUTPUT + ".np_hist.txt"));
            range_gap_hist_out = new PrintStream(new FileOutputStream(OUTPUT + ".range_gap_hist.txt"));
        } catch (FileNotFoundException e) {
            throw new GATKException("Couldn't open output file.");
        }

        Map<Long, Long> zmw_hist = new TreeMap<>();
        for (long zmw : zmw_num_hist.keySet()) {
            increment(zmw_hist, zmw_num_hist.get(zmw));
        }

        prl_counts_out.println(Joiner.on("\t").join( "num_polymerase_reads", "polymerase_read_yield_bp", "polymerase_read_yield_bp_times_passes", "polymerase_read_min_length_bp", "polymerase_read_max_length_bp", "polymerase_read_mean_length_bp", "polymerase_read_sd_length_bp" ));
        prl_counts_out.println(Joiner.on("\t").join( psa.count(), psa.sum(), yield_times_passes, psa.min(), psa.max(), psa.mean(), psa.sampleStandardDeviation() ));

        prl_hist_out.println("polymerase_read_length\tlength_bp");
        for (int i : prl_hist.keySet()) {
            prl_hist_out.println(i + "\t" + prl_hist.get(i));
        }

        prl_yield_hist_out.println("polymerase_read_length\tyield_bp");
        for (int i : prl_yield_hist.keySet()) {
            prl_yield_hist_out.println(i + "\t" + prl_yield_hist.get(i));
        }

        prl_nx_out.println("NX\tvalue");
        rl_nx_out.println("NX\tvalue");
        for (int X = 1; X < 100; X++) {
            prl_nx_out.println(X + "\t" + NX(polymerase_read_lengths, X));
            rl_nx_out.println(X + "\t" + NX(read_lengths, X));
        }

        rl_counts_out.println(Joiner.on("\t").join( "num_reads", "read_yield_bp", "read_yield_bp_times_passes", "read_min_length_bp", "read_max_length_bp", "read_mean_length_bp", "read_sd_length_bp" ));
        rl_counts_out.println(Joiner.on("\t").join( rsa.count(), rsa.sum(), yield_times_passes, rsa.min(), rsa.max(), rsa.mean(), rsa.sampleStandardDeviation() ));

        rl_hist_out.println("read_length_bp\tcount");
        for (int i : rl_hist.keySet()) {
            rl_hist_out.println(i + "\t" + rl_hist.get(i));
        }

        rl_yield_hist_out.println("read_length\tyield_bp");
        for (int i : rl_yield_hist.keySet()) {
            rl_yield_hist_out.println(i + "\t" + rl_yield_hist.get(i));
        }

        zmw_hist_out.println("zmw_count\tcount");
        for (long i : zmw_hist.keySet()) {
            zmw_hist_out.println(i + "\t" + zmw_hist.get(i));
        }

        np_hist_out.println("np\tcount");
        for (int i : np_hist.keySet()) {
            np_hist_out.println(i + "\t" + np_hist.get(i));
        }

        range_gap_hist_out.println("range_gap\tcount");
        for (int i : range_gap_hist.keySet()) {
            range_gap_hist_out.println(i + "\t" + range_gap_hist.get(i));
        }

        prl_counts_out.close();
        prl_hist_out.close();
        prl_yield_hist_out.close();
        prl_nx_out.close();
        rl_counts_out.close();
        rl_hist_out.close();
        rl_yield_hist_out.close();
        rl_nx_out.close();
        zmw_hist_out.close();
        np_hist_out.close();
        range_gap_hist_out.close();
    }
}
