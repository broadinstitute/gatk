package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import com.opencsv.CSVReader;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.util.Precision;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.FlowBasedAlignmentArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowFeatureMapper;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.FlowBasedAlignmentLikelihoodEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.LikelihoodEngineArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.report.GATKReport;
import org.broadinstitute.hellbender.utils.report.GATKReportTable;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.zip.GZIPOutputStream;

@CommandLineProgramProperties(
        summary = "Ground Truth Scorer",
        oneLineSummary = "Score reads against a reference/ground truth",
        programGroup = FlowBasedProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class GroundTruthScorer extends ReadWalker {
    private static final Logger logger = LogManager.getLogger(GroundTruthScorer.class);

    public static final String OUTPUT_CSV_LONG_NAME = "output-csv";
    public static final String REPORT_FILE_LONG_NAME = "report-file";
    public static final String USE_SOFTCLIPPED_BASES_LONG_NAME = "use-softclipped-bases";
    public static final String GENOME_PRIOR_LONG_NAME = "genome-prior";
    public static final String FEATURES_FILE_LONG_NAME = "features-file";
    public static final String NORMALIZED_SCORE_THRESHOLD_LONG_NAME = "normalized-score-threshold";
    public static final String ADD_MEAN_CALL_LONG_NAME = "add-mean-call";

    private static final int QUAL_VALUE_MAX = 60;
    private static final int HMER_VALUE_MAX = 30; //TODO: This should become a parameter
    private static final int BASE_VALUE_MAX = FlowBasedRead.DEFAULT_FLOW_ORDER.length() - 1;

    private static final double NORMALIZED_SCORE_THRESHOLD_DEFAULT = -0.1;

    private static class BooleanAccumulator {
        long falseCount;
        long trueCount;
        BooleanAccumulator[] bins;

        void add(final boolean b) {
            if (b) {
                trueCount++;
            } else {
                falseCount++;
            }
        }

        void add(final boolean b, final int bin) {
            add(b);
            if ( bins != null && bin >= 0 && bin < bins.length ) {
                bins[bin].add(b);
            }
        }
        void add(final boolean b, final int bin, final int bin2) {
            add(b);
            if ( bins != null && bin >= 0 && bin < bins.length ) {
                bins[bin].add(b, bin2);
            }
        }
        void add(final boolean b, final int bin, final int bin2, final int bin3) {
            add(b);
            if ( bins != null && bin >= 0 && bin < bins.length ) {
                bins[bin].add(b, bin2, bin3);
            }
        }

        long getCount() {
            return falseCount + trueCount;
        }

        double getFalseRate() {
            return (getCount() == 0) ? 0.0 : ((double)falseCount / getCount());
        }

        static BooleanAccumulator[] newReport(final int size, final int binCount, final int binCount2, final int binCount3) {
            BooleanAccumulator[]   report = new BooleanAccumulator[size];
            for ( byte i = 0 ; i < report.length ; i++ ) {
                report[i] = new BooleanAccumulator();
                if ( binCount != 0 ) {
                    report[i].bins = new BooleanAccumulator[binCount];
                    for ( int j = 0 ; j < report[i].bins.length ; j++ ) {
                        report[i].bins[j] = new BooleanAccumulator();
                        if ( binCount2 != 0 ) {
                            report[i].bins[j].bins = new BooleanAccumulator[binCount2];
                            for ( int k = 0 ; k < report[i].bins[j].bins.length ; k++ ) {
                                report[i].bins[j].bins[k] = new BooleanAccumulator();
                                if (binCount3 != 0) {
                                    report[i].bins[j].bins[k].bins = new BooleanAccumulator[binCount3];
                                    for (int m = 0; m < report[i].bins[j].bins[k].bins.length; m++) {
                                        report[i].bins[j].bins[k].bins[m] = new BooleanAccumulator();
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return report;
        }

        static GATKReportTable newReportTable(final BooleanAccumulator[] report, final String name, final double probThreshold) {
            final GATKReportTable table = new GATKReportTable(name + "Report", "error rate per " + name, 4);
            table.addColumn(name, "%d");
            table.addColumn("count", "%d");
            table.addColumn("error", "%f");
            table.addColumn("phred", "%d");
            for (int i = 0; i < report.length; i++) {
                final double rate =  report[i].getFalseRate();
                final double phredRate = (rate == 0.0 && report[i].getCount() != 0 && probThreshold != 0.0 ) ? probThreshold : rate;

                table.set(i, 0, i);
                table.set(i, 1, report[i].getCount());
                table.set(i, 2, rate);
                table.set(i, 3, phredRate != 0 ? (int)Math.ceil(-10.0 * Math.log10(phredRate)) : 0);
            }
            return table;
        }

        static GATKReportTable newReportTable(final BooleanAccumulator[] report, final String name1, final String name2) {
            final GATKReportTable table = new GATKReportTable(name1 + "_" + name2 + "Report", "error rate per " + name1 + " by " + name2, 4);
            table.addColumn(name1, "%d");
            table.addColumn(name2, "%d");
            table.addColumn("count", "%d");
            table.addColumn("error", "%f");
            int rowIndex = 0;
            for (int i = 0; i < report.length; i++) {
                for ( int j = 0; j < report[i].bins.length ; j++ ) {
                    table.set(rowIndex, 0, i);
                    table.set(rowIndex, 1, j);
                    table.set(rowIndex, 2, report[i].bins[j].getCount());
                    table.set(rowIndex, 3, report[i].bins[j].getFalseRate());
                    rowIndex++;
                }
            }
            return table;
        }

        static GATKReportTable newReportTable(final BooleanAccumulator[] report, final String name1, final String name2, final String name3, final String name4) {
            final GATKReportTable table = new GATKReportTable(name1 + "_" + name2 + "_" + name3 + "_" + name4 + "_Report", "error rate per " + name1 + " by " + name2 + " and " + name3, 5);
            table.addColumn(name1, "%d");
            table.addColumn(name2, "%d");
            table.addColumn(name3, "%s");
            table.addColumn(name4, "%s");
            table.addColumn("count", "%d");
            int rowIndex = 0;
            for (int i = 0; i < report.length; i++) {
                for ( int j = 0; j < report[i].bins.length ; j++ ) {
                    for ( int k = 0; k < report[i].bins[j].bins.length ; k++ ) {
                        for ( int m = 0; m < report[i].bins[j].bins[k].bins.length ; m++ ) {
                            table.set(rowIndex, 0, i);
                            table.set(rowIndex, 1, j);
                            table.set(rowIndex, 2, binToDeviation(k));
                            table.set(rowIndex, 3, String.format("%c", binToBase(m)));
                            table.set(rowIndex, 4, report[i].bins[j].bins[k].bins[m].getCount());
                            rowIndex++;
                        }
                    }
                }
            }
            return table;
        }
    }


    @Argument(fullName = OUTPUT_CSV_LONG_NAME, doc="main CSV output file. supported file extensions: .csv, .csv.gz.")
    public GATKPath outputCsvPath = null;

    @Argument(fullName = REPORT_FILE_LONG_NAME, doc="report output file.", optional = true)
    public GATKPath reportFilePath = null;

    @ArgumentCollection
    public LikelihoodEngineArgumentCollection likelihoodArgs = new LikelihoodEngineArgumentCollection();

    @ArgumentCollection
    public FlowBasedAlignmentArgumentCollection fbargs = new FlowBasedAlignmentArgumentCollection();

    @Argument(fullName = USE_SOFTCLIPPED_BASES_LONG_NAME, doc="", optional = true)
    public boolean useSoftclippedBases;

    @Argument(fullName = GENOME_PRIOR_LONG_NAME, doc="CSV input file containing genome-prior (one line per base with hmer frequencies).", optional = true)
    public GATKPath genomePriorPath;

    @Argument(fullName = FEATURES_FILE_LONG_NAME, doc="A VCF file containing features to be used as a use for filtering reads.", optional = true)
    public FeatureDataSource<VariantContext> features;

    @Argument(fullName = NORMALIZED_SCORE_THRESHOLD_LONG_NAME, doc="threshold for normalized score, below which reads are ignored", optional = true)
    public double normalizedScoreThreshold = NORMALIZED_SCORE_THRESHOLD_DEFAULT;

    @Argument(fullName = ADD_MEAN_CALL_LONG_NAME, doc="Add ReadMeanCall and ReadProbs columns to output", optional = true)
    public boolean addMeanCalll;

    @Argument(fullName = "gt-no-output", doc = "do not generate output records", optional = true)
    public boolean      noOutput = false;


    // locals
    private FlowBasedAlignmentLikelihoodEngine likelihoodCalculationEngine;
    private PrintWriter                         outputCsv;
    private DecimalFormat                       doubleFormat = new DecimalFormat("0.0#####");
    private GenomePriorDB                       genomePriorDB;
    private BooleanAccumulator[]                qualReport;
    private String[]                            csvFieldOrder;

    // static/const
    static final private String[]       CSV_FIELD_ORDER_BASIC = {
            "ReadName", "ReadKey", "ReadIsReversed", "ReadMQ", "ReadRQ", "GroundTruthKey", "ReadSequence",
            "Score", "NormalizedScore", "ErrorProbability",
            "ReadKeyLength", "GroundTruthKeyLength", "CycleSkipStatus", "Cigar", "LowestQBaseTP"
    };
    static final private String[]       CSV_FIELD_ORDER_MEAN_CALL = {
            "ReadProbs", "ReadMeanCall"
    };

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        // establish csv fields
        List<String>        order = new LinkedList<>(Arrays.asList(CSV_FIELD_ORDER_BASIC));
        if ( addMeanCalll ) {
            order.addAll(Arrays.asList(CSV_FIELD_ORDER_MEAN_CALL));
        }
        csvFieldOrder = order.toArray(new String[0]);

        // create likelihood engine
        final ReadLikelihoodCalculationEngine engine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(likelihoodArgs, false);
        if ( engine instanceof FlowBasedAlignmentLikelihoodEngine ) {
            likelihoodCalculationEngine = (FlowBasedAlignmentLikelihoodEngine)engine;
        } else {
            throw new GATKException("must use a flow based likelihood calculation engine");
        }

        // open genome prior if provided
        if ( genomePriorPath != null ) {
            try {
                genomePriorDB = new GenomePriorDB(genomePriorPath);
            } catch (IOException e) {
                throw new GATKException("failed to open genome-prior file: " + genomePriorPath);
            }
        }

        // open output, write header
        try {
            if (outputCsvPath.toPath().toString().endsWith(".gz")) {
                outputCsv = new PrintWriter(new GZIPOutputStream(outputCsvPath.getOutputStream()));
            } else {
                outputCsv = new PrintWriter(outputCsvPath.getOutputStream());
            }
        } catch (IOException e) {
            throw new GATKException("failed to open csv output: " + outputCsvPath, e);
        }
        emitCsvHeaders();

        // initialize reports
        if ( reportFilePath != null ) {
            qualReport = BooleanAccumulator.newReport(QUAL_VALUE_MAX + 1, HMER_VALUE_MAX + 1, deviationToBin(HMER_VALUE_MAX + 1), BASE_VALUE_MAX + 1);
        }
    }

    @Override
    public void closeTool() {

        // close main output csv
        if ( outputCsv != null ) {
            outputCsv.close();
        }

        // write reports
        if ( reportFilePath != null ) {
            final GATKReport report = new GATKReport(
                    BooleanAccumulator.newReportTable(qualReport, "qual", fbargs.probabilityRatioThreshold),
                    BooleanAccumulator.newReportTable(qualReport, "qual", "hmer"),
                    BooleanAccumulator.newReportTable(qualReport, "qual", "hmer", "deviation", "base"));
            try ( final PrintStream ps = new PrintStream(reportFilePath.getOutputStream()) ) {
                report.print(ps);
            }
        }

        super.closeTool();
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // handle read clipping
        final GATKRead clippedRead;
        if (isSoftClipped(read) ) {
            if (useSoftclippedBases) {
                referenceContext.setWindow(read.getStart() - read.getUnclippedStart(), read.getUnclippedEnd() - read.getEnd());;
                clippedRead = ReadClipper.revertSoftClippedBases(read);
            } else {
                clippedRead = ReadClipper.hardClipSoftClippedBases(read);
            }
        } else {
            clippedRead = read;
        }

        // filter?
        if ( (features != null) && !filter(clippedRead, referenceContext) ) {
            return;
        }

        // create flow read/haplotype
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(getHeaderForReads(), clippedRead);
        final FlowBasedRead flowRead = new FlowBasedRead(clippedRead, rgInfo.flowOrder, rgInfo.maxClass, fbargs);
        final Haplotype haplotype = new Haplotype(referenceContext.getBases(), true);
        final FlowBasedHaplotype flowHaplotype = new FlowBasedHaplotype(haplotype, rgInfo.flowOrder);

        // is this really needed?
        if ( !flowRead.isValid() ) {
           return;
        }

        // compute score
        final int         hapKeyLength = flowHaplotype.getKeyLength();
        final double      score = FlowFeatureMapper.computeLikelihoodLocal(flowRead, flowHaplotype, hapKeyLength, false);
        final double      normalizedScore = score / flowRead.getKeyLength();
        if ( normalizedScore < normalizedScoreThreshold )
            return;

        // compute error probability
        final double[]    errorProb = computeErrorProb(flowRead, genomePriorDB);

        // cycle skip
        final FlowBasedReadUtils.CycleSkipStatus cycleSkipStatus = FlowBasedReadUtils.getCycleSkipStatus(flowRead, referenceContext);

        // accumulate reports
        if ( cycleSkipStatus != FlowBasedReadUtils.CycleSkipStatus.CS && qualReport != null ) {
            addToQualReport(flowRead, referenceContext, errorProb);
        }

        // emit
        try {
            emit(flowRead, flowHaplotype, score, normalizedScore, errorProb, read, cycleSkipStatus);
        } catch (IOException e) {
            throw new GATKException("failed to write output record", e);
        }
    }

    private boolean filter(final GATKRead read, final ReferenceContext referenceContext) {

        // loop on features contained in the read, check that they are in agreement with read data
        Iterator<VariantContext>    iter = features.query(new SimpleInterval(read));
        byte[]                      ref = referenceContext.getBases();
        while ( iter.hasNext() ) {
            final VariantContext vc = iter.next();
            for ( int refCoord = vc.getStart() ; refCoord <= vc.getEnd() ; refCoord++ ) {

                // get byte from read
                Optional<Byte>      readByte = ReadUtils.getReadBaseAtReferenceCoordinate(read, refCoord);
                if ( !readByte.isPresent() ) {
                    return false;
                }

                // get byte from reference
                byte                refByte = ref[refCoord - referenceContext.getWindow().getStart()];

                // compare
                if ( refByte != readByte.get() ) {
                    return false;
                }
            }
        }

        // if here, no interference detected
        return true;
    }

    /*
     * compute error probability vector for a read
     *
     * The vector has one element for each flow key, representing the probability complementing the call-probability to 1
     * This is further complicated by the optional presence of a genome-prior database, which provides factoring for
     * each hmer length (on a base basis)
     */
    private double[] computeErrorProb(final FlowBasedRead flowRead, final GenomePriorDB genomePriorDB) {

        final int[] key = flowRead.getKey();
        final byte[] flowOrder = flowRead.getFlowOrderArray();

        final double[] probCol = new double[flowRead.getMaxHmer() + 1];
        double[] result = new double[key.length];

        for ( int i = 0 ; i < key.length ; i++ ) {

            // step 1 - extract column & sum
            double  sum = 0;
            for ( int j = 0 ; j < probCol.length ; j++ ) {
                sum += (probCol[j] = flowRead.getProb(i, j));
            }

            // step 2 - normalize column
            for ( int j = 0 ; j < probCol.length ; j++ ) {
                probCol[j] /= sum;
            }

            // step 3 - scale according to prior genome?
            if ( genomePriorDB != null ) {

                long[] prior = (genomePriorDB != null) ? genomePriorDB.getPriorForBase(flowOrder[i]) : null;
                sum = 0;
                for ( int j = 0 ; j < probCol.length ; j++ ) {
                    sum += (probCol[j] *= prior[j]);
                }

                // assign normalized result
                result[i] = 1 - (probCol[Math.min(probCol.length - 1, Math.min(key[i], flowRead.getMaxHmer()))] / sum);
            } else {

                // assign normalized result
                result[i] = 1 - probCol[Math.min(key[i], flowRead.getMaxHmer())];
            }
        }

        return result;
    }

    /*
     * compute lowest-quality-base-tp-value vector for a read
     *
     * The vector has one element for each flow key, representing the value of the tp for the hmer base
     * which has the lowest quality value
     *
     * example:
     * Bases:  TTTTT
     * Qs:     ABCBA
     * Tp:     -1 1 0 1 -1
     * Output: -1
     */
    private byte[] computeLowestQBaseTP(final FlowBasedRead flowRead) {

        final int[] key = flowRead.getKey();
        byte[] result = new byte[key.length];
        final byte[] tp = flowRead.getAttributeAsByteArray("tp");
        final byte[] qual = flowRead.getBaseQualitiesNoCopy();

        // loop om key
        int seq_i = 0;
        for ( int i = 0 ; i < key.length ; i++ ) {

            // extract hmer length, zero is easy
            int hmer = key[i];
            if ( hmer == 0 ) {
                result[i] = 0;
                continue;
            }

            // scan qualities for the lowest value, start with first
            // as qualities and tp are symetric, we can scan up to the middle
            // when finding the middle (offset) account for even/odd hmers
            result[i] = tp[seq_i];
            byte lowestQ = qual[seq_i];
            int hmer_scan_length = (hmer + 1) / 2;
            for ( int j = 1 ; j < hmer_scan_length ; j++ ) {
                if ( qual[seq_i + j] < lowestQ ) {
                    result[i] = tp[seq_i + j];
                    lowestQ = qual[seq_i + j];
                }
            }

            // advance
            seq_i += hmer;
        }

        return result;
    }

    private void emitCsvHeaders() {

        outputCsv.println(StringUtils.join(csvFieldOrder, ","));
    }

    private void emit(final FlowBasedRead flowRead, final FlowBasedHaplotype refHaplotype, double score, final double normalizedScore, final double[] errorProb,
                      GATKRead read,
                      FlowBasedReadUtils.CycleSkipStatus cycleSkipStatus) throws IOException {

        // build line columns
        final Map<String,Object> cols = new LinkedHashMap<>();

        // read info
        cols.put("ReadName", flowRead.getName());
        cols.put("ReadIsReversed", flowRead.isReverseStrand() ? 1 : 0);
        cols.put("ReadMQ", flowRead.getMappingQuality());
        cols.put("ReadRQ", flowRead.getAttributeAsFloat("rq"));
        cols.put("CycleSkipStatus", cycleSkipStatus);
        cols.put("Cigar", read.getCigar().toString());


        // keys, seq, etc
        cols.put("ReadKey", "\"" + StringUtils.join(flowRead.getKey(), ',') + "\"");
        cols.put("GroundTruthKey", "\"" + StringUtils.join(refHaplotype.getKey(), ',') + "\"");
        cols.put("ReadSequence", flowRead.getBasesString());
        cols.put("ReadKeyLength", flowRead.getKeyLength());
        cols.put("GroundTruthKeyLength", refHaplotype.getKeyLength());

        // scores
        cols.put("Score", score);
        cols.put("NormalizedScore", normalizedScore);
        cols.put("ErrorProbability", "\"" + StringUtils.join(
                Arrays.stream(errorProb).mapToObj(v -> doubleFormat.format(v)).toArray(),
                ',') + "\"");

        // lowest q base tp
        final byte[]    lowestQBaseTP = computeLowestQBaseTP(flowRead);
        cols.put("LowestQBaseTP", "\"" + StringUtils.join(lowestQBaseTP, ',') + "\"");

        // add read probabilities
        if ( addMeanCalll ) {
            double[][] readProbsAndMeanCall = collectReadProbs(flowRead);
            cols.put("ReadProbs", "\"" + StringUtils.join(
                    Arrays.stream(readProbsAndMeanCall[0]).mapToObj(v -> doubleFormat.format(v)).toArray(),
                    ',') + "\"");
            cols.put("ReadMeanCall", "\"" + StringUtils.join(
                    Arrays.stream(readProbsAndMeanCall[1]).mapToObj(v -> doubleFormat.format(v)).toArray(),
                    ',') + "\"");
        }

        // construct line
        StringBuilder       sb = new StringBuilder();
        int                 colIndex = 0;
        for ( String field : csvFieldOrder ) {
            if ( colIndex++ > 0 ) {
                sb.append(',');
            }
            if ( !cols.containsKey(field) ) {
                throw new GATKException("column missing from csv line: " + field);
            }
            sb.append(cols.get(field));
            cols.remove(field);
        }
        if ( cols.size() > 0 ) {
            throw new GATKException("invalid columns on csv line: " + cols.keySet());
        }

        // output line
        if ( !noOutput ) {
            outputCsv.println(sb);
        }
    }

    private double[][] collectReadProbs(final FlowBasedRead read) {

        final int keyLength = read.getKeyLength();
        final int maxHmer = read.getMaxHmer();
        final double[] probs = new double[keyLength * (maxHmer + 1)];
        final double[] meanCall = new double[keyLength];

        // retrieve probs
        int pos = 0;
        for ( int flow = 0 ; flow < keyLength ; flow++ ) {
            double mc = 0;
            int mc_sum = 0;
            for ( int hmer = 0 ; hmer <= maxHmer ; hmer++ ) {
                final double p = read.getProb(flow, hmer);
                probs[pos++] = p;
                mc += (p * hmer);
                mc_sum += hmer;
            }
            meanCall[flow] = mc / mc_sum;
        }

        return new double[][] {probs, meanCall};
    }

    static private class GenomePriorDB {

        final private Map<Byte, long[]>      db = new LinkedHashMap<>();

        GenomePriorDB(GATKPath path) throws IOException {

            final CSVReader     csvReader = new CSVReader(new InputStreamReader(path.getInputStream()));
            String[]            line;
            while ( (line = csvReader.readNext()) != null ) {
                long[]          prior = new long[HMER_VALUE_MAX + 1];
                Byte            base = line[0].getBytes()[0];
                for ( int i = 0 ; i < prior.length ; i++ ) {
                    if ( i == 0 ){
                        prior[i] = Long.parseLong(line[i+1]);
                    } else {
                        prior[i] = Long.parseLong(line[i]);
                    }
                }
                db.put(base, prior);
            }
        }

        long[]   getPriorForBase(byte base) {
            return db.get(base);
        }
    }

    private void addToQualReport(FlowBasedRead flowRead, ReferenceContext referenceContext, final double[] errorProb) {

        // convert reference to key space
        final Haplotype             haplotype = new Haplotype(referenceContext.getBases(), true);
        final FlowBasedHaplotype    flowHaplotype = new FlowBasedHaplotype(haplotype, flowRead.getFlowOrder());

        // access keys & flow order
        final int[]                 readKey = flowRead.getKey();
        final int[]                 hapKey = flowHaplotype.getKey();
        final byte[]                flowOrder = flowRead.getFlowOrder().getBytes();

        // loop on key positions
        if ( readKey.length != hapKey.length ) {
            return;
        }
        for ( int flow = 0 ; flow < readKey.length ; flow++ ) {

            // determine quality
            final double        prob = Precision.round(errorProb[flow], (QUAL_VALUE_MAX / 10) + 1);
            final int           qual = (int)Math.ceil(-10 * Math.log10(prob));

            // determine if matches reference
            final int           deviation = readKey[flow] - hapKey[flow];
            final boolean       same = (deviation == 0);

            // accumulate
            if ( qual < qualReport.length ) {
                int         baseBin = baseToBin(flowOrder[flow % flowOrder.length], flowRead.isReverseStrand());
                qualReport[qual].add(same, readKey[flow], deviation, baseBin);
            }
        }
    }

    static private int baseToBin(byte base, boolean isReverseStrand) {
        final byte trueBase = !isReverseStrand ? base : BaseUtils.simpleComplement(base);
        return FlowBasedRead.DEFAULT_FLOW_ORDER.indexOf(trueBase);
    }

    static private byte binToBase(int bin) {
        return (byte)FlowBasedRead.DEFAULT_FLOW_ORDER.charAt(bin);
    }

    // 0,-1,1,-2,2... -> 0,1,2,3,4...
    static private int deviationToBin(final int deviation) {
        if ( deviation >= 0 ) {
            return deviation * 2;
        } else {
            return (-deviation * 2) - 1;
        }
    }

    // 0,1,2,3,4... -> 0,-1,1,-2,2...
    static private String binToDeviation(final int bin) {
        if ( bin == 0 ) {
            return "0";
        } else if ( (bin % 2) == 0 ) {
            return String.format("+%d", bin / 2);
        } else {
            return String.format("%d", -((bin+1) / 2));
        }
    }

    static private boolean isSoftClipped( final GATKRead read ) {
        final CigarOperator firstOperator = read.getCigar().getFirstCigarElement().getOperator();
        final CigarOperator lastOperator = read.getCigar().getLastCigarElement().getOperator();
        return (firstOperator == CigarOperator.SOFT_CLIP && lastOperator != CigarOperator.SOFT_CLIP) ||
                (firstOperator != CigarOperator.SOFT_CLIP && lastOperator == CigarOperator.SOFT_CLIP);
    }
}
