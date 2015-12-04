package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;
import java.util.*;

/**
 * Charts quality score distribution within a BAM file.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = "Program to chart " +
                "quality score distributions in a SAM/BAM/CRAM file.",
        oneLineSummary = "Produces metrics for quality score distributions for a SAM/BAM/CRAM file",
        programGroup = QCProgramGroup.class
)
public final class QualityScoreDistribution extends SinglePassSamProgram {
    public static final String R_SCRIPT = "qualityScoreDistribution.R";

    @Argument(shortName="CHART", doc="A file (with .pdf extension) to write the chart to.")
    public File CHART_OUTPUT;

    @Argument(doc="If set to true calculate mean quality over aligned reads only.")
    public boolean ALIGNED_READS_ONLY = false;

    @Argument(shortName="PF", doc="If set to true calculate mean quality over PF reads only.")
    public boolean PF_READS_ONLY = false;

    @Argument(doc="If set to true, include quality for no-call bases in the distribution.")
    public boolean INCLUDE_NO_CALLS = false;

    @Argument(doc = "Should an output plot be created")
    public boolean PRODUCE_PLOT = false;

    private final long[] qCounts  = new long[128];
    private final long[] oqCounts = new long[128];

    /**
     * A subtitle for the plot, usually corresponding to a library.
     */
    private String plotSubtitle = "";

    private final Logger log = LogManager.getLogger(QualityScoreDistribution.class);

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(CHART_OUTPUT);

        // If we're working with a single library, assign that library's name
        // as a suffix to the plot title
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (readGroups.size() == 1) {
            this.plotSubtitle = readGroups.get(0).getLibrary();
            if (null == this.plotSubtitle) this.plotSubtitle = "";
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        // Skip unwanted records
        if (PF_READS_ONLY && rec.getReadFailsVendorQualityCheckFlag()) return;
        if (ALIGNED_READS_ONLY && rec.getReadUnmappedFlag()) return;
        if (rec.isSecondaryOrSupplementary()) return;

        final byte[] bases = rec.getReadBases();
        final byte[] quals = rec.getBaseQualities();
        final byte[] oq    = rec.getOriginalBaseQualities();

        final int length = quals.length;

        for (int i=0; i<length; ++i) {
            if (INCLUDE_NO_CALLS || !SequenceUtil.isNoCall(bases[i])) {
                qCounts[quals[i]]++;
                if (oq != null) oqCounts[oq[i]]++;
            }
        }
    }

    @Override
    protected void finish() {
        // Built the Histograms out of the long[]s
        final Histogram<Byte> qHisto  = new Histogram<>("QUALITY", "COUNT_OF_Q");
        final Histogram<Byte> oqHisto = new Histogram<>("QUALITY", "COUNT_OF_OQ");

        for (int i=0; i< qCounts.length; ++i) {
            if (qCounts[i]  > 0) qHisto.increment( (byte) i, (double) qCounts[i]);
            if (oqCounts[i] > 0) oqHisto.increment((byte) i, (double) oqCounts[i]);
        }

        final MetricsFile<?,Byte> metrics = getMetricsFile();
        metrics.addHistogram(qHisto);
        if (!oqHisto.isEmpty()) metrics.addHistogram(oqHisto);
        metrics.write(OUTPUT);
        if (qHisto.isEmpty() && oqHisto.isEmpty()) {
            log.warn("No valid bases found in input file. No plot will be produced.");
        }
        else if(PRODUCE_PLOT){
            // Now run R to generate a chart
            final RScriptExecutor executor = new RScriptExecutor();
            executor.addScript(new Resource(R_SCRIPT, QualityScoreDistribution.class));
            executor.addArgs(OUTPUT.getAbsolutePath(), CHART_OUTPUT.getAbsolutePath(), INPUT.getName(), plotSubtitle);
            executor.exec();
        }
    }
}
