package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.gene.RefFlatReader.RefFlatColumns;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;

public final class CollectRnaSeqMetricsTest extends CommandLineProgramTest {
    @Override
    public String getTestedClassName() {
        return CollectRnaSeqMetrics.class.getSimpleName();
    }

    @Test
    public void basic() throws Exception {
        final String sequence = "chr1";
        final String ignoredSequence = "chrM";

        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        // Set seed so that strandedness is consistent among runs.
        builder.setRandomSeed(0);
        final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);
        builder.addPair("pair1", sequenceIndex, 45, 475);
        builder.addPair("pair2", sequenceIndex, 90, 225);
        builder.addPair("pair3", sequenceIndex, 120, 600);
        builder.addFrag("frag1", sequenceIndex, 150, true);
        builder.addFrag("frag2", sequenceIndex, 450, true);
        builder.addFrag("frag3", sequenceIndex, 225, false);
        builder.addPair("rrnaPair", sequenceIndex, 400, 500);

        builder.addFrag("ignoredFrag", builder.getHeader().getSequenceIndex(ignoredSequence), 1, false);

        final File samFile = BaseTest.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        try (final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile)) {
            for (final SAMRecord rec : builder.getRecords()) samWriter.addAlignment(rec);
        }

        // Create an interval list with one ribosomal interval.
        final Interval rRnaInterval = new Interval(sequence, 300, 520, true, "rRNA");
        final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
        rRnaIntervalList.add(rRnaInterval);
        final File rRnaIntervalsFile = BaseTest.createTempFile("tmp.rRna.", ".interval_list");
        rRnaIntervalList.write(rRnaIntervalsFile);

        // Generate the metrics.
        final File metricsFile = BaseTest.createTempFile("tmp.", ".rna_metrics");

        final String[] args = new String[] {
                "--input", samFile.getAbsolutePath(),
                "--output", metricsFile.getAbsolutePath(),
                "--REF_FLAT", getRefFlatFile(sequence).getAbsolutePath(),
                "--RIBOSOMAL_INTERVALS", rRnaIntervalsFile.getAbsolutePath(),
                "--STRAND_SPECIFICITY", "SECOND_READ_TRANSCRIPTION_STRAND",
                "--IGNORE_SEQUENCE", ignoredSequence
        };
        runCommandLine(args);

        final MetricsFile<RnaSeqMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(metricsFile));

        final RnaSeqMetrics metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 396);
        Assert.assertEquals(metrics.PF_BASES, 432);
        Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 108L);
        Assert.assertEquals(metrics.CODING_BASES, 136);
        Assert.assertEquals(metrics.UTR_BASES, 51);
        Assert.assertEquals(metrics.INTRONIC_BASES, 50);
        Assert.assertEquals(metrics.INTERGENIC_BASES, 51);
        Assert.assertEquals(metrics.CORRECT_STRAND_READS, 3);
        Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 4);
        Assert.assertEquals(metrics.IGNORED_READS, 1);
    }

    @Test
    public void testMultiLevel() throws Exception {
        final String sequence = "chr1";
        final String ignoredSequence = "chrM";

        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate, false);
        // Set seed so that strandedness is consistent among runs.
        builder.setRandomSeed(0);
        final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("2");
        rg1.setSample("Sample");
        rg1.setLibrary("foo");
        builder.setReadGroup(rg1);
        builder.addPair("pair1", sequenceIndex, 45, 475);
        builder.addPair("pair2", sequenceIndex, 90, 225);
        builder.addFrag("frag1", sequenceIndex, 150, true);
        builder.addFrag("frag2", sequenceIndex, 450, true);

        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("3");
        rg2.setSample("Sample");
        rg2.setLibrary("bar");
        builder.setReadGroup(rg2);
        builder.addPair("pair3", sequenceIndex, 120, 600);
        builder.addFrag("frag3", sequenceIndex, 225, false);
        builder.addPair("rrnaPair", sequenceIndex, 400, 500);

        builder.addFrag("ignoredFrag", builder.getHeader().getSequenceIndex(ignoredSequence), 1, false);

        final File samFile = BaseTest.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        try (final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile)) {
            for (final SAMRecord rec : builder.getRecords()) samWriter.addAlignment(rec);
        }

        // Create an interval list with one ribosomal interval.
        final Interval rRnaInterval = new Interval(sequence, 300, 520, true, "rRNA");
        final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
        rRnaIntervalList.add(rRnaInterval);
        final File rRnaIntervalsFile = BaseTest.createTempFile("tmp.rRna.", ".interval_list");
        rRnaIntervalList.write(rRnaIntervalsFile);

        // Generate the metrics.
        final File metricsFile = BaseTest.createTempFile("tmp.", ".rna_metrics");

        final String[] args = new String[] {
                "--input", samFile.getAbsolutePath(),
                "--output", metricsFile.getAbsolutePath(),
                "--REF_FLAT", getRefFlatFile(sequence).getAbsolutePath(),
                "--RIBOSOMAL_INTERVALS", rRnaIntervalsFile.getAbsolutePath(),
                "--STRAND_SPECIFICITY", "SECOND_READ_TRANSCRIPTION_STRAND",
                "--IGNORE_SEQUENCE", ignoredSequence,
                "--LEVEL", "SAMPLE",
                "--LEVEL", "LIBRARY"
        };
        runCommandLine(args);

        final MetricsFile<RnaSeqMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(metricsFile));

        for (final RnaSeqMetrics metrics : output.getMetrics()) {
            if (metrics.LIBRARY == null) {
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 396);
                Assert.assertEquals(metrics.PF_BASES, 432);
                Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 108L);
                Assert.assertEquals(metrics.CODING_BASES, 136);
                Assert.assertEquals(metrics.UTR_BASES, 51);
                Assert.assertEquals(metrics.INTRONIC_BASES, 50);
                Assert.assertEquals(metrics.INTERGENIC_BASES, 51);
                Assert.assertEquals(metrics.CORRECT_STRAND_READS, 3);
                Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 4);
                Assert.assertEquals(metrics.IGNORED_READS, 1);
            }
            else if (metrics.LIBRARY.equals("foo")) {
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 216);
                Assert.assertEquals(metrics.PF_BASES, 216);
                Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 36L);
                Assert.assertEquals(metrics.CODING_BASES, 89);
                Assert.assertEquals(metrics.UTR_BASES, 51);
                Assert.assertEquals(metrics.INTRONIC_BASES, 25);
                Assert.assertEquals(metrics.INTERGENIC_BASES, 15);
                Assert.assertEquals(metrics.CORRECT_STRAND_READS, 3);
                Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 2);
                Assert.assertEquals(metrics.IGNORED_READS, 0);

            }
            else if (metrics.LIBRARY.equals("bar")) {
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 180);
                Assert.assertEquals(metrics.PF_BASES, 216);
                Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 72L);
                Assert.assertEquals(metrics.CODING_BASES, 47);
                Assert.assertEquals(metrics.UTR_BASES, 0);
                Assert.assertEquals(metrics.INTRONIC_BASES, 25);
                Assert.assertEquals(metrics.INTERGENIC_BASES, 36);
                Assert.assertEquals(metrics.CORRECT_STRAND_READS, 0);
                Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 2);
                Assert.assertEquals(metrics.IGNORED_READS, 1);

            }
        }
    }

    public File getRefFlatFile(String sequence) throws Exception {
        // Create a refFlat file with a single gene containing two targets, one of which is overlapped by the
        // ribosomal interval.
        final String[] refFlatFields = new String[RefFlatColumns.values().length];
        refFlatFields[RefFlatColumns.GENE_NAME.ordinal()] = "myGene";
        refFlatFields[RefFlatColumns.TRANSCRIPT_NAME.ordinal()] = "myTranscript";
        refFlatFields[RefFlatColumns.CHROMOSOME.ordinal()] = sequence;
        refFlatFields[RefFlatColumns.STRAND.ordinal()] = "+";
        refFlatFields[RefFlatColumns.TX_START.ordinal()] = "49";
        refFlatFields[RefFlatColumns.TX_END.ordinal()] = "500";
        refFlatFields[RefFlatColumns.CDS_START.ordinal()] = "74";
        refFlatFields[RefFlatColumns.CDS_END.ordinal()] = "400";
        refFlatFields[RefFlatColumns.EXON_COUNT.ordinal()] = "2";
        refFlatFields[RefFlatColumns.EXON_STARTS.ordinal()] = "49,249";
        refFlatFields[RefFlatColumns.EXON_ENDS.ordinal()] = "200,500";

        final File refFlatFile = BaseTest.createTempFile("tmp.", ".refFlat");
        try (final PrintStream refFlatStream = new PrintStream(refFlatFile)) {
            refFlatStream.println(StringUtil.join("\t", refFlatFields));
        }

        return refFlatFile;
    }
}
