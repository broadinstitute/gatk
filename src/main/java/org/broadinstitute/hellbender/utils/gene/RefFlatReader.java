package org.broadinstitute.hellbender.utils.gene;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.utils.text.parsers.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.util.*;

/**
 * Loads gene annotations from a refFlat file into an OverlapDetector<Gene>. Discards annotations that are not
 * internally consistent, e.g. transcripts on different chromosomes or different strands.
 */
public final class RefFlatReader {
    private static final Logger LOG = LogManager.getLogger();
    // These are in the order that columns appear in refFlat format.
    public enum RefFlatColumns{GENE_NAME, TRANSCRIPT_NAME, CHROMOSOME, STRAND, TX_START, TX_END, CDS_START, CDS_END,
        EXON_COUNT, EXON_STARTS, EXON_ENDS}

    private static final String[] RefFlatColumnLabels = new String[RefFlatColumns.values().length];

    static {
        for (int i = 0; i < RefFlatColumnLabels.length; ++i) {
            RefFlatColumnLabels[i] = RefFlatColumns.values()[i].name();
        }
    }

    private final File refFlatFile;
    private final SAMSequenceDictionary sequenceDictionary;

    RefFlatReader(final File refFlatFile, final SAMSequenceDictionary sequenceDictionary) {
        this.refFlatFile = refFlatFile;
        this.sequenceDictionary = sequenceDictionary;
    }

    static OverlapDetector<Gene> load(final File refFlatFile, final SAMSequenceDictionary sequenceDictionary) {
        return new RefFlatReader(refFlatFile, sequenceDictionary).load();
    }

    OverlapDetector<Gene> load() {
        final OverlapDetector<Gene> overlapDetector = new OverlapDetector<>(0, 0);

        final int expectedColumns = RefFlatColumns.values().length;
        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(refFlatFile, RefFlatColumnLabels);
        final Map<String, List<TabbedTextFileWithHeaderParser.Row>> refFlatLinesByGene =
                new HashMap<>();

        for (final TabbedTextFileWithHeaderParser.Row row : parser) {
            final int lineNumber = parser.getCurrentLineNumber(); // getCurrentLineNumber returns the number of the next line
            if (row.getFields().length != expectedColumns) {
                throw new GeneAnnotationException("Wrong number of fields in refFlat file " + refFlatFile + " at line " + lineNumber);
            }
            final String geneName = row.getField(RefFlatColumns.GENE_NAME.name());
            final String transcriptName = row.getField(RefFlatColumns.TRANSCRIPT_NAME.name());
            final String transcriptDescription = geneName + ":" + transcriptName;
            final String chromosome = row.getField(RefFlatColumns.CHROMOSOME.name());
            if (!isSequenceRecognized(chromosome)) {
                LOG.debug("Skipping " + transcriptDescription + " due to unrecognized sequence " + chromosome);
            } else {
                List<TabbedTextFileWithHeaderParser.Row> transcriptLines = refFlatLinesByGene.get(geneName);
                if (transcriptLines == null) {
                    transcriptLines = new ArrayList<>();
                    refFlatLinesByGene.put(geneName, transcriptLines);
                }
                transcriptLines.add(row);
            }
        }

        int longestInterval = 0;
        int numIntervalsOver1MB = 0;

        for (final List<TabbedTextFileWithHeaderParser.Row> transcriptLines : refFlatLinesByGene.values()) {
            try {
                final Gene gene = makeGeneFromRefFlatLines(transcriptLines);
                overlapDetector.addLhs(gene, gene);
                if (gene.length() > longestInterval) longestInterval = gene.length();
                if (gene.length() > 1000000) ++numIntervalsOver1MB;
            } catch (Exception e) {
                LOG.debug(e.getMessage() + " -- skipping");
            }
        }
        LOG.debug("Longest gene: " + longestInterval + "; number of genes > 1MB: " + numIntervalsOver1MB);
        return overlapDetector;
    }

    private boolean isSequenceRecognized(final String sequence) {
        return (sequenceDictionary.getSequence(sequence) != null);
    }


    private Gene makeGeneFromRefFlatLines(final List<TabbedTextFileWithHeaderParser.Row> transcriptLines) {
        final String geneName = transcriptLines.get(0).getField(RefFlatColumns.GENE_NAME.name());
        final String strandStr = transcriptLines.get(0).getField(RefFlatColumns.STRAND.name());
        final boolean negative = strandStr.equals("-");
        final String chromosome = transcriptLines.get(0).getField(RefFlatColumns.CHROMOSOME.name());

        // Figure out the extend of the gene
        int start = Integer.MAX_VALUE;
        int end = Integer.MIN_VALUE;
        for (final TabbedTextFileWithHeaderParser.Row row: transcriptLines) {
            start = Math.min(start, row.getIntegerField(RefFlatColumns.TX_START.name()) + 1);
            end   = Math.max(end, row.getIntegerField(RefFlatColumns.TX_END.name()));
        }

        final Gene gene = new Gene(chromosome, start, end, negative, geneName);

        for (final TabbedTextFileWithHeaderParser.Row row: transcriptLines) {
            if (!strandStr.equals(row.getField(RefFlatColumns.STRAND.name()))) {
                throw new GeneAnnotationException("Strand disagreement in refFlat file for gene " + geneName);
            }
            if (!chromosome.equals(row.getField(RefFlatColumns.CHROMOSOME.name()))) {
                throw new GeneAnnotationException("Chromosome disagreement(" + chromosome + " != " +
                        row.getField(RefFlatColumns.CHROMOSOME.name()) + ") in refFlat file for gene " + geneName);
            }

            // This adds it to the Gene also
            final Gene.Transcript tx = makeTranscriptFromRefFlatLine(gene, row);
        }

        return gene;
    }

    /**
     * Conversion from 0-based half-open to 1-based inclusive intervals is done here.
     */
    private Gene.Transcript makeTranscriptFromRefFlatLine(final Gene gene, final TabbedTextFileWithHeaderParser.Row row) {
        final String geneName = row.getField(RefFlatColumns.GENE_NAME.name());
        final String transcriptName = row.getField(RefFlatColumns.TRANSCRIPT_NAME.name());
        final String transcriptDescription = geneName + ":" + transcriptName;
        final int exonCount = Integer.parseInt(row.getField(RefFlatColumns.EXON_COUNT.name()));
        final String[] exonStarts = row.getField(RefFlatColumns.EXON_STARTS.name()).split(",");
        final String[] exonEnds = row.getField(RefFlatColumns.EXON_ENDS.name()).split(",");

        if (exonCount != exonStarts.length) {
            throw new GeneAnnotationException("Number of target starts does not agree with number of targets for " + transcriptDescription);
        }
        if (exonCount != exonEnds.length) {
            throw new GeneAnnotationException("Number of target ends does not agree with number of targets for " + transcriptDescription);
        }

        final int transcriptionStart = row.getIntegerField(RefFlatColumns.TX_START.name()) + 1;
        final int transcriptionEnd = row.getIntegerField(RefFlatColumns.TX_END.name());
        final int codingStart = row.getIntegerField(RefFlatColumns.CDS_START.name()) + 1;
        final int codingEnd = row.getIntegerField(RefFlatColumns.CDS_END.name());

        final Gene.Transcript tx = gene.addTranscript(transcriptName, transcriptionStart, transcriptionEnd, codingStart, codingEnd, exonCount);

        for (int i = 0; i < exonCount; ++i) {
            final Gene.Transcript.Exon e = tx.addExon(Integer.parseInt(exonStarts[i]) + 1, Integer.parseInt(exonEnds[i]));

            if (e.start >= e.end) {
                throw new GeneAnnotationException("Exon has 0 or negative extent for " + transcriptDescription);
            }
            if (i > 0 && tx.exons[i-1].end >= tx.exons[i].start) {
                throw new GeneAnnotationException("Exons overlap for " + transcriptDescription);
            }
        }

        return tx;
    }
}
