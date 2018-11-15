package org.broadinstitute.hellbender.tools;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Stream;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "(Experimental) Processes reads from a MITESeq experiment.",
        oneLineSummary = "(EXPERIMENTAL) Processes reads from a MITESeq experiment.",
        programGroup = CoverageAnalysisProgramGroup.class
)
@BetaFeature
public class AnalyzeSaturationMutagenesis extends GATKTool {
    @Argument(doc = "minimum quality score for analyzed portion of read", fullName = "min-q")
    @VisibleForTesting static int minQ = 30;

    @Argument(doc = "minimum size of high-quality portion of read", fullName = "min-length")
    @VisibleForTesting static int minLength = 15;

    @Argument(doc = "minimum number of wt calls flanking variant", fullName = "min-flanking-length")
    private static int minFlankingLength = 18;

    @Argument(doc = "reference indices of the ORF (1-based, closed), for example, '134-180,214-238'", fullName = "orf")
    private static String orfCoords;

    @Argument(doc = "minimum number of observations of reported variants", fullName = "min-variant-obs")
    private static long minVariantObservations = 3;

    @Argument(doc = "codon translation (a string of 64 amino acid codes", fullName = "codon-translation")
    private static String codonTranslation = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF";

    @Argument(doc = "output file prefix", fullName = "output-file-prefix", shortName = "O")
    private static String outputFilePrefix;

    @Argument(doc = "paired mode", fullName = "paired-mode")
    private static boolean pairedMode = true;

    private byte[] refSeq; // the amplicon -- all bytes are converted to upper-case 'A', 'C', 'G', or 'T', no nonsense
    private long[] refCoverage; // number of molecules aligning to each reference position -- same length as above
    private long[] refCoverageSizeHistogram; // number of molecules having a given reference coverage size
    private IntervalCounter intervalCounter; // count of molecules having a particular [start, stop) on reference

    private CodonTracker codonTracker;

    // a map of SNV sets onto number of observations of that set of variations
    private final HopscotchMap<SNVCollectionCount, Long, SNVCollectionCount> variationCounts =
            new HopscotchMap<>(10000000);

    private ReadCounts readCounts;

    // a bunch of mutually exclusive counts of molecules
    private long nWildTypeMolecules = 0; // number of molecules in which no variation from reference was detected
    private long nInconsistentPairs = 0; // number of molecules where mates with conflicting variants in overlap region
    private long nInsufficientFlankMolecules = 0; // number of molecules where variation was too close to end of region
    private long nLowQualityVariantMolecules = 0; // number of molecules where a variation was called with low quality
    private long nCalledVariantMolecules = 0; // number of molecules with a least one variant

    // a place to stash the first read of a pair during pairwise processing of the read stream
    private GATKRead read1 = null;

    private static final int UPPERCASE_MASK = 0xDF; // e.g., 'a' & UPPERCASE_MASK == 'A'
    private final static byte NO_CALL = (byte)'-';

    private static final int N_REGULAR_CODONS = 64;
    @VisibleForTesting static final int FRAME_PRESERVING_INDEL_INDEX = 64;
    @VisibleForTesting static final int FRAME_SHIFTING_INDEL_INDEX = 65;
    private static final int CODON_COUNT_ROW_SIZE = 66;

    private static final String[] labelForCodonValue = {
        "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
        "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
        "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
        "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
    };

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        initializeRefSeq();
        codonTracker = new CodonTracker(orfCoords, refSeq);
        if ( codonTranslation.length() != N_REGULAR_CODONS ) {
            throw new UserException("codon-translation string must contain exactly 64 characters");
        }
    }

    @Override
    public void traverse() {
        // ignore non-primary alignments
        final Stream<GATKRead> reads = getTransformedReadStream(ReadFilterLibrary.PRIMARY_LINE);

        if ( !pairedMode ) {
            reads.forEach(read -> {
                try {
                    applyReport(getReadReport(read, refSeq, readCounts));
                } catch ( final Exception e ) {
                    final String readName = read.getName();
                    throw new GATKException("Caught unexpected exception on read " +
                            readCounts.getNReads() + ": " + readName, e);
                }
            });
        } else {
            read1 = null;
            reads.forEach(read -> {
                try {
                    if ( !read.isPaired() ) {
                        applyReport(getReadReport(read, refSeq, readCounts));
                    } else if ( read1 == null ) {
                        read1 = read;
                    } else if ( !read1.getName().equals(read.getName()) ) {
                        System.err.println("Read " + read1.getName() + " has no mate.");
                        applyReport(getReadReport(read1, refSeq, readCounts));
                        read1 = read;
                    } else {
                        applyReport(combineReports(getReadReport(read1, refSeq, readCounts),
                                                    getReadReport(read, refSeq, readCounts)));
                        read1 = null;
                    }
                } catch ( final Exception e ) {
                    final String readName = read.getName();
                    throw new GATKException("Caught unexpected exception on read " +
                            readCounts.getNReads() + ": " + readName, e);
                }
            });
            if ( read1 != null ) {
                System.err.println("Read " + read1.getName() + " has no mate.");
                try {
                    applyReport(getReadReport(read1, refSeq, readCounts));
                } catch ( final Exception e ) {
                    final String readName = read1.getName();
                    throw new GATKException("Caught unexpected exception on read " +
                            readCounts.getNReads() + ": " + readName, e);
                }
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        writeVariationCounts(getVariationEntries());
        writeRefCoverage();
        writeCodonCounts();
        writeCodonFractions();
        writeAACounts();
        writeAAFractions();
        writeReadCounts();
        writeCoverageSizeHistogram();
        return null;
    }

    private List<SNVCollectionCount> getVariationEntries() {
        final long outputSize =
                variationCounts.stream().filter(entry -> entry.getCount() >= minVariantObservations).count();
        final List<SNVCollectionCount> variationEntries = new ArrayList<>((int)outputSize);
        for ( final SNVCollectionCount entry : variationCounts ) {
            final long count = entry.getCount();
            if ( count >= minVariantObservations ) {
                variationEntries.add(entry);
            }
        }

        // natural order for SNVCollectionCount is lexically on the SNV list.
        variationEntries.sort(Comparator.naturalOrder());
        return variationEntries;
    }

    private void writeVariationCounts( final List<SNVCollectionCount> variationEntries ) {
        final String variantsFile = outputFilePrefix + ".variantCounts";
        try ( final BufferedWriter writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(variantsFile))) ) {
            final DecimalFormat formatter = new DecimalFormat("0.0");
            for ( final SNVCollectionCount entry : variationEntries ) {
                writer.write(Long.toString(entry.getCount()));
                writer.write('\t');
                final List<SNV> snvs = entry.getSNVs();
                final int start = Math.max(0, snvs.get(0).getRefIndex() - minFlankingLength);
                final int end = Math.min(refSeq.length, snvs.get(snvs.size() - 1).getRefIndex() + minFlankingLength);
                writer.write(Long.toString(intervalCounter.countSpanners(start, end)));
                writer.write('\t');
                writer.write(formatter.format(entry.getMeanRefCoverage()));
                writer.write('\t');
                writer.write(Integer.toString(snvs.size()));
                String sep = "\t";
                for ( final SNV snv : snvs ) {
                    writer.write(sep);
                    sep = ", ";
                    writer.write(snv.toString());
                }
                describeVariantsAsCodons(writer, snvs);
                writer.newLine();
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + variantsFile, ioe);
        }
    }

    private void describeVariantsAsCodons( final BufferedWriter writer, final List<SNV> snvs ) throws IOException {
        final List<CodonVariation> codonVariations = codonTracker.encodeSNVsAsCodons(snvs);
        final int[] refCodonValues = codonTracker.getRefCodonValues();

        writer.write('\t');
        writer.write(Long.toString(codonVariations.size()));

        String sep = "\t";
        for ( final CodonVariation variation : codonVariations ) {
            writer.write(sep);
            sep = ", ";
            final int codonId = variation.getCodonId();
            writer.write(Integer.toString(codonId+1));
            writer.write(':');
            if ( variation.isFrameshift() ) {
                writer.write("FS");
            } else {
                writer.write(variation.isInsertion() ? "---" : labelForCodonValue[refCodonValues[codonId]]);
                writer.write('>');
                writer.write(variation.isDeletion() ? "---" : labelForCodonValue[variation.getCodonValue()]);
            }
        }

        sep = "\t";
        for ( final CodonVariation variation : codonVariations ) {
            writer.write(sep);
            sep = ", ";
            final int codonId = variation.getCodonId();
            if ( variation.isFrameshift() ) {
                writer.write("FS");
            } else if ( variation.isInsertion() ) {
                writer.write("I:->");
                writer.write(codonTranslation.charAt(variation.getCodonValue()));
            } else if ( variation.isDeletion() ) {
                writer.write("D:");
                writer.write(codonTranslation.charAt(refCodonValues[codonId]));
                writer.write(":-");
            } else {
                final char fromAA = codonTranslation.charAt(refCodonValues[codonId]);
                final char toAA = codonTranslation.charAt(variation.getCodonValue());
                final char label = fromAA == toAA ? 'S' : CodonTracker.isStop(variation.getCodonValue()) ? 'N' : 'M';
                writer.write(label);
                writer.write(':');
                writer.write(fromAA);
                writer.write('>');
                writer.write(toAA);
            }
        }
    }

    private void writeRefCoverage() {
        final String refCoverageFile = outputFilePrefix + ".refCoverage";
        try ( final BufferedWriter writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(refCoverageFile))) ) {
            writer.write("RefPos\tCoverage");
            writer.newLine();
            int refPos = 1;
            for ( final long coverage : refCoverage ) {
                writer.write(Integer.toString(refPos++));
                writer.write('\t');
                writer.write(Long.toString(coverage));
                writer.newLine();
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + refCoverageFile, ioe);
        }
    }

    private void writeCodonCounts() {
        final String codonCountsFile = outputFilePrefix + ".codonCounts";
        try ( final BufferedWriter writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(codonCountsFile))) ) {
            final long[][] codonCounts = codonTracker.getCodonCounts();
            final int nCodons = codonCounts.length;
            writer.write("AAA\tAAC\tAAG\tAAT\tACA\tACC\tACG\tACT\tAGA\tAGC\tAGG\tAGT\tATA\tATC\tATG\tATT\t" +
                    "CAA\tCAC\tCAG\tCAT\tCCA\tCCC\tCCG\tCCT\tCGA\tCGC\tCGG\tCGT\tCTA\tCTC\tCTG\tCTT\t" +
                    "GAA\tGAC\tGAG\tGAT\tGCA\tGCC\tGCG\tGCT\tGGA\tGGC\tGGG\tGGT\tGTA\tGTC\tGTG\tGTT\t" +
                    "TAA\tTAC\tTAG\tTAT\tTCA\tTCC\tTCG\tTCT\tTGA\tTGC\tTGG\tTGT\tTTA\tTTC\tTTG\tTTT\t" +
                    "NFS\tFS\tTotal");
            writer.newLine();
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                final long[] rowCounts = codonCounts[codonId];
                long total = 0;
                for ( final long count : rowCounts ) {
                    writer.write(Long.toString(count));
                    writer.write('\t');
                    total += count;
                }
                writer.write(Long.toString(total));
                writer.newLine();
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + codonCountsFile, ioe);
        }
    }

    private void writeCodonFractions() {
        final String codonFractFile = outputFilePrefix + ".codonFractions";
        try ( final BufferedWriter writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(codonFractFile))) ) {
            final long[][] codonCounts = codonTracker.getCodonCounts();
            final int nCodons = codonCounts.length;
            writer.write("Codon" +
                    "   AAA   AAC   AAG   AAT   ACA   ACC   ACG   ACT" +
                    "   AGA   AGC   AGG   AGT   ATA   ATC   ATG   ATT" +
                    "   CAA   CAC   CAG   CAT   CCA   CCC   CCG   CCT" +
                    "   CGA   CGC   CGG   CGT   CTA   CTC   CTG   CTT" +
                    "   GAA   GAC   GAG   GAT   GCA   GCC   GCG   GCT" +
                    "   GGA   GGC   GGG   GGT   GTA   GTC   GTG   GTT" +
                    "   TAA   TAC   TAG   TAT   TCA   TCC   TCG   TCT" +
                    "   TGA   TGC   TGG   TGT   TTA   TTC   TTG   TTT" +
                    "   NFS    FS    Total");
            writer.newLine();
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                writer.write(String.format("%5d", codonId + 1));
                final long[] rowCounts = codonCounts[codonId];
                final long total = Arrays.stream(rowCounts).sum();
                for ( final long count : rowCounts ) {
                    writer.write(String.format("%6.2f", 100. * count / total));
                }
                writer.write(String.format("%9d", total));
                writer.newLine();
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + codonFractFile, ioe);
        }
    }

    private void writeAACounts() {
        final String aaCountsFile = outputFilePrefix + ".aaCounts";
        try ( final BufferedWriter writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(aaCountsFile))) ) {
            final long[][] codonCounts = codonTracker.getCodonCounts();
            final int nCodons = codonCounts.length;
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                final long[] rowCounts = codonCounts[codonId];
                final SortedMap<Character, Long> aaCounts = new TreeMap<>();
                for ( int codonValue = 0; codonValue != N_REGULAR_CODONS; ++codonValue ) {
                    aaCounts.merge(codonTranslation.charAt(codonValue), rowCounts[codonValue], Long::sum);
                }
                if ( codonId == 0 ) {
                    String prefix = "";
                    for ( final char chr : aaCounts.keySet() ) {
                        writer.write(prefix);
                        prefix = "\t";
                        writer.write(chr);
                    }
                    writer.newLine();
                }
                String prefix = "";
                for ( final long count : aaCounts.values() ) {
                    writer.write(prefix);
                    prefix = "\t";
                    writer.write(Long.toString(count));
                }
                writer.newLine();
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + aaCountsFile, ioe);
        }
    }

    private void writeAAFractions() {
        final String aaFractFile = outputFilePrefix + ".aaFractions";
        try ( final BufferedWriter writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(aaFractFile))) ) {
            final long[][] codonCounts = codonTracker.getCodonCounts();
            final int nCodons = codonCounts.length;
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                final long[] rowCounts = codonCounts[codonId];
                final SortedMap<Character, Long> aaCounts = new TreeMap<>();
                for ( int codonValue = 0; codonValue != N_REGULAR_CODONS; ++codonValue ) {
                    aaCounts.merge(codonTranslation.charAt(codonValue), rowCounts[codonValue], Long::sum);
                }
                if ( codonId == 0 ) {
                    writer.write("Codon");
                    for ( final char chr : aaCounts.keySet() ) {
                        writer.write("     ");
                        writer.write(chr);
                    }
                    writer.write("    Total");
                    writer.newLine();
                }
                writer.write(String.format("%5d", codonId + 1));
                final long total = Arrays.stream(rowCounts).sum();
                for ( final long count : aaCounts.values() ) {
                    writer.write(String.format("%6.2f", 100. * count / total));
                }
                writer.write(String.format("%9d", total));
                writer.newLine();
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + aaFractFile, ioe);
        }
    }

    private void writeReadCounts() {
        final String readCountsFile = outputFilePrefix + ".readCounts";
        try ( final BufferedWriter writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(readCountsFile))) ) {
            final DecimalFormat df = new DecimalFormat("0.000");
            final long totalReads = readCounts.getNReads();
            writer.write("Total Reads:\t" + totalReads);
            writer.newLine();
            final long nUnmappedReads = readCounts.getNReadsUnmapped();
            writer.write("Unmapped Reads:\t" + nUnmappedReads + "\t" +
                    df.format(100. * nUnmappedReads / totalReads) + "%");
            writer.newLine();
            final long nLowQualityReads = readCounts.getNReadsLowQuality();
            writer.write("LowQ Reads:\t" + nLowQualityReads + "\t" +
                    df.format(100. * nLowQualityReads / totalReads) + "%");
            writer.newLine();
            final long totalMolecules = nInconsistentPairs + nWildTypeMolecules + nInsufficientFlankMolecules +
                    nLowQualityVariantMolecules + nCalledVariantMolecules;
            writer.write("Number of inconsistent pair molecules:\t" + nInconsistentPairs + "\t" +
                    df.format(100. * nInconsistentPairs / totalMolecules) + "%");
            writer.newLine();
            writer.write("Number of wild type molecules:\t" + nWildTypeMolecules + "\t" +
                    df.format(100. * nWildTypeMolecules / totalMolecules) + "%");
            writer.newLine();
            writer.write("Number of insufficient flank molecules:\t" + nInsufficientFlankMolecules + "\t" +
                    df.format(100. * nInsufficientFlankMolecules / totalMolecules) + "%");
            writer.newLine();
            writer.write("Number of low quality variation molecules:\t" + nLowQualityVariantMolecules + "\t" +
                    df.format(100. * nLowQualityVariantMolecules / totalMolecules) + "%");
            writer.newLine();
            writer.write("Number of called variant molecules:\t" + nCalledVariantMolecules + "\t" +
                    df.format(100. * nCalledVariantMolecules / totalMolecules) + "%");
            writer.newLine();
            writer.write("Base calls evaluated for variants:\t" +
                    df.format(100. * Arrays.stream(refCoverage).sum() / readCounts.getNTotalBaseCalls()) + "%");
            writer.newLine();
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + readCountsFile, ioe);
        }
    }

    private void writeCoverageSizeHistogram() {
        if ( refCoverageSizeHistogram == null ) {
            return;
        }

        final String trimCountsFile = outputFilePrefix + ".coverageLengthCounts";
        try ( final BufferedWriter writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(trimCountsFile))) ) {
            int len = refCoverageSizeHistogram.length;
            // peel off the empty tail -- but always print the 1st 100 entries
            while ( len > 101 ) {
                if ( refCoverageSizeHistogram[len - 1] > 0 ) break;
                len -= 1;
            }
            for ( int idx = 1; idx != len; ++idx ) {
                writer.write(Integer.toString(idx));
                writer.write('\t');
                writer.write(Long.toString(refCoverageSizeHistogram[idx]));
                writer.newLine();
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + trimCountsFile, ioe);
        }
    }

    @VisibleForTesting final static class ReadCounts {
        private long nReadsTotal = 0; // total number of reads in input bam
        private long nReadsUnmapped = 0; // number of unmapped reads in bam
        private long nReadsLowQuality = 0; // number of reads where trimming made the read disappear
        private long nTotalBaseCalls = 0; // number of base calls over all reads in bam (mapped or not, call trimmed or not)

        public void bumpNReads() { nReadsTotal += 1; }
        public long getNReads() { return nReadsTotal; }

        public void bumpNReadsUnmapped() { nReadsUnmapped += 1; }
        public long getNReadsUnmapped() { return nReadsUnmapped; }

        public void bumpNLowQualityReads() { nReadsLowQuality += 1; }
        public long getNReadsLowQuality() { return nReadsLowQuality; }

        public void addBaseCalls( final long nBases ) { nTotalBaseCalls += nBases; }
        public long getNTotalBaseCalls() { return nTotalBaseCalls; }
    }

    // describes an interval on some sequence as a pair of offsets (0-based, half-open).
    @VisibleForTesting final static class Interval {
        private final int start;
        private final int end;

        public Interval( final int start, final int end ) {
            this.start = start;
            this.end = end;
        }

        public int getStart() { return start; }
        public int getEnd() { return end; }
        public int size() { return end - start; }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof Interval && equals((Interval)obj);
        }

        public boolean equals( final Interval that ) {
            return this.start == that.start && this.end == that.end;
        }

        @Override
        public int hashCode() {
            return 47 * (47 * start + end);
        }

        public static Interval nullInterval = new Interval(0, 0);
    }

    // a description of a single-base deviation from reference
    @VisibleForTesting static final class SNV implements Comparable<SNV> {
        private final int refIndex;
        private final byte refCall;
        private final byte variantCall;
        private final byte qual; // not a part of equality or hash

        public SNV( final int refIndex, final byte refCall, final byte variantCall, final byte qual ) {
            this.refIndex = refIndex;
            this.refCall = refCall;
            this.variantCall = variantCall;
            this.qual = qual;
        }

        public int getRefIndex() { return refIndex; }
        public byte getRefCall() { return refCall; }
        public byte getVariantCall() { return variantCall; }
        public byte getQuality() { return qual; }

        @Override
        public int hashCode() {
            return 47 * (47 * (47 * refIndex + refCall) + variantCall);
        }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof SNV && equals((SNV)obj);
        }

        public boolean equals( final SNV that ) {
            return this.refIndex == that.refIndex &&
                    this.refCall == that.refCall &&
                    this.variantCall == that.variantCall;
        }

        @Override
        public int compareTo( final SNV that ) {
            int result = Integer.compare(this.refIndex, that.refIndex);
            if ( result == 0 ) result = Byte.compare(this.refCall, that.refCall);
            if ( result == 0 ) result = Byte.compare(this.variantCall, that.variantCall);
            return result;
        }

        @Override
        public String toString() {
            return (refIndex + 1) + ":" + (char)refCall + ">" + (char)variantCall;
        }
    }

    // a count of molecules that start and stop at particular places on the reference
    // this allows us to calculate the number of molecules that span any given interval
    @VisibleForTesting static final class IntervalCounter {
        // triangular matrix indexed first by starting position on reference, and second by the size of the interval
        final long[][] counts;

        public IntervalCounter( final int refLen ) {
            counts = new long[refLen][];
            for ( int rowIndex = 0; rowIndex != refLen; ++rowIndex ) {
                counts[rowIndex] = new long[refLen - rowIndex + 1];
            }
        }

        public void addCount( final Interval refCoverage ) {
            final int refStart = refCoverage.getStart();
            counts[refStart][refCoverage.getEnd() - refStart] += 1;
        }

        public long countSpanners( final int refStart, final int refEnd ) {
            long total = 0;
            for ( int rowIndex = 0; rowIndex <= refStart; ++rowIndex ) {
                final long[] row = counts[rowIndex];
                for ( int spanIndex = refEnd - rowIndex; spanIndex < row.length; ++spanIndex ) {
                    total += row[spanIndex];
                }
            }
            return total;
        }
    }

    // an array of SNVs that serves as a key (comparison, equality, and hashCode depend only on this part)
    // and a count of the number of observations of those SNVs, plus the total reference coverage over all observations
    // implements a Map.Entry as a single object to conserve memory
    @VisibleForTesting static final class SNVCollectionCount
            implements Map.Entry<SNVCollectionCount, Long>, Comparable<SNVCollectionCount> {
        private static final SNV[] emptyArray = new SNV[0];
        private final SNV[] snvs;
        private long count; // number of observations of this set of SNVs -- not part of key
        private int totalRefCoverage; // the sum of the reference coverage over all observations -- not part of key
        private final int hash; // depends only on the array of SNVs

        public SNVCollectionCount( final List<SNV> snvs, final int refCoverage ) {
            this.snvs = snvs.toArray(emptyArray);
            this.count = 1;
            this.totalRefCoverage = refCoverage;
            int hashVal = 0;
            for ( final SNV snv : snvs ) {
                hashVal = 47 * hashVal + snv.hashCode();
            }
            hash = 47 * hashVal;
        }

        @Override
        public SNVCollectionCount getKey() { return this; }

        @Override
        public Long getValue() { return count; }

        @Override
        public Long setValue( final Long value ) {
            final Long result = count;
            count = value;
            return result;
        }

        public List<SNV> getSNVs() { return Arrays.asList(snvs); }
        public long getCount() { return count; }

        public void bumpCount( final int refCoverage ) {
            count += 1;
            totalRefCoverage += refCoverage;
        }

        public float getMeanRefCoverage() {
            return 1.f * totalRefCoverage / count;
        }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof SNVCollectionCount && equals((SNVCollectionCount)obj);
        }

        public boolean equals( final SNVCollectionCount that ) {
            return this.hash == that.hash && Arrays.equals(this.snvs, that.snvs);
        }

        @Override
        public int hashCode() {
            return hash;
        }

        @Override
        public int compareTo( final SNVCollectionCount that ) {
            final int minSize = Math.min(this.snvs.length, that.snvs.length);
            for ( int idx = 0; idx != minSize; ++idx ) {
                final int result = this.snvs[idx].compareTo(that.snvs[idx]);
                if ( result != 0 ) {
                    return result;
                }
            }
            return Integer.compare(this.snvs.length, that.snvs.length);
        }
    }

    @VisibleForTesting enum CodonVariationType {
        FRAMESHIFT,
        INSERTION,
        DELETION,
        MODIFICATION
    }

    @VisibleForTesting static final class CodonVariation {
        private final int codonId;
        private final int codonValue; // ignored for FRAMESHIFT and DELETION
        private final CodonVariationType variationType;

        private CodonVariation( final int codonId, final int codonValue, final CodonVariationType variationType ) {
            this.codonId = codonId;
            this.codonValue = codonValue;
            this.variationType = variationType;
        }

        public int getCodonId() { return codonId; }
        public int getCodonValue() { return codonValue; }
        public CodonVariationType getVariationType() { return variationType; }
        public boolean isFrameshift() { return variationType == CodonVariationType.FRAMESHIFT; }
        public boolean isInsertion() { return variationType == CodonVariationType.INSERTION; }
        public boolean isDeletion() { return variationType == CodonVariationType.DELETION; }
        public boolean isModification() { return variationType == CodonVariationType.MODIFICATION; }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof CodonVariation && equals((CodonVariation)obj);
        }

        public boolean equals( final CodonVariation that ) {
            return this.codonId == that.codonId &&
                    this.codonValue == that.codonValue &&
                    this.variationType == that.variationType;
        }

        @Override
        public int hashCode() {
            return 47 * (47 * (47 * codonId + codonValue) + variationType.ordinal());
        }

        public static CodonVariation createFrameshift( final int codonId ) {
            return new CodonVariation(codonId, -1, CodonVariationType.FRAMESHIFT);
        }
        public static CodonVariation createInsertion( final int codonId, final int codonValue ) {
            return new CodonVariation(codonId, codonValue, CodonVariationType.INSERTION);
        }
        public static CodonVariation createDeletion( final int codonId ) {
            return new CodonVariation(codonId, -1, CodonVariationType.DELETION);
        }
        public static CodonVariation createModification( final int codonId, final int codonValue ) {
            return new CodonVariation(codonId, codonValue, CodonVariationType.MODIFICATION);
        }
    }

    @VisibleForTesting static final class CodonTracker {
        private final byte[] refSeq;
        private final List<Interval> exonList;
        private final long[][] codonCounts;
        private final int[] refCodonValues;
        @VisibleForTesting static int NO_FRAME_SHIFT_CODON = -1;

        public CodonTracker( final String orfCoords, final byte[] refSeq ) {
            this.refSeq = refSeq;
            exonList = getExons(orfCoords, refSeq.length);

            codonCounts = new long[exonList.stream().mapToInt(Interval::size).sum() / 3][];
            for ( int codonId = 0; codonId != codonCounts.length; ++codonId ) {
                codonCounts[codonId] = new long[CODON_COUNT_ROW_SIZE];
            }

            refCodonValues = parseReferenceIntoCodons(refSeq, exonList);
        }

        public int[] getRefCodonValues() { return refCodonValues; }
        public long[][] getCodonCounts() { return codonCounts; }

        public List<CodonVariation> encodeSNVsAsCodons( final List<SNV> snvs ) {
            final List<CodonVariation> codonVariations = new ArrayList<>();
            final Iterator<SNV> snvIterator = snvs.iterator();
            SNV snv = null;
            while ( snvIterator.hasNext() ) {
                final SNV testSNV = snvIterator.next();
                if ( isExonic(testSNV.getRefIndex()) ) {
                    snv = testSNV;
                    break;
                }
            }

            int frameShiftCodonId = findFrameShift(snvs);

            final int lastExonEnd = exonList.get(exonList.size() - 1).getEnd();
            while ( snv != null ) {
                int refIndex = snv.getRefIndex();
                if ( refIndex >= lastExonEnd ) {
                    break;
                }
                Interval curExon = null;
                final Iterator<Interval> exonIterator = exonList.iterator();
                while ( exonIterator.hasNext() ) {
                    final Interval testExon = exonIterator.next();
                    if ( testExon.getStart() <= refIndex && testExon.getEnd() > refIndex ) {
                        curExon = testExon;
                        break;
                    }
                }
                if ( curExon == null ) {
                    throw new GATKException("can't find current exon, even though refIndex should be exonic.");
                }

                int codonId = exonicBaseIndex(refIndex);
                int codonPhase = codonId % 3;
                codonId /= 3;

                int codonValue = refCodonValues[codonId];
                if ( codonPhase == 0 ) codonValue = 0;
                else if ( codonPhase == 1 ) codonValue >>= 4;
                else codonValue >>= 2;

                int leadLag = 0;
                do {
                    boolean codonValueAltered = false;
                    boolean bumpRefIndex = false;
                    if ( snv == null || snv.getRefIndex() != refIndex ) {
                        codonValue = (codonValue << 2) | "ACGT".indexOf(refSeq[refIndex]);
                        codonValueAltered = true;
                        bumpRefIndex = true;
                    } else {
                        if ( snv.getVariantCall() == NO_CALL ) {
                            if ( codonId == frameShiftCodonId ) {
                                codonVariations.add(CodonVariation.createFrameshift(codonId));
                                frameShiftCodonId = NO_FRAME_SHIFT_CODON;
                            }
                            if ( --leadLag == -3 ) {
                                codonVariations.add(CodonVariation.createDeletion(codonId));
                                if ( ++codonId == refCodonValues.length ) {
                                    return codonVariations;
                                }
                                leadLag = 0;
                            }
                            bumpRefIndex = true;
                        } else if ( snv.getRefCall() == NO_CALL ) {
                            leadLag += 1;
                            codonValue = (codonValue << 2) | "ACGT".indexOf(snv.getVariantCall());
                            codonValueAltered = true;
                        } else {
                            codonValue = (codonValue << 2) | "ACGT".indexOf(snv.getVariantCall());
                            codonValueAltered = true;
                            bumpRefIndex = true;
                        }
                        snv = null;
                        while ( snvIterator.hasNext() ) {
                            final SNV testSNV = snvIterator.next();
                            if ( testSNV.getRefIndex() >= lastExonEnd || isExonic(testSNV.getRefIndex()) ) {
                                snv = testSNV;
                                break;
                            }
                        }
                    }
                    if ( bumpRefIndex ) {
                        if ( ++refIndex == curExon.getEnd() ) {
                            if ( exonIterator.hasNext() ) {
                                curExon = exonIterator.next();
                                refIndex = curExon.getStart();
                            }
                        }
                        if ( refIndex == refSeq.length ) {
                            return codonVariations;
                        }
                    }
                    if ( codonValueAltered ) {
                        if ( ++codonPhase == 3 ) {
                            if ( codonId == frameShiftCodonId ) {
                                codonVariations.add(CodonVariation.createFrameshift(codonId));
                                frameShiftCodonId = NO_FRAME_SHIFT_CODON;
                            }
                            if ( leadLag == 3 ) {
                                codonVariations.add(CodonVariation.createInsertion(codonId, codonValue));
                                leadLag = 0;
                                codonId -= 1;
                            } else if ( codonValue != refCodonValues[codonId] ) {
                                codonVariations.add(CodonVariation.createModification(codonId, codonValue));
                            }
                            if ( isStop(codonValue) ) {
                                return codonVariations;
                            }
                            if ( ++codonId == refCodonValues.length ) {
                                return codonVariations;
                            }
                            codonPhase = 0;
                            codonValue = 0;
                        }
                    }
                } while ( leadLag != 0 || codonPhase != 0 );
            }

            return codonVariations;
        }

        public void reportVariantCodonCounts( final Interval refCoverage,
                                              final List<CodonVariation> variantCodons ) {
            final int startingCodonId = (exonicBaseIndex(refCoverage.getStart()) + 2) / 3;
            final int endingCodonId = exonicBaseIndex(refCoverage.getEnd()) / 3;
            final Iterator<CodonVariation> variantCodonIterator = variantCodons.iterator();
            CodonVariation codonVariation = variantCodonIterator.hasNext() ? variantCodonIterator.next() : null;
            for ( int codonId = startingCodonId; codonId < endingCodonId; ++codonId ) {
                while ( codonVariation != null && codonVariation.getCodonId() < codonId ) {
                    codonVariation = variantCodonIterator.hasNext() ? variantCodonIterator.next() : null;
                }
                if ( codonVariation == null || codonVariation.getCodonId() != codonId ) {
                    codonCounts[codonId][refCodonValues[codonId]] += 1;
                } else {
                    boolean framePreservingIndel = false;
                    do {
                        switch ( codonVariation.getVariationType() ) {
                            case FRAMESHIFT:
                                codonCounts[codonId][FRAME_SHIFTING_INDEL_INDEX] += 1;
                                break;
                            case DELETION:
                            case INSERTION:
                                framePreservingIndel = true;
                                break;
                            case MODIFICATION:
                                codonCounts[codonId][codonVariation.getCodonValue()] += 1;
                                break;
                        }
                        codonVariation = variantCodonIterator.hasNext() ? variantCodonIterator.next() : null;
                    } while ( codonVariation != null && codonVariation.getCodonId() == codonId );
                    if ( framePreservingIndel ) {
                        codonCounts[codonId][FRAME_PRESERVING_INDEL_INDEX] += 1;
                    }
                }
            }
        }

        public void reportWildCodonCounts( final Interval refCoverage ) {
            final int startingCodonId = (exonicBaseIndex(refCoverage.getStart()) + 2) / 3;
            final int endingCodonId = exonicBaseIndex(refCoverage.getEnd()) / 3;
            for ( int codonId = startingCodonId; codonId < endingCodonId; ++codonId ) {
                codonCounts[codonId][refCodonValues[codonId]] += 1;
            }
        }

        @VisibleForTesting static List<Interval> getExons( final String orfCoords, final int refLen ) {
            final List<Interval> exonList = new ArrayList<>();
            for ( final String coordPair : orfCoords.split(",") ) {
                final String[] coords = coordPair.split("-");
                if ( coords.length != 2 ) {
                    throw new UserException("Can't interpret ORF as list of pairs of coords: " + orfCoords);
                }
                try {
                    final int start = Integer.valueOf(coords[0]);
                    if ( start < 1 ) {
                        throw new UserException("Coordinates of ORF are 1-based.");
                    }
                    final int end = Integer.valueOf(coords[1]);
                    if ( end < start ) {
                        throw new UserException("Found ORF end coordinate less than start: " + orfCoords);
                    }
                    if ( end > refLen ) {
                        throw new UserException("Found ORF end coordinate larger than reference length: " + orfCoords);
                    }
                    // convert 1-based, inclusive intervals to 0-based, half-open
                    final Interval exon = new Interval(start - 1, end);
                    exonList.add(exon);
                } catch ( final NumberFormatException nfe ) {
                    throw new UserException("Can't interpret ORF coords as integers: " + orfCoords);
                }
                for ( int idx = 1; idx < exonList.size(); ++idx ) {
                    if ( exonList.get(idx - 1).getEnd() >= exonList.get(idx).getStart() ) {
                        throw new UserException("ORF coordinates are not sorted: " + orfCoords);
                    }
                }
            }

            final int orfLen = exonList.stream().mapToInt(Interval::size).sum();
            if ( (orfLen % 3) != 0 ) {
                throw new UserException("ORF length must be divisible by 3.");
            }

            return exonList;
        }

        @VisibleForTesting static int[] parseReferenceIntoCodons( final byte[] refSeq, final List<Interval> exonList ) {
            final int nCodons = exonList.stream().mapToInt(Interval::size).sum() / 3;
            final int[] refCodonValues = new int[nCodons];
            int codonId = 0;
            int codonPhase = 0;
            int codonValue = 0;
            for ( final Interval exon : exonList ) {
                final int exonEnd = exon.getEnd();
                for ( int refIndex = exon.getStart(); refIndex != exonEnd; ++refIndex ) {
                    codonValue = (codonValue << 2) | "ACGT".indexOf(refSeq[refIndex]);
                    if ( ++codonPhase == 3 ) {
                        if ( isStop(codonValue) && codonId != nCodons - 1 ) {
                            final int idx = refIndex + 1;
                            throw new UserException("There is an upstream stop codon at reference index " + idx + ".");
                        }
                        refCodonValues[codonId] = codonValue;
                        codonValue = 0;
                        codonPhase = 0;
                        codonId += 1;
                    }
                }
            }

            final int START_CODON = 0x0E;
            if ( refCodonValues[0] != START_CODON ) {
                System.err.println("WARNING:  Your ORF does not start with the expected ATG codon.");
            }
            final int lastCodon = refCodonValues[nCodons - 1];
            if ( !isStop(lastCodon) ) {
                System.err.println("WARNING:  Your ORF does not end with the expected stop codon.");
            }

            return refCodonValues;
        }

        @VisibleForTesting int findFrameShift( final List<SNV> snvs ) {
            int codonId = NO_FRAME_SHIFT_CODON;
            int leadLag = 0;
            for ( final SNV snv : snvs ) {
                if ( isExonic(snv.getRefIndex()) ) {
                    if ( snv.getVariantCall() == NO_CALL ) {
                        if ( leadLag == 0 ) {
                            codonId = exonicBaseIndex(snv.getRefIndex()) / 3;
                        }
                        if ( --leadLag == -3 ) {
                            leadLag = 0;
                        }
                    } else if ( snv.getRefCall() == NO_CALL ) {
                        if ( leadLag == 0 ) {
                            codonId = exonicBaseIndex(snv.getRefIndex()) / 3;
                        }
                        if ( ++leadLag == 3 ) {
                            leadLag = 0;
                        }
                    }
                    if ( leadLag == 0 ) {
                        codonId = NO_FRAME_SHIFT_CODON;
                    }
                }
            }
            return codonId;
        }

        static boolean isStop( final int codonValue ) {
            final int STOP_OCH = 0x30;
            final int STOP_AMB = 0x32;
            final int STOP_OPA = 0x38;
            return codonValue == STOP_OCH || codonValue == STOP_AMB || codonValue == STOP_OPA;
        }

        @VisibleForTesting boolean isExonic( final int refIndex ) {
            for ( final Interval exon : exonList ) {
                if ( exon.getStart() > refIndex ) return false;
                if ( exon.getEnd() > refIndex ) return true;
            }
            return false;
        }

        @VisibleForTesting int exonicBaseIndex( final int refIndex ) {
            int baseCount = 0;
            for ( final Interval exon : exonList ) {
                if ( refIndex >= exon.getEnd() ) {
                    baseCount += exon.size();
                } else {
                    if ( refIndex > exon.getStart() ) {
                        baseCount += refIndex - exon.getStart();
                    }
                    break;
                }
            }
            return baseCount;
        }
    }

    @VisibleForTesting interface ReadReport {
        List<Interval> getRefCoverage();
        List<SNV> getVariations();
        boolean hasCleanFlanks( final int minFlankingLength, final int refLength );
        void apply( final CodonTracker codonTracker,
                    final HopscotchMap<SNVCollectionCount, Long, SNVCollectionCount> variationCounts,
                    final IntervalCounter intervalCounter,
                    final int coverage );

        static void updateVariationCounts(
                final List<SNV> variations,
                final HopscotchMap<SNVCollectionCount, Long, SNVCollectionCount> variationCounts,
                final int coverage ) {
            final SNVCollectionCount newVal = new SNVCollectionCount(variations, coverage);
            final SNVCollectionCount oldVal = variationCounts.find(newVal);
            if ( oldVal != null ) {
                oldVal.bumpCount(coverage);
            } else {
                variationCounts.add(newVal);
            }
        }
    }

    @VisibleForTesting final static class SingleReadReport implements ReadReport {
        final List<Interval> refCoverage;
        final List<SNV> snvList;

        public SingleReadReport( final List<Interval> refCoverage,
                           final List<SNV> snvList ) {
            this.refCoverage = refCoverage;
            this.snvList = snvList;
        }

        @Override
        public List<Interval> getRefCoverage() {
            return refCoverage;
        }

        public int getFirstRefIndex() { return refCoverage.get(0).getStart(); }
        public int getLastRefIndex() { return refCoverage.get(refCoverage.size() - 1).getEnd(); }

        @Override
        public List<SNV> getVariations() {
            return snvList;
        }

        @Override
        public boolean hasCleanFlanks( final int minFlankingLength, final int refLength ) {
            return hasCleanLeftFlank(minFlankingLength) && hasCleanRightFlank(minFlankingLength, refLength);
        }

        public boolean hasCleanLeftFlank( final int minFlankingLength ) {
            return snvList.isEmpty() ||
                    Math.max(0, snvList.get(0).getRefIndex() - minFlankingLength) >= getFirstRefIndex();
        }

        public boolean hasCleanRightFlank( final int minFlankingLength, final int refLength ) {
            return snvList.isEmpty() ||
                    Math.min(refLength - 1, snvList.get(snvList.size() - 1).getRefIndex() + minFlankingLength) <
                            getLastRefIndex();
        }

        @Override
        public void apply( final CodonTracker codonTracker,
                           final HopscotchMap<SNVCollectionCount, Long, SNVCollectionCount> variationCounts,
                           final IntervalCounter intervalCounter,
                           final int coverage ) {
            final Interval totalCoverage = new Interval(getFirstRefIndex(), getLastRefIndex());
            intervalCounter.addCount(totalCoverage);

            if ( snvList.isEmpty() ) {
                codonTracker.reportWildCodonCounts(totalCoverage);
            } else {
                codonTracker.reportVariantCodonCounts(totalCoverage, codonTracker.encodeSNVsAsCodons(snvList));
                ReadReport.updateVariationCounts(snvList, variationCounts, coverage);
            }
        }

        public static SingleReadReport nullReport = new SingleReadReport(new ArrayList<>(), new ArrayList<>());
    }

    @VisibleForTesting final static class PairedReadReport implements ReadReport {
        final SingleReadReport report1;
        final SingleReadReport report2;
        final List<SNV> combinedSNVList;

        public PairedReadReport( final SingleReadReport report1, final SingleReadReport report2 ) {
            this.report1 = report1;
            this.report2 = report2;
            this.combinedSNVList = combineVariations(report1, report2);
        }

        @Override
        public List<Interval> getRefCoverage() {
            final List<Interval> refCoverage1 = report1.getRefCoverage();
            final List<Interval> refCoverage2 = report2.getRefCoverage();
            if ( refCoverage1.isEmpty() ) return refCoverage2;
            if ( refCoverage2.isEmpty() ) return refCoverage1;

            final List<Interval> combinedCoverage = new ArrayList<>(refCoverage1.size() + refCoverage2.size());
            final Iterator<Interval> refCoverageItr1 = refCoverage1.iterator();
            final Iterator<Interval> refCoverageItr2 = refCoverage2.iterator();
            Interval refInterval1 = refCoverageItr1.next();
            Interval refInterval2 = refCoverageItr2.next();
            Interval curInterval;
            if ( refInterval1.getStart() < refInterval2.getStart() ) {
                curInterval = refInterval1;
                refInterval1 = refCoverageItr1.hasNext() ? refCoverageItr1.next() : null;
            } else {
                curInterval = refInterval2;
                refInterval2 = refCoverageItr2.hasNext() ? refCoverageItr2.next() : null;
            }
            while ( refInterval1 != null || refInterval2 != null ) {
                final Interval testInterval;
                if ( refInterval1 == null ) {
                    testInterval = refInterval2;
                    refInterval2 = refCoverageItr2.hasNext() ? refCoverageItr2.next() : null;
                } else if ( refInterval2 == null ) {
                    testInterval = refInterval1;
                    refInterval1 = refCoverageItr1.hasNext() ? refCoverageItr1.next() : null;
                } else if ( refInterval1.getStart() < refInterval2.getStart() ) {
                    testInterval = refInterval1;
                    refInterval1 = refCoverageItr1.hasNext() ? refCoverageItr1.next() : null;
                } else {
                    testInterval = refInterval2;
                    refInterval2 = refCoverageItr2.hasNext() ? refCoverageItr2.next() : null;
                }
                if ( curInterval.getEnd() < testInterval.getStart() ) {
                    combinedCoverage.add(curInterval);
                    curInterval = testInterval;
                } else {
                    curInterval =
                            new Interval(curInterval.getStart(), Math.max(curInterval.getEnd(), testInterval.getEnd()));
                }
            }
            combinedCoverage.add(curInterval);

            return combinedCoverage;
        }

        @Override
        public List<SNV> getVariations() { return combinedSNVList; }

        @Override
        public boolean hasCleanFlanks( final int minFlankingLength, final int refLength ) {
            // return true if the earlier-starting report is clean, and the later-ending report is clean
            return (report1.getFirstRefIndex() < report2.getFirstRefIndex() ? report1 : report2)
                    .hasCleanLeftFlank(minFlankingLength)
                    &&
                   (report1.getLastRefIndex() > report2.getLastRefIndex() ? report1 : report2)
                    .hasCleanRightFlank(minFlankingLength, refLength);
        }

        @Override
        public void apply( final CodonTracker codonTracker,
                           final HopscotchMap<SNVCollectionCount, Long, SNVCollectionCount> variationCounts,
                           final IntervalCounter intervalCounter,
                           final int coverage ) {
            final int start1 = report1.getFirstRefIndex();
            final int end1 = report1.getLastRefIndex();
            final int start2 = report2.getFirstRefIndex();
            final int end2 = report2.getLastRefIndex();

            final int overlapStart = Math.max(start1, start2);
            final int overlapEnd = Math.min(end1, end2);
            if ( overlapStart > overlapEnd ) { // disjoint -- report each independently
                report1.apply(codonTracker, variationCounts, intervalCounter, coverage);
                report2.apply(codonTracker, variationCounts, intervalCounter, coverage);
            } else {
                final Interval totalCoverage = new Interval(Math.min(start1, start2), Math.max(end1, end2));
                intervalCounter.addCount(totalCoverage);

                if ( combinedSNVList.isEmpty() ) {
                    codonTracker.reportWildCodonCounts(totalCoverage);
                } else {
                    final List<CodonVariation> codonVariations = codonTracker.encodeSNVsAsCodons(combinedSNVList);
                    codonTracker.reportVariantCodonCounts(totalCoverage, codonVariations);
                    ReadReport.updateVariationCounts(combinedSNVList, variationCounts, coverage);
                }
            }
        }

        private static List<SNV> combineVariations( final SingleReadReport report1, final SingleReadReport report2 ) {
            if ( report1.getVariations().isEmpty() ) return report2.getVariations();
            if ( report2.getVariations().isEmpty() ) return report1.getVariations();

            final int overlapStart = Math.max(report1.getFirstRefIndex(), report2.getFirstRefIndex());
            final int overlapEnd = Math.min(report1.getLastRefIndex(), report2.getLastRefIndex());

            final List<SNV> combinedSNVs = new ArrayList<>();
            final Iterator<SNV> iterator1 = report1.getVariations().iterator();
            final Iterator<SNV> iterator2 = report2.getVariations().iterator();
            SNV snv1 = iterator1.hasNext() ? iterator1.next() : null;
            SNV snv2 = iterator2.hasNext() ? iterator2.next() : null;
            while ( snv1 != null || snv2 != null ) {
                final SNV next;
                if ( snv1 == null ) {
                    next = snv2;
                    snv2 = iterator2.hasNext() ? iterator2.next() : null;
                    final int refIndex = next.getRefIndex();
                    if ( refIndex >= overlapStart && refIndex < overlapEnd ) {
                        return null;
                    }
                } else if ( snv2 == null ) {
                    next = snv1;
                    snv1 = iterator1.hasNext() ? iterator1.next() : null;
                    final int refIndex = next.getRefIndex();
                    if ( refIndex >= overlapStart && refIndex < overlapEnd ) {
                        return null;
                    }
                } else {
                    final int refIndex1 = snv1.getRefIndex();
                    final int refIndex2 = snv2.getRefIndex();
                    if ( refIndex1 < refIndex2 ) {
                        next = snv1;
                        snv1 = iterator1.hasNext() ? iterator1.next() : null;
                        if ( refIndex1 >= overlapStart && refIndex1 < overlapEnd ) {
                            return null;
                        }
                    } else if ( refIndex2 < refIndex1 ) {
                        next = snv2;
                        snv2 = iterator2.hasNext() ? iterator2.next() : null;
                        if ( refIndex2 >= overlapStart && refIndex2 < overlapEnd ) {
                            return null;
                        }
                    } else if ( !snv1.equals(snv2) ) {
                        return null;
                    } else {
                        next = snv1.getQuality() > snv2.getQuality() ? snv1 : snv2;
                        snv1 = iterator1.hasNext() ? iterator1.next() : null;
                        snv2 = iterator2.hasNext() ? iterator2.next() : null;
                    }
                }
                combinedSNVs.add(next);
            }
            return combinedSNVs;
        }
    }

    private void initializeRefSeq() {
        final ReferenceDataSource reference = ReferenceDataSource.of(referenceArguments.getReferencePath());
        final SAMSequenceDictionary seqDict = reference.getSequenceDictionary();
        if ( seqDict.size() != 1 ) {
            throw new UserException("Expecting a reference with a single contig. " +
                    "The supplied reference has " + seqDict.size() + " contigs.");
        }
        final SAMSequenceRecord tig0 = seqDict.getSequence(0);
        final int refSeqLen = tig0.getSequenceLength();
        final SimpleInterval wholeTig = new SimpleInterval(tig0.getSequenceName(), 1, refSeqLen);
        refSeq = Arrays.copyOf(reference.queryAndPrefetch(wholeTig).getBases(), refSeqLen);
        for ( int idx = 0; idx < refSeqLen; ++idx ) {
            switch ( refSeq[idx] &= UPPERCASE_MASK ) { // make into upper case
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                    break;
                default:
                    throw new UserException("Reference sequence contains something other than A, C, G, and T.");
            }
        }
        refCoverage = new long[refSeq.length];
        refCoverageSizeHistogram = new long[refSeq.length + 1];
        intervalCounter = new IntervalCounter(refSeq.length);
    }

    @VisibleForTesting static SingleReadReport getReadReport( final GATKRead read,
                                                              final byte[] refSeq,
                                                              final ReadCounts readCounts ) {
        readCounts.bumpNReads();
        readCounts.addBaseCalls(read.getLength());

        if ( read.isUnmapped() ) {
            readCounts.bumpNReadsUnmapped();
            return SingleReadReport.nullReport;
        }

        final Interval trim = calculateTrim(read.getBaseQualitiesNoCopy(), readCounts);
        if ( trim.size() == 0 ) {
            return SingleReadReport.nullReport;
        }

        final int readStart = trim.getStart();
        final int readEnd = trim.getEnd();
        final Cigar cigar = read.getCigar();
        final Iterator<CigarElement> cigarIterator = cigar.getCigarElements().iterator();
        CigarElement cigarElement = cigarIterator.next();
        CigarOperator cigarOperator = cigarElement.getOperator();
        int cigarElementCount = cigarElement.getLength();

        final byte[] readSeq = read.getBasesNoCopy();
        final byte[] readQuals = read.getBaseQualitiesNoCopy();

        final List<SNV> variations = new ArrayList<>();

        int refIndex = read.getStart() - 1; // 0-based numbering
        int readIndex = 0;

        // pretend that soft-clips are matches
        if ( cigarOperator == CigarOperator.S ) {
            refIndex -= cigarElementCount;
        }

        final List<Interval> refCoverageList = new ArrayList<>();
        int refCoverageBegin = -1;
        int refCoverageEnd = -1;

        while ( true ) {
            if ( readIndex >= readStart && refIndex >= 0 ) {
                if ( refCoverageBegin == -1 ) {
                    refCoverageBegin = refIndex;
                    refCoverageEnd = refIndex;
                }
                if ( cigarOperator == CigarOperator.D ) {
                    variations.add(new SNV(refIndex, refSeq[refIndex], NO_CALL, readQuals[readIndex]));
                } else if ( cigarOperator == CigarOperator.I ) {
                    final byte call = (byte)(readSeq[readIndex] & UPPERCASE_MASK);
                    variations.add(new SNV(refIndex, NO_CALL, call, readQuals[readIndex]));
                } else if ( cigarOperator == CigarOperator.M || cigarOperator == CigarOperator.S ) {
                    final byte call = (byte)(readSeq[readIndex] & UPPERCASE_MASK);
                    if ( call != refSeq[refIndex] ) {
                        variations.add(new SNV(refIndex, refSeq[refIndex], call, readQuals[readIndex]));
                    }
                    if ( refIndex == refCoverageEnd ) {
                        refCoverageEnd += 1;
                    } else {
                        refCoverageList.add(new Interval(refCoverageBegin, refCoverageEnd));
                        refCoverageBegin = refIndex;
                        refCoverageEnd = refIndex + 1;
                    }
                } else {
                    throw new GATKException("unanticipated cigar operator: " + cigarOperator.toString());
                }
            }

            if ( cigarOperator != CigarOperator.D ) {
                if ( ++readIndex == readEnd ) {
                    break;
                }
            }

            if ( cigarOperator != CigarOperator.I ) {
                if ( ++refIndex == refSeq.length )
                    break;
            }

            if ( --cigarElementCount == 0 ) {
                if ( !cigarIterator.hasNext() ) {
                    throw new GATKException("unexpectedly exhausted cigar iterator");
                }
                cigarElement = cigarIterator.next();
                cigarOperator = cigarElement.getOperator();
                cigarElementCount = cigarElement.getLength();
            }
        }

        if ( refCoverageBegin < refCoverageEnd ) {
            refCoverageList.add(new Interval(refCoverageBegin, refCoverageEnd));
        }

        return new SingleReadReport(refCoverageList, variations);
    }

    @VisibleForTesting static Interval calculateTrim( final byte[] quals, final ReadCounts readCounts ) {
        // find initial end-trim
        int readStart = 0;
        int hiQCount = 0;
        while ( readStart < quals.length ) {
            if ( quals[readStart] < minQ ) {
                hiQCount = 0;
            } else if ( ++hiQCount == minLength ) {
                break;
            }
            readStart += 1;
        }
        if ( readStart == quals.length ) {
            readCounts.bumpNLowQualityReads();
            return Interval.nullInterval;
        }
        readStart -= minLength - 1;

        // find final end-trim
        int readEnd = quals.length - 1;
        hiQCount = 0;
        while ( readEnd >= 0 ) {
            if ( quals[readEnd] < minQ ) {
                hiQCount = 0;
            } else if ( ++hiQCount == minLength ) {
                break;
            }
            readEnd -= 1;
        }
        readEnd += minLength;

        return new Interval(readStart, readEnd);
    }

    private static ReadReport combineReports( final SingleReadReport report1, final SingleReadReport report2 ) {
        if ( report1.getRefCoverage().isEmpty() ) return report2;
        if ( report2.getRefCoverage().isEmpty() ) return report1;
        return new PairedReadReport(report1, report2);
    }

    private void applyReport( final ReadReport readReport ) {
        final List<Interval> refCoverageList = readReport.getRefCoverage();
        if ( refCoverageList.isEmpty() ) {
            return;
        }

        final List<SNV> variations = readReport.getVariations();
        if ( variations == null ) {
            nInconsistentPairs += 1;
        } else if ( variations.isEmpty() ) {
            nWildTypeMolecules += 1;
            final int coverage = updateRefCoverage(refCoverageList);
            readReport.apply(codonTracker, variationCounts, intervalCounter, coverage);
        } else if ( variations.stream().anyMatch(
                snv -> snv.getQuality() < minQ || "-ACGT".indexOf(snv.getVariantCall()) == -1) ) {
            nLowQualityVariantMolecules += 1;
        } else if ( !readReport.hasCleanFlanks(minFlankingLength, refSeq.length) ) {
            nInsufficientFlankMolecules += 1;
        } else {
            nCalledVariantMolecules += 1;
            final int coverage = updateRefCoverage(refCoverageList);
            readReport.apply(codonTracker, variationCounts, intervalCounter, coverage);
        }
    }

    private int updateRefCoverage( final List<Interval> refCoverageList ) {
        int coverage = 0;
        for ( final Interval refInterval : refCoverageList ) {
            final int refIntervalStart = refInterval.getStart();
            final int refIntervalEnd = refInterval.getEnd();
            coverage += refIntervalEnd - refIntervalStart;
            for ( int idx = refInterval.getStart(); idx != refIntervalEnd; ++idx ) {
                refCoverage[idx] += 1;
            }
        }
        refCoverageSizeHistogram[coverage] += 1;
        return coverage;
    }
}
