package org.broadinstitute.hellbender.tools;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
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

/**
 * <p>Process reads from a saturation mutagenesis experiment.</p>
 * <p>This tool processes reads from an experiment that systematically perturbs a mini-gene to ascertain which
 * amino-acid variations are tolerable at each codon of the open reading frame.  It's main job is to discover
 * variations from wild-type sequence among the reads, and to summarize the variations observed.<p>
 *
 * <h3>Input</h3>
 * <ul>
 * <li>A BAM file. Any combination of paired or unpaired reads can be presented.  In paired mode, the default mode
 * of operation of the program where overlapping pairs are combined bedore variant calling, the BAM must be
 * in unsorted or query-name sorted order (so that read pairs are adjacent).  The BAM file should be
 * aligned to the amplicon under scrutiny.  You may also specify a text file listing several BAMs.</li>
 * <li>A fasta-formatted reference file with the wild-type sequence of the gene as a single contig.  It's best
 * to include the entire amplicon including any expected 5' and 3' UTRs.</li>
 * <li>A description of the location of the open reading frame (ORF) within the reference.  For example,
 * <b>--orf 128-1285</b>
 * would describe an ORF that begins at reference position 128 (1-based coordinates) and ends at
 * reference position 1285 (inclusive).
 * It's possible (though unlikely) to have an ORF with multiple exons:  List the exon ranges separated by commas.
 * For example,
 * <b>--orf 128-1284,1384-1433</b>
 * would describe an ORF with two exons.
 * The total length of the ORF must be divisible by 3, and should begin with a start codon and end with
 * a stop codon.</li>
 * <li>An outputPathAndPrefix that specifies where to write the output reports.</li>
 * </ul>
 * <h3>Output</h3>
 * <ul>
 * <li>The most important report is a tab-delimited text file named outputPathAndPrefix.variantCounts which
 * describes the observed variations from reference, the number of times each was observed, the effect of the
 * observed variations on the codons, and their translation into amino acids.</li>
 * <li>There are a number of additional reports that summarize this information for each codon.</li>
 * </ul>
 * <h3>Usage example</h3>
 * <pre>
 *     gatk AnalyzeSaturationMutagenesis \
 *     -I input_reads.bam \
 *     -R referenceGene.fasta \
 *     --orf 128-1285 \
 *     -O /path/to/output/and/reportNamePrefix
 * </pre>
 */

@DocumentedFeature
@CommandLineProgramProperties(
        summary =
"(Experimental) Processes reads from a MITESeq or other saturation mutagenesis experiment.\n" +
"Main output is a tab-delimited text file reportNamePrefix.variantCounts.\n" +
"Columns are:\n" +
"1 - Number of times the described variant was observed\n" +
"2 - Number of molecules that covered the variant region and its flanks\n" +
"3 - Mean length of trimmed, aligned reference coverage in the observed variant molecules\n" +
"4 - The number of base calls that varied\n" +
"5 - The variant base calls as a comma-delimited string, e.g., 17:A>T says reference base A was called as T at reference position 17\n" +
"    A missing base call is represented as a hyphen.\n" +
"6 - The number of variant codons\n" +
"7 - The variant codons as a comma-delimited string, e.g., 3:AAG>CGA says the 3rd codon went from AAG to CGA\n" +
"8 - The variant amino-acids that result from the variant codons, e.g., M:K>R indicates a missense variation from Lysine to Arginine\n" +
"All reference coordinates are 1-based.",
        oneLineSummary = "(EXPERIMENTAL) Processes reads from a MITESeq or other saturation mutagenesis experiment.",
        usageExample = "gatk AnalyzeSaturationMutagenesis -I input_reads.bam -R referenceGene.fasta --orf 128-1285 -O /path/to/output/and/reportNamePrefix",
        programGroup = CoverageAnalysisProgramGroup.class
)
@BetaFeature
public final class AnalyzeSaturationMutagenesis extends GATKTool {
    @Argument(doc = "minimum quality score for analyzed portion of read", fullName = "min-q")
    @VisibleForTesting static int minQ = 30;

    @Argument(doc = "minimum size of high-quality portion of read", fullName = "min-length")
    @VisibleForTesting static int minLength = 15;

    @Argument(doc = "minimum number of wt calls flanking variant", fullName = "min-flanking-length")
    @VisibleForTesting static int minFlankingLength = 18;

    @Argument(doc = "reference interval(s) of the ORF (1-based, inclusive), for example, '134-180,214-238' (no spaces)",
            fullName = "orf")
    private static String orfCoords;

    @Argument(doc = "minimum number of observations of reported variants", fullName = "min-variant-obs")
    private static long minVariantObservations = 3;

    @Argument(doc = "codon translation (a string of 64 amino acid codes", fullName = "codon-translation")
    private static String codonTranslation = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF";

    @Argument(doc = "output file prefix", fullName = "output-file-prefix",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME)
    private static String outputFilePrefix;

    @Argument(doc = "paired mode evaluation of variants (combine mates, when possible)", fullName = "paired-mode")
    private static boolean pairedMode = true;

    private static Reference reference;
    private static CodonTracker codonTracker; // device for turning SNVs into CodonVariations

    // a map of SNV sets onto number of observations of that set of variations
    private final HopscotchMap<SNVCollectionCount, Long, SNVCollectionCount> variationCounts =
            new HopscotchMap<>(10000000);

    private final ReadCounts readCounts = new ReadCounts();
    private final MoleculeCounts moleculeCounts = new MoleculeCounts();

    // a place to stash the first read of a pair during pairwise processing of the read stream
    private GATKRead read1 = null;

    private static final int UPPERCASE_MASK = 0xDF; // e.g., 'a' & UPPERCASE_MASK == 'A'
    private final static byte NO_CALL = (byte)'-';

    private static final String[] LABEL_FOR_CODON_VALUE = {
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
        final SAMFileHeader header = getHeaderForReads();
        if ( pairedMode && header != null && header.getSortOrder() == SortOrder.coordinate ) {
            throw new UserException("In paired mode the BAM cannot be coordinate sorted.  Mates must be adjacent.");
        }
        if ( codonTranslation.length() != CodonTracker.N_REGULAR_CODONS ) {
            throw new UserException("codon-translation string must contain exactly 64 characters");
        }
        reference = new Reference(ReferenceDataSource.of(referenceArguments.getReferencePath()));
        codonTracker = new CodonTracker(orfCoords, reference.getRefSeq(), logger);
    }

    @Override
    public void traverse() {
        final byte[] refSeq = reference.getRefSeq();

        // ignore non-primary alignments
        final Stream<GATKRead> reads = getTransformedReadStream(ReadFilterLibrary.PRIMARY_LINE);

        if ( !pairedMode ) {
            reads.forEach(read -> {
                try {
                    final ReadReport report = getReadReport(read, refSeq, readCounts, reference);
                    report.updateCounts(moleculeCounts, codonTracker, variationCounts, reference);
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
                        final ReadReport report = getReadReport(read, refSeq, readCounts, reference);
                        report.updateCounts(moleculeCounts, codonTracker, variationCounts, reference);
                    } else if ( read1 == null ) {
                        read1 = read;
                    } else if ( !read1.getName().equals(read.getName()) ) {
                        logger.warn("Read " + read1.getName() + " has no mate.");
                        final ReadReport report = getReadReport(read1, refSeq, readCounts, reference);
                        report.updateCounts(moleculeCounts, codonTracker, variationCounts, reference);
                        read1 = read;
                    } else {
                        final ReadReport report = combineReports(
                                getReadReport(read1, refSeq, readCounts, reference),
                                getReadReport(read, refSeq, readCounts, reference));
                        report.updateCounts(moleculeCounts, codonTracker, variationCounts, reference);
                        read1 = null;
                    }
                } catch ( final Exception e ) {
                    final String readName = read.getName();
                    throw new GATKException("Caught unexpected exception on read " +
                            readCounts.getNReads() + ": " + readName, e);
                }
            });
            if ( read1 != null ) {
                logger.warn("Read " + read1.getName() + " has no mate.");
                try {
                    final ReadReport report = getReadReport(read1, refSeq, readCounts, reference);
                    report.updateCounts(moleculeCounts, codonTracker, variationCounts, reference);
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
        return super.onTraversalSuccess();
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
        final int refSeqLength = reference.getRefSeqLength();
        final String variantsFile = outputFilePrefix + ".variantCounts";
        try ( final BufferedWriter writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(variantsFile))) ) {
            final DecimalFormat formatter = new DecimalFormat("0.0");
            for ( final SNVCollectionCount entry : variationEntries ) {
                writer.write(Long.toString(entry.getCount()));
                writer.write('\t');
                final List<SNV> snvs = entry.getSNVs();
                final int start = Math.max(0, snvs.get(0).getRefIndex() - minFlankingLength);
                final int end = Math.min(refSeqLength, snvs.get(snvs.size() - 1).getRefIndex() + minFlankingLength);
                writer.write(Long.toString(reference.countSpanners(start, end)));
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
                writer.write(variation.isInsertion() ? "---" : LABEL_FOR_CODON_VALUE[refCodonValues[codonId]]);
                writer.write('>');
                writer.write(variation.isDeletion() ? "---" : LABEL_FOR_CODON_VALUE[variation.getCodonValue()]);
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
            for ( final long coverage : reference.getCoverage() ) {
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
            for ( final String codonLabel : LABEL_FOR_CODON_VALUE ) {
                writer.write(codonLabel);
                writer.write('\t');
            }
            writer.write("NFS\tFS\tTotal");
            writer.newLine();
            final int nCodons = codonCounts.length;
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
            writer.write("Codon");
            for ( final String codonLabel : LABEL_FOR_CODON_VALUE ) {
                writer.write("   ");
                writer.write(codonLabel);
            }
            writer.write("   NFS    FS    Total");
            writer.newLine();
            final int nCodons = codonCounts.length;
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
                for ( int codonValue = 0; codonValue != CodonTracker.N_REGULAR_CODONS; ++codonValue ) {
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
                for ( int codonValue = 0; codonValue != CodonTracker.N_REGULAR_CODONS; ++codonValue ) {
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
            final long nInconsistentPairs = moleculeCounts.getNInconsistentPairs();
            final long nInsufficientFlankMolecules = moleculeCounts.getnInsufficientFlank();
            final long nWildTypeMolecules = moleculeCounts.getNWildType();
            final long nLowQualityVariantMolecules = moleculeCounts.getNLowQualityVariant();
            final long nCalledVariantMolecules = moleculeCounts.getCalledVariant();
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
                    df.format(100. * reference.getTotalCoverage() / readCounts.getNTotalBaseCalls()) + "%");
            writer.newLine();
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write " + readCountsFile, ioe);
        }
    }

    private void writeCoverageSizeHistogram() {
        final long[] refCoverageSizeHistogram = reference.getCoverageSizeHistogram();

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

    @VisibleForTesting final static class Reference {
        // the amplicon -- all bytes are converted to upper-case 'A', 'C', 'G', or 'T', no nonsense
        private final byte[] refSeq;
        // number of molecules aligning to each reference position -- same length as above
        private final long[] coverage;
        // number of molecules having a given reference coverage size
        private final long[] coverageSizeHistogram;
        // count of molecules having a particular [start, stop) on reference
        private final IntervalCounter intervalCounter;

        public Reference( final ReferenceDataSource refSource ) {
            final SAMSequenceDictionary seqDict = refSource.getSequenceDictionary();
            if ( seqDict.size() != 1 ) {
                throw new UserException("Expecting a reference with a single contig. " +
                        "The supplied reference has " + seqDict.size() + " contigs.");
            }
            final SAMSequenceRecord tig0 = seqDict.getSequence(0);
            final int refSeqLen = tig0.getSequenceLength();
            final SimpleInterval wholeTig = new SimpleInterval(tig0.getSequenceName(), 1, refSeqLen);
            refSeq = Arrays.copyOf(refSource.queryAndPrefetch(wholeTig).getBases(), refSeqLen);
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
            coverage = new long[refSeq.length];
            coverageSizeHistogram = new long[refSeq.length + 1];
            intervalCounter = new IntervalCounter(refSeq.length);
        }

        @VisibleForTesting Reference( final byte[] refSeq ) {
            this.refSeq = refSeq;
            coverage = new long[refSeq.length];
            coverageSizeHistogram = new long[refSeq.length + 1];
            intervalCounter = new IntervalCounter(refSeq.length);
        }

        public int getRefSeqLength() { return refSeq.length; }
        public byte[] getRefSeq() { return refSeq; }
        public long getTotalCoverage() { return Arrays.stream(coverage).sum(); }
        public long[] getCoverage() { return coverage; }
        public long[] getCoverageSizeHistogram() { return coverageSizeHistogram; }
        public long countSpanners( final int refStart, final int refEnd ) {
            return intervalCounter.countSpanners(refStart, refEnd);
        }
        public void updateSpan( final Interval refSpan ) {
            intervalCounter.addCount(refSpan);
        }

        // returns sum of interval sizes
        private int updateCoverage( final List<Interval> refCoverageList ) {
            int coverageLen = 0;
            for ( final Interval refInterval : refCoverageList ) {
                final int refIntervalStart = refInterval.getStart();
                final int refIntervalEnd = refInterval.getEnd();
                coverageLen += refIntervalEnd - refIntervalStart;
                for ( int idx = refInterval.getStart(); idx != refIntervalEnd; ++idx ) {
                    coverage[idx] += 1;
                }
            }
            coverageSizeHistogram[coverageLen] += 1;
            return coverageLen;
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

    // a bunch of mutually exclusive counts of molecules
    @VisibleForTesting final static class MoleculeCounts {
        private long nWildType = 0; // number of molecules in which no variation from reference was detected
        private long nInconsistentPairs = 0; // number of molecules where mates with conflicting variants in overlap region
        private long nInsufficientFlank = 0; // number of molecules where variation was too close to end of region
        private long nLowQualityVariant = 0; // number of molecules where a variation was called with low quality
        private long nCalledVariant = 0; // number of molecules with a least one variant

        public void bumpNWildType() { nWildType += 1; }
        public long getNWildType() { return nWildType; }
        public void bumpInconsistentPairs() { nInconsistentPairs += 1; }
        public long getNInconsistentPairs() { return nInconsistentPairs; }
        public void bumpInsufficientFlank() { nInsufficientFlank += 1; }
        public long getnInsufficientFlank() { return nInsufficientFlank; }
        public void bumpNLowQualityVariant() { nLowQualityVariant += 1; }
        public long getNLowQualityVariant() { return nLowQualityVariant; }
        public void bumpCalledVariant() { nCalledVariant += 1; }
        public long getCalledVariant() { return nCalledVariant; }
    }

    // describes an interval on some sequence as a pair of offsets (0-based, half-open).
    @VisibleForTesting final static class Interval {
        private final int start;
        private final int end;

        public Interval( final int start, final int end ) {
            if ( start < 0 || end < start ) {
                throw new GATKException("Illegal interval: [" + start + "," + end + ")");
            }
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

        public static Interval NULL_INTERVAL = new Interval(0, 0);
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
            final int refEnd = refCoverage.getEnd();
            if ( refEnd > counts.length ) {
                throw new GATKException("illegal span: [" + refStart + "," + refEnd + ")");
            }
            counts[refStart][refCoverage.getEnd() - refStart] += 1;
        }

        public long countSpanners( final int refStart, final int refEnd ) {
            if ( refStart < 0 || refEnd < refStart || refEnd > counts.length ) {
                throw new GATKException("illegal span: [" + refStart + "," + refEnd + ")");
            }
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

    // A class to translate a sequence of base-call variants into a sequence of codon variants.
    // also has some counters to track the number of observations of different codon values at each codon
    @VisibleForTesting static final class CodonTracker {
        private final byte[] refSeq; // wild-type base calls for the reference
        private final List<Interval> exonList; // the ORF
        private final long[][] codonCounts; // for each codon, the number of times we've observed each codon value
        private final int[] refCodonValues; // wild-type codon value for each codon (parsed from refSeq and ORF)

        @VisibleForTesting static int NO_FRAME_SHIFT_CODON = -1;
        @VisibleForTesting static final int N_REGULAR_CODONS = 64;
        @VisibleForTesting static final int FRAME_PRESERVING_INDEL_INDEX = 64;
        @VisibleForTesting static final int FRAME_SHIFTING_INDEL_INDEX = 65;
        private static final int CODON_COUNT_ROW_SIZE = 66;



        public CodonTracker( final String orfCoords, final byte[] refSeq, final Logger logger ) {
            this.refSeq = refSeq;
            exonList = getExons(orfCoords, refSeq.length);

            codonCounts = new long[exonList.stream().mapToInt(Interval::size).sum() / 3][];
            for ( int codonId = 0; codonId != codonCounts.length; ++codonId ) {
                codonCounts[codonId] = new long[CODON_COUNT_ROW_SIZE];
            }

            refCodonValues = parseReferenceIntoCodons(refSeq, exonList, logger);
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

        @VisibleForTesting static int[] parseReferenceIntoCodons( final byte[] refSeq,
                                                                  final List<Interval> exonList,
                                                                  final Logger logger ) {
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
               logger.warn("WARNING:  Your ORF does not start with the expected ATG codon.");
            }
            final int lastCodon = refCodonValues[nCodons - 1];
            if ( !isStop(lastCodon) ) {
                logger.warn("WARNING:  Your ORF does not end with the expected stop codon.");
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

        // is the given reference index a part of the ORF?
        @VisibleForTesting boolean isExonic( final int refIndex ) {
            for ( final Interval exon : exonList ) {
                if ( exon.getStart() > refIndex ) return false;
                if ( exon.getEnd() > refIndex ) return true;
            }
            return false;
        }

        // the codonId for a given reference position.
        // if the position is intronic, it will be the id of the next codon.
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

    @VisibleForTesting static class ReadReport {
        final List<Interval> refCoverage;
        final List<SNV> snvList;

        public ReadReport( final GATKRead read, final Interval trim, final byte[] refSeq ) {
            snvList = new ArrayList<>();
            refCoverage = new ArrayList<>();

            // figure out what parts of the reference are covered by the read, and which base calls are variant
            // with respect to reference.  process only the high-quality part of the read, as directed by the trim
            // interval.  (note that the trim interval refers to an interval of the read, not the reference.)
            final int readStart = trim.getStart();
            final int readEnd = trim.getEnd();
            final Cigar cigar = read.getCigar();
            final Iterator<CigarElement> cigarIterator = cigar.getCigarElements().iterator();
            CigarElement cigarElement = cigarIterator.next();
            CigarOperator cigarOperator = cigarElement.getOperator();
            int cigarElementCount = cigarElement.getLength();

            final byte[] readSeq = read.getBasesNoCopy();
            final byte[] readQuals = read.getBaseQualitiesNoCopy();

            int refIndex = read.getStart() - 1; // 0-based numbering
            int readIndex = 0;

            // pretend that soft-clips are matches:  when the read starts just before the reference the first few bases
            // of reference sequence might be clipped.  we want to capture any variation in those first few bases.
            if ( cigarOperator == CigarOperator.S ) {
                refIndex -= cigarElementCount; // this can go negative and that's ok.
            }

            int refCoverageBegin = -1;
            int refCoverageEnd = -1;

            while ( true ) {
                if ( readIndex >= readStart && refIndex >= 0 ) {
                    if ( refCoverageBegin == -1 ) {
                        refCoverageBegin = refIndex;
                        refCoverageEnd = refIndex;
                    }
                    if ( cigarOperator == CigarOperator.D ) {
                        snvList.add(new SNV(refIndex, refSeq[refIndex], NO_CALL, readQuals[readIndex]));
                    } else if ( cigarOperator == CigarOperator.I ) {
                        final byte call = (byte)(readSeq[readIndex] & UPPERCASE_MASK);
                        snvList.add(new SNV(refIndex, NO_CALL, call, readQuals[readIndex]));
                    } else if ( cigarOperator == CigarOperator.M || cigarOperator == CigarOperator.S ) {
                        final byte call = (byte)(readSeq[readIndex] & UPPERCASE_MASK);
                        if ( call != refSeq[refIndex] ) {
                            snvList.add(new SNV(refIndex, refSeq[refIndex], call, readQuals[readIndex]));
                        }
                        if ( refIndex == refCoverageEnd ) {
                            refCoverageEnd += 1;
                        } else {
                            refCoverage.add(new Interval(refCoverageBegin, refCoverageEnd));
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
                refCoverage.add(new Interval(refCoverageBegin, refCoverageEnd));
            }
        }

        public List<Interval> getRefCoverage() {
            return refCoverage;
        }

        public int getFirstRefIndex() { return refCoverage.get(0).getStart(); }
        public int getLastRefIndex() { return refCoverage.get(refCoverage.size() - 1).getEnd(); }

        public List<SNV> getVariations() {
            return snvList;
        }

        public boolean hasCleanFlanks( final int minFlankingLength, final int refLength ) {
            return hasCleanLeftFlank(minFlankingLength) && hasCleanRightFlank(minFlankingLength, refLength);
        }

        public boolean hasCleanLeftFlank( final int minFlankingLength ) {
            return snvList.isEmpty() ||
                    Math.max(0, snvList.get(0).getRefIndex() - minFlankingLength) >= getFirstRefIndex();
        }

        public boolean hasCleanRightFlank( final int minFlankingLength, final int refLength ) {
            return snvList.isEmpty() ||
                    Math.min(refLength - 1, snvList.get(snvList.size() - 1).getRefIndex() + minFlankingLength)
                            < getLastRefIndex();
        }

        // apply the report: update reference coverage, translate SNVs to codon values, and count 'em all up
        public void updateCounts( final MoleculeCounts moleculeCounts,
                                  final CodonTracker codonTracker,
                                  final HopscotchMap<SNVCollectionCount, Long, SNVCollectionCount> variationCounts,
                                  final Reference reference ) {
            final List<Interval> refCoverage = getRefCoverage();
            if ( refCoverage.isEmpty() ) return;

            if ( snvList == null ) {
                moleculeCounts.bumpInconsistentPairs();
                return;
            }
            if ( snvList.stream().anyMatch(
                    snv -> snv.getQuality() < minQ || "-ACGT".indexOf(snv.getVariantCall()) == -1) ) {
                moleculeCounts.bumpNLowQualityVariant();
                return;
            }
            if ( !hasCleanFlanks(minFlankingLength, reference.getRefSeqLength()) ) {
                moleculeCounts.bumpInsufficientFlank();
                return;
            }

            final int coverage = reference.updateCoverage(refCoverage);
            final Interval totalCoverage = new Interval(getFirstRefIndex(), getLastRefIndex());
            reference.updateSpan(totalCoverage);

            if ( snvList.isEmpty() ) {
                moleculeCounts.bumpNWildType();
                codonTracker.reportWildCodonCounts(totalCoverage);
            } else {
                moleculeCounts.bumpCalledVariant();
                codonTracker.reportVariantCodonCounts(totalCoverage, codonTracker.encodeSNVsAsCodons(snvList));
                final SNVCollectionCount newVal = new SNVCollectionCount(snvList, coverage);
                final SNVCollectionCount oldVal = variationCounts.find(newVal);
                if ( oldVal != null ) {
                    oldVal.bumpCount(coverage);
                } else {
                    variationCounts.add(newVal);
                }
            }
        }

        protected ReadReport( final List<Interval> refCoverage, final List<SNV> snvList ) {
            this.refCoverage = refCoverage;
            this.snvList = snvList;
        }

        public static ReadReport NULL_REPORT = new ReadReport(new ArrayList<>(), new ArrayList<>());
    }

    @VisibleForTesting final static class PairedReadReport extends ReadReport {
        final ReadReport report1;
        final ReadReport report2;

        public PairedReadReport( final ReadReport report1, final ReadReport report2 ) {
            super(combineCoverage(report1, report2), combineVariations(report1, report2));
            this.report1 = report1;
            this.report2 = report2;
        }

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
        public void updateCounts( final MoleculeCounts moleculeCounts,
                                  final CodonTracker codonTracker,
                                  final HopscotchMap<SNVCollectionCount, Long, SNVCollectionCount> variationCounts,
                                  final Reference reference ) {
            final int start1 = report1.getFirstRefIndex();
            final int end1 = report1.getLastRefIndex();
            final int start2 = report2.getFirstRefIndex();
            final int end2 = report2.getLastRefIndex();

            final int overlapStart = Math.max(start1, start2);
            final int overlapEnd = Math.min(end1, end2);
            if ( overlapStart > overlapEnd ) { // disjoint -- report each independently
                report1.updateCounts(moleculeCounts, codonTracker, variationCounts, reference);
                report2.updateCounts(moleculeCounts, codonTracker, variationCounts, reference);
            } else { // overlapping -- report molecule as one
                super.updateCounts(moleculeCounts, codonTracker, variationCounts, reference);
            }
        }

        private static List<Interval> combineCoverage( final ReadReport report1, final ReadReport report2 ) {
            final List<Interval> refCoverage1 = report1.getRefCoverage();
            final List<Interval> refCoverage2 = report2.getRefCoverage();
            if ( refCoverage1.isEmpty() ) return refCoverage2;
            if ( refCoverage2.isEmpty() ) return refCoverage1;

            final List<Interval> combinedCoverage =
                    new ArrayList<>(refCoverage1.size() + refCoverage2.size());
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

        // if the reports have overlapping coverage, check that the SNVs agree in the overlapping part,
        // and if they don't, return a null snvList.
        // if there's no overlap, just glue the two lists together
        private static List<SNV> combineVariations( final ReadReport report1, final ReadReport report2 ) {
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
                        return null; // list 2 has a SNV in the overlapping region that isn't present in list 1
                    }
                } else if ( snv2 == null ) {
                    next = snv1;
                    snv1 = iterator1.hasNext() ? iterator1.next() : null;
                    final int refIndex = next.getRefIndex();
                    if ( refIndex >= overlapStart && refIndex < overlapEnd ) {
                        return null; // list 1 has a SNV in the overlapping region that isn't present in list 2
                    }
                } else {
                    final int refIndex1 = snv1.getRefIndex();
                    final int refIndex2 = snv2.getRefIndex();
                    if ( refIndex1 < refIndex2 ) {
                        next = snv1;
                        snv1 = iterator1.hasNext() ? iterator1.next() : null;
                        if ( refIndex1 >= overlapStart && refIndex1 < overlapEnd ) {
                            return null; // list 1 has an earlier SNV in the overlapping region than isn't in list 2
                        }
                    } else if ( refIndex2 < refIndex1 ) {
                        next = snv2;
                        snv2 = iterator2.hasNext() ? iterator2.next() : null;
                        if ( refIndex2 >= overlapStart && refIndex2 < overlapEnd ) {
                            return null; // list 2 has an earlier SNV in the overlapping region than isn't in list 1
                        }
                    } else if ( !snv1.equals(snv2) ) {
                        return null; // the SNVs at this position aren't equivalent -- the lists disagree
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

    @VisibleForTesting static ReadReport getReadReport( final GATKRead read,
                                                        final byte[] refSeq,
                                                        final ReadCounts readCounts,
                                                        final Reference reference ) {
        readCounts.bumpNReads();
        readCounts.addBaseCalls(read.getLength());

        if ( read.isUnmapped() || read.isDuplicate() || read.failsVendorQualityCheck() ) {
            readCounts.bumpNReadsUnmapped();
            return ReadReport.NULL_REPORT;
        }

        final Interval trim = calculateTrim(read.getBaseQualitiesNoCopy(), readCounts);
        if ( trim.size() == 0 ) {
            return ReadReport.NULL_REPORT;
        }

        return new ReadReport(read, trim, reference.getRefSeq());
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
            return Interval.NULL_INTERVAL;
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

    private static ReadReport combineReports( final ReadReport report1, final ReadReport report2 ) {
        if ( report1.getRefCoverage().isEmpty() ) return report2;
        if ( report2.getRefCoverage().isEmpty() ) return report1;
        return new PairedReadReport(report1, report2);
    }
}
