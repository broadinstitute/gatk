package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.dragstr.BufferedReferenceBases;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

/**
 * Collects training data for a model of polymerase slippage, the PCR and sequencing errors that change the number of
 * repeat units in short tandem repeats (STRs).
 *
 * <p>
 *     An STR is a run of at least {@code --min-repeat-units} whole copies of a repeat unit of at most
 *     {@code --max-repeat-unit-length} bases, at least {@code --min-repeat-length} bases long in total.  Units are made of
 *     A, C, G and T and are not themselves copies of a shorter unit.  STRs never overlap: the reference is scanned from left
 *     to right, the longest STR starting at each position is taken (the one with the shorter unit in a tie), and scanning
 *     resumes after its end.  A partial copy of the unit at the end of a run is not part of the STR.
 * </p>
 * <p>
 *     The tool makes a single pass over coordinate-sorted reads, finding the STRs in the reference just ahead of them.  A
 *     fragment counts toward an STR if one of its reads is aligned to the reference bases {@code --flanking-bases} away
 *     from the STR on each side.  Its change in length is the number of read bases between those two aligned bases minus
 *     the number of reference bases between them, so it includes every indel in the STR or its flanks wherever the aligner
 *     put it.  A fragment counts once even if both its reads span the STR, using whichever read comes first.  Reads are
 *     filtered as in Mutect2.
 * </p>
 * <p>
 *     The output table has a row for each STR spanned by at least one fragment: its location, repeat unit and number of
 *     units, and the number of fragments that change its length by at most -4 units, -3, -2, -1, 0, 1, 2, 3 and at least 4
 *     units, and by a length that isn't a whole number of units.
 * </p>
 * <p>
 *     Each row also gives the two non-zero changes of a whole number of units with the most fragments, with ties going to
 *     the smaller change and then to the deletion, and their population allele frequencies in the
 *     {@code --germline-resource}.  These help to tell slippage from germline variation.  A change of 0 means there is no
 *     such allele, and the allele frequency is NaN if the resource doesn't have the allele or no resource was given.  A
 *     resource allele matches if it changes the STR's length by the same number of units, however the resource
 *     represents it.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk CollectSTRSlippageCounts \
 *   -R reference.fasta \
 *   -I normal.bam \
 *   --germline-resource af-only-gnomad.vcf.gz \
 *   -O slippage-counts.table
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Counts the fragments that change the length of each short tandem repeat (STR) by each number of repeat units, " +
                "as training data for a model of polymerase slippage",
        oneLineSummary = "Collects changes in STR length in reads for training a polymerase slippage model",
        programGroup = CoverageAnalysisProgramGroup.class)
@DocumentedFeature
@ExperimentalFeature
public class CollectSTRSlippageCounts extends ReadWalker {
    public static final String MIN_REPEAT_LENGTH_LONG_NAME = "min-repeat-length";
    public static final String MIN_REPEAT_UNITS_LONG_NAME = "min-repeat-units";
    public static final String MAX_REPEAT_UNIT_LENGTH_LONG_NAME = "max-repeat-unit-length";
    public static final String FLANKING_BASES_LONG_NAME = "flanking-bases";

    public static final int DEFAULT_MIN_REPEAT_LENGTH = 7;
    public static final int DEFAULT_MIN_REPEAT_UNITS = 2;
    public static final int DEFAULT_MAX_REPEAT_UNIT_LENGTH = 5;
    public static final int DEFAULT_FLANKING_BASES = 5;

    // changes in length are binned by number of units, from at most -MAX_BINNED_CHANGE to at least +MAX_BINNED_CHANGE,
    // followed by one bin for changes that are not a whole number of units
    public static final int MAX_BINNED_CHANGE = 4;
    public static final int NOT_WHOLE_UNITS_BIN = 2 * MAX_BINNED_CHANGE + 1;
    public static final List<String> COUNT_COLUMNS = List.of("minus_4_or_more", "minus_3", "minus_2", "minus_1",
            "no_change", "plus_1", "plus_2", "plus_3", "plus_4_or_more", "not_whole_units");

    // the most common non-zero changes of whole units, and their allele frequencies in the germline resource
    public static final int ANNOTATED_ALTS = 2;
    public static final List<String> ALT_COLUMNS = List.of("top_alt_change", "top_alt_af", "second_alt_change", "second_alt_af");

    private static final int REFERENCE_BUFFER_SIZE = 100_000;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "The output table")
    private GATKPath outputTable;

    @Argument(fullName = M2ArgumentCollection.GERMLINE_RESOURCE_LONG_NAME, optional = true,
            doc = "Population vcf of germline sequencing containing allele frequencies, for annotating the most common changes in each STR's length")
    private FeatureInput<VariantContext> germlineResource;

    @Argument(fullName = MIN_REPEAT_LENGTH_LONG_NAME, minValue = 1, optional = true,
            doc = "Minimum total length in bases of an STR")
    private int minRepeatLength = DEFAULT_MIN_REPEAT_LENGTH;

    @Argument(fullName = MIN_REPEAT_UNITS_LONG_NAME, minValue = 2, optional = true,
            doc = "Minimum number of copies of the repeat unit in an STR")
    private int minRepeatUnits = DEFAULT_MIN_REPEAT_UNITS;

    @Argument(fullName = MAX_REPEAT_UNIT_LENGTH_LONG_NAME, minValue = 1, optional = true,
            doc = "Maximum length in bases of an STR's repeat unit")
    private int maxRepeatUnitLength = DEFAULT_MAX_REPEAT_UNIT_LENGTH;

    @Argument(fullName = FLANKING_BASES_LONG_NAME, minValue = 1, optional = true,
            doc = "Number of reference bases on each side of an STR to which a read must be aligned for its fragment to count")
    private int flankingBases = DEFAULT_FLANKING_BASES;

    private ReferenceDataSource reference;
    private FeatureDataSource<VariantContext> germlineResourceSource;
    private SlippageCountsWriter writer;

    // the current contig's reference and scanner, and the STRs the scanner has found that reads might yet span, in order
    private BufferedReferenceBases contigReference;
    private RepeatScanner scanner;
    private final Deque<Repeat> pending = new ArrayDeque<>();
    private final Set<String> finishedContigs = new HashSet<>();
    private int lastReadStart;
    private long repeatsWritten = 0;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Mutect2Engine.makeStandardMutect2ReadFilters();
    }

    @Override
    public void onTraversalStart() {
        // our own readers, since walkers may not use the engine's outside of apply()
        reference = ReferenceDataSource.of(referenceArguments.getReferencePath());
        if (germlineResource != null) {
            germlineResourceSource = new FeatureDataSource<>(germlineResource, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES, VariantContext.class);
        }
        try {
            writer = new SlippageCountsWriter(outputTable.toPath());
            writer.writeMetadata(TableUtils.SAMPLE_METADATA_TAG, String.join(",", ReadUtils.getSamplesFromHeader(getHeaderForReads())));
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputTable, "could not create the STR slippage counts table", e);
        }
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if (scanner == null || !read.getContig().equals(scanner.getContig())) {
            startContig(read.getContig());
        } else if (read.getStart() < lastReadStart) {
            throw new UserException.BadInput(String.format("Reads must be sorted by coordinate, but read %s at %s:%d follows a read starting at %d",
                    read.getName(), read.getContig(), read.getStart(), lastReadStart));
        }
        lastReadStart = read.getStart();

        // find every STR that could fit in this read with its flanks
        scanner.scanThrough(read.getEnd() - flankingBases - minRepeatLength + 1, pending::addLast);

        // reads come in order of start, so once a read starts too late to flank an STR no later read can flank it
        while (!pending.isEmpty() && pending.peekFirst().getStart() - flankingBases < read.getStart()) {
            write(pending.pollFirst());
        }

        // pending STRs don't overlap, so their ends are in order too
        for (final Repeat repeat : pending) {
            if (repeat.getEnd() + flankingBases > read.getEnd()) {
                break;
            }
            repeat.addRead(read, flankingBases);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (scanner != null) {
            finishContig();
        }
        return repeatsWritten;
    }

    @Override
    public void closeTool() {
        if (reference != null) {
            reference.close();
        }
        if (germlineResourceSource != null) {
            germlineResourceSource.close();
        }
        if (writer != null) {
            try {
                writer.close();
            } catch (final IOException e) {
                throw new UserException.CouldNotCreateOutputFile(outputTable, "could not close the STR slippage counts table", e);
            }
        }
    }

    private void startContig(final String contig) {
        if (scanner != null) {
            finishContig();
        }
        if (finishedContigs.contains(contig)) {
            throw new UserException.BadInput("Reads must be sorted by coordinate, but reads on contig " + contig + " are not together");
        }
        contigReference = BufferedReferenceBases.of(reference, contig, REFERENCE_BUFFER_SIZE);
        scanner = new RepeatScanner(contigReference, maxRepeatUnitLength, minRepeatUnits, minRepeatLength);
        lastReadStart = 0;
    }

    private void finishContig() {
        pending.forEach(this::write);
        pending.clear();
        finishedContigs.add(scanner.getContig());
    }

    private void write(final Repeat repeat) {
        if (repeat.getFragmentCount() == 0) {
            return;
        }
        final List<Integer> alts = repeat.mostCommonAlts(ANNOTATED_ALTS);
        final List<VariantContext> resourceVariants = germlineResourceSource == null || alts.isEmpty() ? Collections.emptyList() :
                germlineResourceSource.queryAndPrefetch(resourceInterval(contigReference, repeat));
        final List<Double> alleleFrequencies = alts.stream()
                .map(alt -> resourceAlleleFrequency(contigReference, repeat, resourceVariants, alt))
                .collect(Collectors.toList());
        try {
            writer.writeRecord(new Row(repeat, alts, alleleFrequencies));
            repeatsWritten++;
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputTable, "could not write to the STR slippage counts table", e);
        }
    }

    /**
     * The bin of a change in an STR's length: from 0 for at most -{@link #MAX_BINNED_CHANGE} units to
     * 2 * {@link #MAX_BINNED_CHANGE} for at least +{@link #MAX_BINNED_CHANGE} units, or {@link #NOT_WHOLE_UNITS_BIN}
     */
    @VisibleForTesting
    static int bin(final int changeInBases, final int unitLength) {
        if (changeInBases % unitLength != 0) {
            return NOT_WHOLE_UNITS_BIN;
        }
        final int changeInUnits = Math.max(-MAX_BINNED_CHANGE, Math.min(MAX_BINNED_CHANGE, changeInBases / unitLength));
        return changeInUnits + MAX_BINNED_CHANGE;
    }

    /**
     * The interval in which a germline resource would put a change in an STR's length: resources put indels as far left
     * as they can go, which is before the STR if the bases before it continue its pattern, after an anchor base
     */
    @VisibleForTesting
    static SimpleInterval resourceInterval(final BufferedReferenceBases contigReference, final Repeat repeat) {
        int anchor = Math.max(repeat.getStart() - 1, 1);
        while (anchor > 1 && Nucleotide.same(contigReference.byteAt(anchor), contigReference.byteAt(anchor + repeat.getUnit().length()))) {
            anchor--;
        }
        return new SimpleInterval(repeat.getContig(), anchor, repeat.getEnd());
    }

    /**
     * The population allele frequency in a germline resource of the allele that changes an STR's length by
     * {@code changeInUnits}, or NaN if the resource doesn't have it.  A resource allele matches if it turns the reference
     * into the same sequence, however the resource represents it.
     */
    @VisibleForTesting
    static double resourceAlleleFrequency(final BufferedReferenceBases contigReference, final Repeat repeat,
                                          final List<VariantContext> resourceVariants, final int changeInUnits) {
        if (resourceVariants.isEmpty() || -changeInUnits > repeat.getUnits()) {
            return Double.NaN;
        }

        // compare sequences over a window covering the STR and every resource variant, since multiallelic records can
        // extend past the STR
        final int windowStart = Math.min(repeat.getStart(), resourceVariants.stream().mapToInt(VariantContext::getStart).min().getAsInt());
        final int windowEnd = Math.max(repeat.getEnd(), resourceVariants.stream().mapToInt(VariantContext::getEnd).max().getAsInt());
        final byte[] window = new byte[windowEnd - windowStart + 1];
        contigReference.copyBytesAt(windowStart, window, 0, window.length);
        final byte[] upperCaseWindow = new String(window).toUpperCase().getBytes();

        final int offsetOfRepeat = repeat.getStart() - windowStart;
        final String unit = repeat.getUnit();
        final byte[] expected = changeInUnits < 0 ? replace(upperCaseWindow, offsetOfRepeat, -changeInUnits * unit.length(), new byte[0]) :
                replace(upperCaseWindow, offsetOfRepeat, 0, unit.repeat(changeInUnits).getBytes());
        for (final VariantContext variant : resourceVariants) {
            final int offset = variant.getStart() - windowStart;
            final byte[] ref = variant.getReference().getBases();
            if (!Arrays.equals(ref, 0, ref.length, upperCaseWindow, offset, offset + ref.length)) {
                continue;   // the resource disagrees with the reference
            }
            final List<Double> alleleFrequencies = variant.getAttributeAsDoubleList(VCFConstants.ALLELE_FREQUENCY_KEY, Double.NaN);
            for (int n = 0; n < variant.getAlternateAlleles().size() && n < alleleFrequencies.size(); n++) {
                if (Arrays.equals(replace(upperCaseWindow, offset, ref.length, variant.getAlternateAllele(n).getBases()), expected)) {
                    return alleleFrequencies.get(n);
                }
            }
        }
        return Double.NaN;
    }

    // the bases with length of them starting at offset replaced by the replacement
    private static byte[] replace(final byte[] bases, final int offset, final int length, final byte[] replacement) {
        final byte[] result = new byte[bases.length - length + replacement.length];
        System.arraycopy(bases, 0, result, 0, offset);
        System.arraycopy(replacement, 0, result, offset, replacement.length);
        System.arraycopy(bases, offset + length, result, offset + replacement.length, bases.length - offset - length);
        return result;
    }

    /**
     * Finds the STRs of one contig from left to right
     */
    @VisibleForTesting
    static final class RepeatScanner {
        private final BufferedReferenceBases sequence;
        private final int maxUnitLength;
        private final int minUnits;
        private final int minLength;
        private int position = 1;   // the next position at which an STR may start

        RepeatScanner(final BufferedReferenceBases sequence, final int maxUnitLength, final int minUnits, final int minLength) {
            this.sequence = sequence;
            this.maxUnitLength = maxUnitLength;
            this.minUnits = minUnits;
            this.minLength = minLength;
        }

        String getContig() {
            return sequence.contigID();
        }

        /**
         * Give every STR starting at or before {@code target} that hasn't been given already to {@code consumer}, in order
         */
        void scanThrough(final int target, final Consumer<Repeat> consumer) {
            final long stop = Math.min(target, sequence.length());
            while (position <= stop) {
                final Repeat repeat = repeatAt(position);
                if (repeat == null) {
                    position++;
                } else {
                    consumer.accept(repeat);
                    position = repeat.getEnd() + 1;
                }
            }
        }

        /**
         * The longest qualifying STR starting at {@code start}, or null if there is none: for each unit length up to the
         * maximum, the run of whole copies of the unit starting there, preferring the shorter unit in a tie.  A unit that is
         * copies of a shorter unit, like AA, never wins, because the shorter unit's run is at least as long.
         */
        @VisibleForTesting
        Repeat repeatAt(final int start) {
            int bestUnitLength = 0;
            int bestUnits = 0;
            for (int unitLength = 1; unitLength <= maxUnitLength && start + unitLength - 1 <= sequence.length(); unitLength++) {
                // units must be made of A, C, G and T, and longer units contain this base too
                if (!Nucleotide.decode(sequence.byteAt(start + unitLength - 1)).isStandard()) {
                    break;
                }

                long runEnd = start + unitLength;   // exclusive
                while (runEnd <= sequence.length() && Nucleotide.same(sequence.byteAt(runEnd), sequence.byteAt(runEnd - unitLength))) {
                    runEnd++;
                }

                final int units = (int) ((runEnd - start) / unitLength);
                if (units >= minUnits && units * unitLength >= minLength && units * unitLength > bestUnits * bestUnitLength) {
                    bestUnitLength = unitLength;
                    bestUnits = units;
                }
            }
            if (bestUnits == 0) {
                return null;
            }
            final byte[] unit = new byte[bestUnitLength];
            sequence.copyBytesAt(start, unit, 0, bestUnitLength);
            return new Repeat(sequence.contigID(), start, new String(unit).toUpperCase(), bestUnits);
        }
    }

    /**
     * An STR in the reference, and the fragments that span it binned by their change in its length
     */
    @VisibleForTesting
    static final class Repeat {
        private final String contig;
        private final int start;    // 1-based
        private final String unit;
        private final int units;
        private final int[] counts = new int[COUNT_COLUMNS.size()];
        private final Map<Integer, Integer> altCounts = new HashMap<>();    // fragments by non-zero change of whole units
        private final Set<String> countedFragments = new HashSet<>();

        Repeat(final String contig, final int start, final String unit, final int units) {
            this.contig = contig;
            this.start = start;
            this.unit = unit;
            this.units = units;
        }

        String getContig() { return contig; }

        int getStart() { return start; }

        int getEnd() { return start + unit.length() * units - 1; }

        String getUnit() { return unit; }

        int getUnits() { return units; }

        int[] getCounts() { return counts; }

        int getFragmentCount() { return countedFragments.size(); }

        /**
         * Count the read's fragment if the read is aligned to the reference bases {@code flankingBases} away from the STR on
         * each side, and if the fragment's other read hasn't been counted already.
         *
         * @return whether the fragment was counted
         */
        boolean addRead(final GATKRead read, final int flankingBases) {
            final int leftAnchor = start - flankingBases;
            final int rightAnchor = getEnd() + flankingBases;
            final Pair<Integer, CigarOperator> left = ReadUtils.getReadIndexForReferenceCoordinate(read, leftAnchor);
            final Pair<Integer, CigarOperator> right = ReadUtils.getReadIndexForReferenceCoordinate(read, rightAnchor);
            if (left.getRight() == null || right.getRight() == null || !left.getRight().isAlignment() ||
                    !right.getRight().isAlignment() || !countedFragments.add(read.getName())) {
                return false;
            }
            final int changeInBases = (right.getLeft() - left.getLeft()) - (rightAnchor - leftAnchor);
            counts[bin(changeInBases, unit.length())]++;
            if (changeInBases != 0 && changeInBases % unit.length() == 0) {
                altCounts.merge(changeInBases / unit.length(), 1, Integer::sum);
            }
            return true;
        }

        /**
         * Up to {@code n} non-zero changes of whole units with the most fragments, most first, with ties going to the
         * smaller change and then to the deletion
         */
        List<Integer> mostCommonAlts(final int n) {
            return altCounts.entrySet().stream()
                    .sorted(Comparator.<Map.Entry<Integer, Integer>>comparingInt(entry -> -entry.getValue())
                            .thenComparingInt(entry -> Math.abs(entry.getKey()))
                            .thenComparingInt(Map.Entry::getKey))
                    .limit(n)
                    .map(Map.Entry::getKey)
                    .collect(Collectors.toList());
        }
    }

    // an STR with its most common non-zero changes of whole units and their allele frequencies in the germline resource
    private record Row(Repeat repeat, List<Integer> alts, List<Double> alleleFrequencies) { }

    private static final class SlippageCountsWriter extends TableWriter<Row> {
        private static final TableColumnCollection COLUMNS = new TableColumnCollection(columnNames());

        SlippageCountsWriter(final Path output) throws IOException {
            super(output, COLUMNS);
        }

        private static List<String> columnNames() {
            final List<String> names = new ArrayList<>(List.of("contig", "start", "end", "unit", "units"));
            names.addAll(COUNT_COLUMNS);
            names.addAll(ALT_COLUMNS);
            return names;
        }

        @Override
        protected void composeLine(final Row row, final DataLine dataLine) {
            final Repeat repeat = row.repeat();
            dataLine.set("contig", repeat.getContig())
                    .set("start", repeat.getStart())
                    .set("end", repeat.getEnd())
                    .set("unit", repeat.getUnit())
                    .set("units", repeat.getUnits());
            for (int bin = 0; bin < COUNT_COLUMNS.size(); bin++) {
                dataLine.set(COUNT_COLUMNS.get(bin), repeat.getCounts()[bin]);
            }
            for (int n = 0; n < ANNOTATED_ALTS; n++) {
                final boolean exists = n < row.alts().size();
                dataLine.set(ALT_COLUMNS.get(2 * n), exists ? row.alts().get(n) : 0)
                        .set(ALT_COLUMNS.get(2 * n + 1), exists ? row.alleleFrequencies().get(n) : Double.NaN);
            }
        }
    }
}
