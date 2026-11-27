package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;

/**
 * Tabulates {@link PileupSummary} records at common germline variant sites one locus at a time, for {@link GetPileupSummaries}
 * and for Mutect2, which collects the pileup summaries of its tumor and normal as it traverses their reads.
 */
public final class PileupSummaryCollector implements AutoCloseable {
    private static final Logger logger = LogManager.getLogger(PileupSummaryCollector.class);

    private final FeatureInput<VariantContext> variants;
    private final double minPopulationAlleleFrequency;
    private final double maxPopulationAlleleFrequency;

    private final List<Table> tables = new ArrayList<>();

    private boolean sawVariantsWithoutAlleleFrequency = false;
    private boolean sawVariantsWithAlleleFrequency = false;

    // The variants overlapping a window of the genome, windowContig:[windowStart, windowStart + VARIANT_WINDOW_SIZE),
    // queried once for all the loci in the window, which come in order, rather than at every locus, and the positions of
    // the window that they overlap, relative to its start.  Most loci overlap no variants.
    private static final int VARIANT_WINDOW_SIZE = 10_000;
    private String windowContig = null;
    private int windowStart = 0;
    private List<VariantContext> windowVariants = Collections.emptyList();
    private final BitSet windowPositionsWithVariants = new BitSet(VARIANT_WINDOW_SIZE);

    public PileupSummaryCollector(final FeatureInput<VariantContext> variants, final double minPopulationAlleleFrequency,
                                  final double maxPopulationAlleleFrequency) {
        this.variants = variants;
        this.minPopulationAlleleFrequency = minPopulationAlleleFrequency;
        this.maxPopulationAlleleFrequency = maxPopulationAlleleFrequency;
    }

    /**
     * Write the pileup summaries of the reads that pass a filter, for example the reads of one sample, to a table
     *
     * @param sample the sample recorded in the table's metadata
     */
    public void addTable(final File output, final String sample, final Predicate<GATKRead> readFilter) {
        try {
            final PileupSummary.PileupSummaryTableWriter writer = new PileupSummary.PileupSummaryTableWriter(IOUtils.fileToPath(output));
            writer.writeMetadata(TableUtils.SAMPLE_METADATA_TAG, sample);
            tables.add(new Table(output, writer, readFilter));
        } catch (IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(output, ex);
        }
    }

    /**
     * Write a record to each table if this locus is a biallelic SNP among the variants, with population allele frequency in range
     */
    public void process(final AlignmentContext alignmentContext, final FeatureContext featureContext) {
        final List<VariantContext> vcs = variantsAt(featureContext);
        if (vcs.isEmpty()) {
            return;
        }
        final VariantContext vc = vcs.get(0);

        if ( vc.isBiallelic() && vc.isSNP() && alleleFrequencyInRange(vc) ) {
            final ReadPileup pileup = alignmentContext.getBasePileup();
            for (int table = 0; table < tables.size(); table++) {
                final Predicate<GATKRead> readFilter = tables.get(table).readFilter();
                final ReadPileup tablePileup = pileup.makeFilteredPileup(pe -> readFilter.test(pe.getRead()));
                if (!tablePileup.isEmpty()) {
                    writeRecord(table, new PileupSummary(vc, tablePileup));
                }
            }
        }
    }

    /**
     * The variants at which {@link #process} might write records, for callers that count the bases of each table's
     * pileups themselves, as Mutect2 does one read at a time: at each locus of an interval, the first variant overlapping
     * it, which process would take, if it is a biallelic SNP.  Whether a site's population allele frequency is in range
     * is left to {@link #summarizes}, because process checks it, and warns when it is missing, only at the loci it is given.
     *
     * @param intervalContext the features of an interval
     * @return the sites in the interval, in order of position
     */
    public List<VariantContext> candidateSites(final FeatureContext intervalContext) {
        final SimpleInterval interval = intervalContext.getInterval();
        final List<VariantContext> sites = new ArrayList<>();
        // the variants come in order of their starts, and the last locus that the variants so far overlap
        int lastOverlappedLocus = Integer.MIN_VALUE;
        for (final VariantContext vc : intervalContext.getValues(variants)) {
            if (vc.getStart() > lastOverlappedLocus && vc.getStart() >= interval.getStart() && vc.isBiallelic() && vc.isSNP()) {
                sites.add(vc);
            }
            lastOverlappedLocus = Math.max(lastOverlappedLocus, vc.getEnd());
        }
        return sites;
    }

    /**
     * Whether a site from {@link #candidateSites} has pileup summaries, that is, whether its population allele frequency
     * is in range.  Like process, ask only at loci that are traversed.
     */
    public boolean summarizes(final VariantContext site) {
        return alleleFrequencyInRange(site);
    }

    public int numberOfTables() {
        return tables.size();
    }

    /**
     * Whether the pileups of a table include a read
     */
    public boolean tableIncludes(final int table, final GATKRead read) {
        return tables.get(table).readFilter().test(read);
    }

    /**
     * Write a table's record at a site from {@link #candidateSites} that it {@link #summarizes}, from the counts of A, C,
     * G and T in the table's pileup there, which process would write if the pileup has elements, including deletions
     */
    public void writeRecord(final int table, final VariantContext site, final int[] baseCounts) {
        writeRecord(table, new PileupSummary(site, baseCounts));
    }

    private void writeRecord(final int table, final PileupSummary record) {
        try {
            tables.get(table).writer().writeRecord(record);
        } catch (final IOException ex) {
            throw new UserException(String.format("Encountered an IO exception while writing to %s", tables.get(table).output()), ex);
        }
    }

    /**
     * The variants overlapping a FeatureContext's locus, as {@code featureContext.getValues(variants)} would give them,
     * in the order of the window's query, which filters the same file-ordered cache of the variants
     */
    private List<VariantContext> variantsAt(final FeatureContext featureContext) {
        final SimpleInterval locus = featureContext.getInterval();
        if (locus == null || locus.size() != 1) {
            return featureContext.getValues(variants);
        }

        final int position = locus.getStart();
        if (!locus.getContig().equals(windowContig) || position < windowStart || position >= windowStart + VARIANT_WINDOW_SIZE) {
            windowContig = locus.getContig();
            windowStart = position;
            windowVariants = featureContext.getValues(variants, 0, VARIANT_WINDOW_SIZE - 1);
            windowPositionsWithVariants.clear();
            for (final VariantContext vc : windowVariants) {
                final int from = Math.max(vc.getStart(), windowStart) - windowStart;
                final int to = Math.min(vc.getEnd(), windowStart + VARIANT_WINDOW_SIZE - 1) - windowStart;
                if (from <= to) {
                    windowPositionsWithVariants.set(from, to + 1);
                }
            }
        }

        if (!windowPositionsWithVariants.get(position - windowStart)) {
            return Collections.emptyList();
        }
        final List<VariantContext> result = new ArrayList<>(1);
        for (final VariantContext vc : windowVariants) {
            if (vc.getStart() <= position && position <= vc.getEnd()) {
                result.add(vc);
            }
        }
        return result;
    }

    /**
     * Fail if the variants we saw lack allele frequencies entirely, as opposed to only at some sites, which merely warrants a warning
     */
    public void validateAlleleFrequencies() {
        if (sawVariantsWithoutAlleleFrequency && !sawVariantsWithAlleleFrequency) {
            throw new UserException.BadInput("No variants in population vcf had an allele frequency (AF) field.");
        }
    }

    /**
     * Fail if the variants' header lacks an allele frequency (AF) INFO field
     */
    public static void validateHeader(final VCFHeader header) {
        final boolean alleleFrequencyInHeader = header.getInfoHeaderLines().stream()
                .anyMatch(line -> line.getID().equals(VCFConstants.ALLELE_FREQUENCY_KEY));
        if (!alleleFrequencyInHeader) {
            throw new UserException.BadInput("Population vcf does not have an allele frequency (AF) info field in its header.");
        }
    }

    @Override
    public void close() {
        for (final Table table : tables) {
            try {
                table.writer().close();
            } catch (IOException ex) {
                throw new UserException(String.format("Encountered an IO exception while closing %s", table.output()), ex);
            }
        }
    }

    private boolean alleleFrequencyInRange(final VariantContext vc) {
        if (!vc.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
            if (!sawVariantsWithoutAlleleFrequency) {
                logger.warn(String.format("Variant context at %s:%d lacks allele frequency (AF) field.", vc.getContig(), vc.getStart()));
                sawVariantsWithoutAlleleFrequency = true;
            }
            return false;
        } else {
            sawVariantsWithAlleleFrequency = true;
            final double alleleFrequency = vc.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, -1.0);
            return minPopulationAlleleFrequency < alleleFrequency && alleleFrequency < maxPopulationAlleleFrequency;
        }
    }

    private record Table(File output, PileupSummary.PileupSummaryTableWriter writer, Predicate<GATKRead> readFilter) { }
}
