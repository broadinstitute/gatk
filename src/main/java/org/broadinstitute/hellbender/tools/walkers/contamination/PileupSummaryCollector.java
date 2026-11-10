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
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
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
        final List<VariantContext> vcs = featureContext.getValues(variants);
        if (vcs.isEmpty()) {
            return;
        }
        final VariantContext vc = vcs.get(0);

        if ( vc.isBiallelic() && vc.isSNP() && alleleFrequencyInRange(vc) ) {
            final ReadPileup pileup = alignmentContext.getBasePileup();
            for (final Table table : tables) {
                final ReadPileup tablePileup = pileup.makeFilteredPileup(pe -> table.readFilter().test(pe.getRead()));
                if (!tablePileup.isEmpty()) {
                    try {
                        table.writer().writeRecord(new PileupSummary(vc, tablePileup));
                    } catch (final IOException ex) {
                        throw new UserException(String.format("Encountered an IO exception while writing to %s", table.output()), ex);
                    }
                }
            }
        }
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
