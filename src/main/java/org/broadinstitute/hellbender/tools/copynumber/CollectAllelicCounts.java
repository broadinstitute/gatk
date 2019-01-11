package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.datacollection.AllelicCountCollector;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.MetadataUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.utils.Nucleotide;

import java.io.File;
import java.util.List;

/**
 * Collects reference and alternate allele counts at specified sites. The alt count is defined as the
 * total count minus the ref count, and the alt nucleotide is defined as the non-ref base with the highest count,
 * with ties broken by the order of the bases in {@link AllelicCountCollector#BASES}. Only reads that pass the
 * specified read filters and bases that exceed the specified {@code minimum-base-quality} will be counted.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SAM format read data
 *     </li>
 *     <li>
 *         Reference FASTA file
 *     </li>
 *     <li>
 *         Sites at which allelic counts will be collected
 *     </li>
 * </ul>
 *
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Allelic-counts file.
 *         This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link AllelicCountCollection.AllelicCountTableColumn},
 *         and the corresponding entry rows.
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk CollectAllelicCounts \
 *          -I sample.bam \
 *          -R reference.fa \
 *          -L sites.interval_list \
 *          -O sample.allelicCounts.tsv
 * </pre>
 *
 * @author Lee Lichtenstein &lt;lichtens@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Collects reference and alternate allele counts at specified sites",
        oneLineSummary = "Collects reference and alternate allele counts at specified sites",
        programGroup = CoverageAnalysisProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class CollectAllelicCounts extends LocusWalker {
    private static final Logger logger = LogManager.getLogger(CollectAllelicCounts.class);

    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 30;

    public static final String MINIMUM_BASE_QUALITY_LONG_NAME = "minimum-base-quality";

    @Argument(
            doc = "Output file for allelic counts.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputAllelicCountsFile;

    @Argument(
            doc = "Minimum base quality.  Base calls with lower quality will be filtered out of pileups.",
            fullName = MINIMUM_BASE_QUALITY_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int minimumBaseQuality = 20;

    private AllelicCountCollector allelicCountCollector;

    @Override
    public boolean emitEmptyLoci() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        final SampleLocatableMetadata metadata = MetadataUtils.fromHeader(getHeaderForReads(), Metadata.Type.SAMPLE_LOCATABLE);
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        if (!CopyNumberArgumentValidationUtils.isSameDictionary(metadata.getSequenceDictionary(), sequenceDictionary)) {
            logger.warn("Sequence dictionary in BAM does not match the master sequence dictionary.");
        }
        allelicCountCollector = new AllelicCountCollector(metadata);
        logger.info("Collecting allelic counts...");
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = super.getDefaultReadFilters();
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));
        return filters;
    }

    @Override
    public Object onTraversalSuccess() {
        allelicCountCollector.getAllelicCounts().write(outputAllelicCountsFile);
        logger.info("Allelic counts written to " + outputAllelicCountsFile);
        return("SUCCESS");
    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final byte refAsByte = referenceContext.getBase();
        allelicCountCollector.collectAtLocus(Nucleotide.decode(refAsByte), alignmentContext.getBasePileup(), alignmentContext.getLocation(), minimumBaseQuality);
    }
}
