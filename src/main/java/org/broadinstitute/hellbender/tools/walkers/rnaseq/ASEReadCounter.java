package org.broadinstitute.hellbender.tools.walkers.rnaseq;


import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Calculate read counts per allele for allele-specific expression analysis of RNAseq data
 * <p>
 * <p>
 * This tool calculates allele counts at a set of positions after applying filters that are tuned for enabling
 * allele-specific expression (ASE) analysis of RNAseq data. The filters operate on mapping quality, base quality, depth of coverage,
 * overlapping paired reads and deletions overlapping the position. All thresholds and options are controlled by
 * command-line arguments.
 * </p>
 * <p>
 * <h3>Input</h3>
 * <ul>
 * <li>BAM files (with proper headers) to be analyzed for ASE</li>
 * <li>A VCF file with specific sites to process.</li>
 * </ul>
 * <p>
 * <h3>Output</h3>
 * <p>
 * A table of allele counts at the given sites. By default, it is formatted as a tab-delimited text file
 * that is readable by R and compatible with <a href="http://www.well.ox.ac.uk/~rivas/mamba/">Mamba</a>,
 * a downstream tool developed for allele-specific expression analysis.
 * </p>
 * <h3>Usage Example</h3>
 * <p>
 *     gatk ASEReadCounter \
 *     -R Homo_sapiens_assembly38.fasta \
 *     -I input.bam \
 *     -V sites.vcf.gz \
 *     -O output.table
 * </p>
 * <p>
 * <h3>Note</h3>
 * <ul>
 * <li>Like most GATK tools, this tools filters out duplicate reads by default. However, some ASE methods
 * recommend including duplicate reads in the analysis, so the DuplicateRead filter can be disabled using the
 * "-DF NotDuplicateReadFilter" flag in the command-line.</li>
 * </ul>
 * <h3>Caveat</h3>
 * <ul>
 * <li>This tool will only process biallelic het SNP sites. If your callset contains multiallelic sites, they will be ignored.
 * Optionally, you can subset your callset to just biallelic variants using e.g.
 * SelectVariants
 * with the option "-restrictAllelesTo BIALLELIC".</li>
 * </ul>
 * <p>
 * For more details see <a href="http://biorxiv.org/content/biorxiv/early/2015/03/05/016097.full.pdf">Castel, S. et al. Tools and Best Practices for allelic expression analysis.</a>
 *
 * @author Ami Levy Moonshine
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = "Counts filtered reads at het sites for allele specific expression estimate",
        oneLineSummary = "Generates table of filtered base counts at het sites for allele specific expression",
        programGroup = CoverageAnalysisProgramGroup.class
)
public class ASEReadCounter extends LocusWalker {

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", optional = true)
    private File outputFile = null;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "One or more VCF files")
    private List<FeatureInput<VariantContext>> variants;

    /**
     * If this argument is enabled, loci with total depth lower than this threshold after all filters have been applied
     * will be skipped. This can be set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "min-depth-of-non-filtered-base", shortName = "min-depth", doc = "Minimum number of bases that pass filters", optional = true)
    public int minDepthOfNonFilteredBases = -1;

    /**
     * If this argument is enabled, reads with mapping quality values lower than this threshold will not be counted.
     * This can be set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "min-mapping-quality", shortName = "mmq", doc = "Minimum read mapping quality", optional = true)
    public int minMappingQuality = 0;

    /**
     * If this argument is enabled, bases with quality scores lower than this threshold will not be counted.
     * This can be set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "min-base-quality", shortName = "mbq", doc = "Minimum base quality", optional = true)
    public byte minBaseQuality = 0;

    /**
     * These options modify how the tool deals with overlapping read pairs. The default value is COUNT_FRAGMENTS_REQUIRE_SAME_BASE.
     */
    @Argument(fullName = "count-overlap-reads-handling", shortName = "overlap", doc = "Handling of overlapping reads from the same fragment", optional = true)
    public CountPileupType countType = CountPileupType.COUNT_FRAGMENTS_REQUIRE_SAME_BASE;

    /**
     * Available options are csv, table, rtable. By default, the format is rtable (an r-readable table).
     */
    @Argument(fullName = "output-format", doc = "Format of the output file", optional = true)
    public OUTPUT_FORMAT outputFormat = OUTPUT_FORMAT.RTABLE;

    public String separator = "\t";

    public enum OUTPUT_FORMAT {
        TABLE,
        RTABLE,
        CSV
    }

    public enum CountPileupType {
        /**
         * Count all reads independently (even if from the same fragment).
         */
        COUNT_READS,
        /**
         * Count all fragments (even if the reads that compose the fragment are not consistent at that base).
         */
        COUNT_FRAGMENTS,
        /**
         * Count all fragments (but only if the reads that compose the fragment are consistent at that base).
         */
        COUNT_FRAGMENTS_REQUIRE_SAME_BASE
    }

    private PrintStream outputStream = null;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> defaultFilters = new ArrayList<>();
        defaultFilters.add(ReadFilterLibrary.VALID_ALIGNMENT_START);
        defaultFilters.add(ReadFilterLibrary.VALID_ALIGNMENT_END);
        defaultFilters.add(ReadFilterLibrary.HAS_READ_GROUP);
        defaultFilters.add(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS);
        defaultFilters.add(ReadFilterLibrary.SEQ_IS_STORED);
        defaultFilters.add(ReadFilterLibrary.NOT_DUPLICATE);
        defaultFilters.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
        defaultFilters.add(new ReadFilterLibrary.MappedReadFilter());
        return defaultFilters;
    }

    @Override
    public void onTraversalStart() {
        try {
            outputStream = outputFile != null ? new PrintStream(outputFile) : System.out;
        } catch (FileNotFoundException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }

        if (outputFormat.equals(OUTPUT_FORMAT.CSV)) {
            separator = ",";
        }
        final String header = "contig" + separator + "position" + separator + "variantID" + separator + "refAllele" + separator + "altAllele" + separator + "refCount" + separator + "altCount" + separator + "totalCount" + separator + "lowMAPQDepth" + separator + "lowBaseQDepth" + separator + "rawDepth" + separator + "otherBases" + separator + "improperPairs";
        outputStream.println(header);

        variants.stream()
                .filter(vFeature -> ((VCFHeader)getHeaderForFeatures(vFeature)).getGenotypeSamples().isEmpty())
                .forEach(vFeature -> logger.warn("\n#####################################################################################################\nVariant input file "+vFeature.getName()+" lacks genotype fields. Variants from this file will produce no results in the output.\n#####################################################################################################"));

    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final String contig = alignmentContext.getContig();
        final long position = alignmentContext.getPosition();

        final char refAllele = (char) referenceContext.getBase();

        final List<VariantContext> VCs = featureContext.getValues(variants);
        if (VCs != null && VCs.size() > 1) {
            throw new UserException("More then one variant context at position: " + contig + ":" + position);
        }
        if (VCs == null || VCs.isEmpty()) {
            return;
        }

        final VariantContext vc = VCs.get(0);
        if (!vc.isBiallelic()) {
            logger.warn("Ignoring site: cannot run ASE on non-biallelic sites: " + vc.toString());
            return;
        }

        if (vc.getHetCount() < 1) {
            logger.warn("Ignoring site: variant is not het at postion: " + contig + ":" + position);
            return;
        }

        if (vc.getNAlleles() == 1 || vc.getAlternateAllele(0).getBases().length == 0) {
            throw new UserException("The file of variant sites must contain heterozygous sites and cannot be a GVCF file containing <NON_REF> alleles.");
        }

        final char altAllele = (char) vc.getAlternateAllele(0).getBases()[0];

        final String siteID = vc.getID();
        final ReadPileup pileup = filterPileup(alignmentContext.getBasePileup(), countType);

        // count up the depths of all and QC+ bases
        final String line = calculateLineForSite(pileup, siteID, refAllele, altAllele);
        if (line != null) {
            outputStream.println(line);
        }
    }

    @Override
    public void closeTool() {
        if (outputStream != null)
            outputStream.close();
    }

    private ReadPileup filterPileup(final ReadPileup originalPileup, final CountPileupType countType) {
        SAMFileHeader header = getHeaderForReads();

        final ReadPileup pileupWithDeletions;
        switch (countType) {
            case COUNT_FRAGMENTS_REQUIRE_SAME_BASE: pileupWithDeletions = originalPileup.getOverlappingFragmentFilteredPileup(true, ReadPileup.baseQualTieBreaker, header);
                break;
            case COUNT_READS: pileupWithDeletions = originalPileup;
                break;
            case COUNT_FRAGMENTS: pileupWithDeletions = originalPileup.getOverlappingFragmentFilteredPileup(false, ReadPileup.baseQualTieBreaker, header);
                break;
            default: throw new UserException("Must use valid CountPileupType");
        }

        return pileupWithDeletions.makeFilteredPileup(p -> !p.isDeletion());

    }

    private String calculateLineForSite(final ReadPileup pileup, final String siteID, final char refAllele, final char altAllele) {

        int rawDepth = 0, lowBaseQDepth = 0, lowMAPQDepth = 0, refCount = 0, altCount = 0, totalNonFilteredCount = 0, otherBasesCount = 0, improperPairsCount = 0;

        for (final PileupElement base : pileup) {
            rawDepth++;

            if (base.getRead().isPaired() && (base.getRead().mateIsUnmapped() || !base.getRead().isProperlyPaired())) {
                improperPairsCount++;
                continue;
            }
            if (base.getMappingQual() < minMappingQuality) {
                lowMAPQDepth++;
                continue;
            }

            if (base.getQual() < minBaseQuality) {
                lowBaseQDepth++;
                continue;
            }

            if (base.getBase() == refAllele) {
                refCount++;
            }
            else if (base.getBase() == altAllele) {
                altCount++;
            }
            else {
                otherBasesCount++;
                continue;
            }
            totalNonFilteredCount++;
        }

        if (totalNonFilteredCount < minDepthOfNonFilteredBases) {
            return null;
        }

        final StringBuilder line = new StringBuilder();
        line.append(pileup.getLocation().getContig()).append(separator);
        line.append(pileup.getLocation().getStart()).append(separator);
        line.append(siteID).append(separator);
        line.append(refAllele).append(separator);
        line.append(altAllele).append(separator);
        line.append(refCount).append(separator);
        line.append(altCount).append(separator);
        line.append(totalNonFilteredCount).append(separator);
        line.append(lowMAPQDepth).append(separator);
        line.append(lowBaseQDepth).append(separator);
        line.append(rawDepth).append(separator);
        line.append(otherBasesCount).append(separator);
        line.append(improperPairsCount);

        return line.toString();
    }
}
