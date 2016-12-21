/*
* Copyright 2012-2016 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.rnaseq;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.walkers.Downsample;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.tools.walkers.coverage.CoverageUtils;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.io.PrintStream;
import java.util.List;

/**
 * Calculate read counts per allele for allele-specific expression analysis
 *
 * <p>
 * This tool calculates allele counts at a set of positions after applying filters that are tuned for enabling
 * allele-specific expression (ASE) analysis. The filters operate on mapping quality, base quality, depth of coverage,
 * overlapping paired reads and deletions overlapping the position. All thresholds and options are controlled by
 * command-line arguments.
 * </p>
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>BAM files (with proper headers) to be analyzed for ASE</li>
 *     <li>A VCF file with specific sites to process.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <p>
 * A table of allele counts at the given sites. By default, it is formatted as a tab-delimited text file
 * that is readable by R and compatible with <a href="http://www.well.ox.ac.uk/~rivas/mamba/">Mamba</a>,
 * a downstream tool developed for allele-specific expression analysis.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R reference.fasta \
 *   -T ASEReadCounter \
 *   -o file_name.csv \
 *   -I input.bam \
 *   -sites sites.vcf \
 *   -U ALLOW_N_CIGAR_READS \
 *   [-minDepth 10] \
 *   [--minMappingQuality 10] \
 *   [--minBaseQuality 2] \
 *   [-drf DuplicateRead]
 * </pre>
 *
 * <h3>Note</h3>
 * <ul>
 *     <li>Like most GATK tools, this tools filters out duplicate reads by default. However, some ASE methods
 *     recommend including duplicate reads in the analysis, so the DuplicateRead filter can be disabled using the
 *     "-drf DuplicateRead" flag in the command-line.</li>
 * </ul>
 * <h3>Caveat</h3>
 * <ul>
 *     <li>This tool will only process biallelic SNP sites. If your callset contains multiallelic sites, they will be ignored.
 *     Optionally, you can subset your callset to just biallelic variants using e.g.
 *     <a href="org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php">SelectVariants</a>
 *     with the option "-restrictAllelesTo BIALLELIC".</li>
 * </ul>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@Downsample(by = DownsampleType.BY_SAMPLE, toCoverage = 10000)
//@DisabledReadFilters({DuplicateReadFilter.class})  //currently can be disabled using the command line argument -drf DuplicateRead
public class ASEReadCounter extends LocusWalker<String, Integer> {

    @Output
    public PrintStream out;

    @Input (fullName = "sitesVCFFile",shortName = "sites")
    public RodBinding<VariantContext> sites;

    /**
     * If this argument is enabled, loci with total depth lower than this threshold after all filters have been applied
     * will be skipped. This can be set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minDepthOfNonFilteredBase", shortName = "minDepth", doc = "Minimum number of bases that pass filters", required = false, minValue = -1, maxValue = Integer.MAX_VALUE)
    public int minDepthOfNonFilteredBases = -1;

    /**
     * If this argument is enabled, reads with mapping quality values lower than this threshold will not be counted.
     * This can be set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum read mapping quality", required = false, minValue = -1, maxValue = Integer.MAX_VALUE)
    public int minMappingQuality = 0;

    /**
     * If this argument is enabled, bases with quality scores lower than this threshold will not be counted.
     * This can be set to -1 by default to disable the evaluation and ignore this threshold.
     */
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum base quality", required = false, minValue = -1, maxValue = Byte.MAX_VALUE)
    public byte minBaseQuality = 0;

    /**
     * These options modify how the tool deals with overlapping read pairs. The default value is COUNT_FRAGMENTS_REQUIRE_SAME_BASE.
     */
    @Argument(fullName = "countOverlapReadsType", shortName = "overlap", doc = "Handling of overlapping reads from the same fragment", required = false)
    public CoverageUtils.CountPileupType countType = CoverageUtils.CountPileupType.COUNT_FRAGMENTS_REQUIRE_SAME_BASE;

    /**
     * Available options are csv, table, rtable. By default, the format is rtable (an r-readable table).
     */
    @Argument(fullName = "outputFormat", doc = "Format of the output file, can be CSV, TABLE, RTABLE", required = false)
    public OUTPUT_FORMAT outputFormat = OUTPUT_FORMAT.RTABLE;

    // Hiding these argument pending reevaluation (currently don't seem to work and aren't tested)
    /**
     * Consider a spanning deletion as contributing to coverage. Also enables deletion counts in per-base output.
     */
    @Hidden
    @Argument(fullName = "includeDeletions", shortName = "dels", doc = "Include information on deletions", required = false)
    public boolean includeDeletions = false;

    @Hidden
    @Argument(fullName = "ignoreDeletionSites", doc = "Ignore sites consisting only of deletions", required = false)
    public boolean ignoreDeletionSites = false;

    public String separator = "\t";

    public enum OUTPUT_FORMAT{
        TABLE,
        RTABLE,
        CSV
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // STANDARD WALKER METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    public boolean includeReadsWithDeletionAtLoci() { return includeDeletions && ! ignoreDeletionSites; }

    public void initialize() {

        // Check the output format
        boolean goodOutputFormat = false;
        for ( final OUTPUT_FORMAT f : OUTPUT_FORMAT.values()) {
            goodOutputFormat = goodOutputFormat || f.equals(outputFormat);
        }

        if ( ! goodOutputFormat ) {
            throw new IllegalArgumentException("Improper output format. Can be one of TABLE, RTABLE, CSV. Was "+outputFormat);
        }

        if ( outputFormat.equals(OUTPUT_FORMAT.CSV) ) {
            separator = ",";
        }
        final String header = "contig"+separator+"position"+separator+"variantID"+separator+"refAllele"+separator+"altAllele"+separator+"refCount"+separator+"altCount"+separator+"totalCount"+separator+"lowMAPQDepth"+separator+"lowBaseQDepth"+separator+"rawDepth"+separator+"otherBases"+separator+"improperPairs";
        out.println(header);

    }


    @Override
    public String map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        if ( tracker == null )
            return null;
        final String contig = context.getLocation().getContig();
        final long position = context.getPosition();

        if(position==2633936){
            System.out.println("stophere");
        }
        final char refAllele = (char)ref.getBase();

        final List<VariantContext> VCs =  tracker.getValues(sites, context.getLocation());
        if(VCs != null && VCs.size() > 1)
            throw new UserException("More then one variant context at position: "+contig+":"+position);
        if(VCs == null || VCs.isEmpty())
            return null;

        final VariantContext vc = VCs.get(0);
        if(!vc.isBiallelic()) {
            logger.warn("Ignoring site: cannot run ASE on non-biallelic sites: " + vc.toString());
            return null;
        }

        if ( vc.getNAlleles() == 1 || vc.getAlternateAllele(0).getBases().length == 0 )
            throw new UserException("The file of variant sites must contain heterozygous sites and cannot be a GVCF file containing <NON_REF> alleles.");

        final char altAllele = (char)vc.getAlternateAllele(0).getBases()[0];

        final String siteID = vc.getID();
        final ReadBackedPileup pileup = filterPileup(context.getBasePileup(), countType, includeReadsWithDeletionAtLoci());

        // count up the depths of all and QC+ bases
        return calculateLineForSite(pileup, contig, position, siteID, refAllele, altAllele);

    }

    protected ReadBackedPileup filterPileup(final ReadBackedPileup originalPileup, final CoverageUtils.CountPileupType countType, final boolean includeDeletions){

        ReadBackedPileup pileupWithDeletions;
        if(countType.equals(CoverageUtils.CountPileupType.COUNT_FRAGMENTS_REQUIRE_SAME_BASE))
            pileupWithDeletions = originalPileup.getOverlappingFragmentFilteredPileup(true,true);
        else if(countType.equals(CoverageUtils.CountPileupType.COUNT_READS))
            pileupWithDeletions = originalPileup;
        else if(countType.equals(CoverageUtils.CountPileupType.COUNT_FRAGMENTS))
            pileupWithDeletions = originalPileup.getOverlappingFragmentFilteredPileup(false,true);
        else
            throw new UserException("Must use valid CountPileupType");

        return includeDeletions ? pileupWithDeletions: pileupWithDeletions.getPileupWithoutDeletions();

    }

    protected String calculateLineForSite(final ReadBackedPileup pileup, final String contig, final long position, final String siteID, final char refAllele, final char altAllele){

        int rawDepth = 0, lowBaseQDepth = 0, lowMAPQDepth = 0, refCount = 0, altCount = 0, totalNonFilteredCount = 0, otherBasesCount = 0, improperPairsCount = 0 ;

        for (final PileupElement base : pileup) {
            rawDepth++;

            if (base.getRead().getReadPairedFlag() && (base.getRead().getMateUnmappedFlag() || !base.getRead().getProperPairFlag())){
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

            if(base.getBase() == refAllele)
                refCount++;
            else if(base.getBase() == altAllele)
                altCount++;
            else {
                otherBasesCount++;
                continue;
            }
            totalNonFilteredCount++;
        }

        if(totalNonFilteredCount < minDepthOfNonFilteredBases)
            return null;

        return contig +separator+
                position +separator+
                siteID +separator+
                refAllele +separator+
                altAllele +separator+
                refCount +separator+
                altCount +separator+
                totalNonFilteredCount +separator+
                lowMAPQDepth +separator+
                lowBaseQDepth +separator+
                rawDepth +separator+
                otherBasesCount +separator+
                improperPairsCount;
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(String results, Integer sum) {
        if(results!= null)
            out.println(results);
        return ++sum;
    }

    @Override
    public void onTraversalDone(Integer sum) {
        logger.info("Done processing "+sum+" loci");
        out.close();
    }





}
