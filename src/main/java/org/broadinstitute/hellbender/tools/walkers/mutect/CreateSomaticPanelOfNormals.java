package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;

/**
 * Create a panel of normals (PoN) containing germline and artifactual sites for use with Mutect2.
 *
 * <p>
 *     The tool takes multiple normal sample callsets produced by {@link Mutect2}'s tumor-only mode and collates sites present in two or more samples
 *     into a sites-only VCF. The PoN captures common artifactual and germline variant sites.
 *     Mutect2 then uses the PoN to filter variants at the site-level.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 * <h3>Example workflow</h3>
 *
 * <h4>Step 1. Run Mutect2 in tumor-only mode for each normal sample.</h4>
 * <pre>
 * gatk Mutect2 \
 *   -R reference.fa \
 *   -I normal1.bam \
 *   -tumor normal1_sample_name \
 *   --germline-resource af-only-gnomad.vcf.gz \
 *   -O normal1_for_pon.vcf.gz
 * </pre>
 *
 * <h4>Step 2. Create a file ending with .args or .list extension with the paths to the VCFs from step 1, one per line.</h4>
 * <p>This approach is optional. Other extensions will error the run. </p>
 *
 * <pre>
 *     normal1_for_pon.vcf.gz
 *     normal2_for_pon.vcf.gz
 *     normal3_for_pon.vcf.gz
 * </pre>
 *
 * <h4>Step 3. Combine the normal calls using CreateSomaticPanelOfNormals.</h4>
 *
 * <pre>
 * gatk CreateSomaticPanelOfNormals \
 *   -vcfs normals_for_pon_vcf.args \
 *   -O pon.vcf.gz
 * </pre>
 *
 * <p>The tool also accepts multiple .args files. Pass each in with the -vcfs option.
 * Alternatively, provide each normal's VCF as separate arguments.</p>
 * <pre>
 * gatk CreateSomaticPanelOfNormals \
 *   -vcfs normal1_for_pon_vcf.gz \
 *   -vcfs normal2_for_pon_vcf.gz \
 *   -vcfs normal3_for_pon_vcf.gz \
 *   -O pon.vcf.gz
 * </pre>
 *
 * <p>The resulting VCF will be an eight-column sites-only VCF lacking annotations.</p>
 *
 * <p>By default the tool fails if multiple vcfs have the same sample name, but the --duplicate-sample-strategy argument can be changed to
 *  ALLOW_ALL to allow duplicates or CHOOSE_FIRST to use only the first vcf with a given sample name.</p>
 *
 */
@CommandLineProgramProperties(
        summary = "Make a panel of normals (PoN) for use with Mutect2",
        oneLineSummary = "Make a panel of normals for use with Mutect2",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class CreateSomaticPanelOfNormals extends CommandLineProgram {

    public static final String INPUT_VCFS_LIST_LONG_NAME = "vcfs";
    public static final String INPUT_VCFS_LIST_SHORT_NAME = "vcfs";

    public static final String MIN_SAMPLE_COUNT_LONG_NAME = "min-sample-count";

    public static final int DEFAULT_MIN_SAMPLE_COUNT = 2;

    public static final String DUPLICATE_SAMPLE_STRATEGY_LONG_NAME = "duplicate-sample-strategy";

    public enum DuplicateSampleStrategy {
        THROW_ERROR, CHOOSE_FIRST, ALLOW_ALL
    }

    /**
     * The VCFs can be input as either one or more .args file(s) containing one VCF per line, or VCFs can be
     * specified explicitly on the command line.
     */
    @Argument(fullName = INPUT_VCFS_LIST_LONG_NAME,
            shortName = INPUT_VCFS_LIST_SHORT_NAME,
            doc="VCFs for samples to include. May be specified either one at a time, or as one or more .args file containing multiple VCFs, one per line.", optional = false)
    private Set<File> vcfs = new LinkedHashSet<>(0);

    /**
     * Number of samples containing a variant site required to include it in the panel of normals.
     */
    @Argument(fullName = MIN_SAMPLE_COUNT_LONG_NAME,
            doc="Number of samples containing a variant site required to include it in the panel of normals.", optional = true)
    private int minSampleCount = DEFAULT_MIN_SAMPLE_COUNT;

    /**
     * How to handle duplicate samples: THROW_ERROR to fail, CHOOSE_FIRST to use the first vcf with each sample name, ALLOW_ALL to use all samples regardless of duplicate sample names."
     */
    @Argument(fullName = DUPLICATE_SAMPLE_STRATEGY_LONG_NAME,
            doc="How to handle duplicate samples: THROW_ERROR to fail, CHOOSE_FIRST to use the first vcf with each sample name, ALLOW_ALL to use all samples regardless of duplicate sample names.", optional = false)
    private DuplicateSampleStrategy duplicateSampleStrategy = DuplicateSampleStrategy.THROW_ERROR;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Output vcf", optional = false)
    private File outputVcf = null;

    public Object doWork() {
        final List<File> inputVcfs = new ArrayList<>(vcfs);
        if (!inputVcfs.stream().map(File::getAbsolutePath).allMatch(path -> path.endsWith(".vcf") || path.endsWith(".vcf.gz") || path.endsWith(".args") )) {
            logger.warn("Some input files don't seem to be .vcf or .args files.  Make sure that any input vcfs list end in .args.");
        }
        final Collection<CloseableIterator<VariantContext>> iterators = new ArrayList<>(inputVcfs.size());
        final Collection<VCFHeader> headers = new HashSet<>(inputVcfs.size());
        final VCFHeader headerOfFirstVcf = new VCFFileReader(inputVcfs.get(0), false).getFileHeader();
        final SAMSequenceDictionary sequenceDictionary = headerOfFirstVcf.getSequenceDictionary();
        final VariantContextComparator comparator = headerOfFirstVcf.getVCFRecordComparator();

        final Set<String> samples = new HashSet<>();
        for (final File vcf : inputVcfs) {
            final VCFFileReader reader = new VCFFileReader(vcf, false);

            final VCFHeader header = reader.getFileHeader();
            final String sample = header.getGenotypeSamples().get(0);
            if (duplicateSampleStrategy == DuplicateSampleStrategy.THROW_ERROR && samples.contains(sample)) {
                throw new UserException.BadInput(String.format("Duplicate sample name %s found in multiple input vcfs, the second one being %s.  Consider changing the %s argument", sample, vcf.getAbsolutePath(), DUPLICATE_SAMPLE_STRATEGY_LONG_NAME));
            } else if (duplicateSampleStrategy == DuplicateSampleStrategy.CHOOSE_FIRST && samples.contains(sample)) {
                logger.info(String.format("Skipping input vcf %s because a different vcf with the same sample %s has already been seen", vcf.getAbsolutePath(), sample));
                continue;
            } else {
                Utils.validateArg(comparator.isCompatible(header.getContigLines()), () -> vcf.getAbsolutePath() + " has incompatible contigs.");
                headers.add(header);
                iterators.add(reader.iterator());
                samples.add(sample);
            }
        }

        final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(outputVcf, sequenceDictionary, false, Options.INDEX_ON_THE_FLY);
        writer.writeHeader(new VCFHeader(VCFUtils.smartMergeHeaders(headers, false)));

        final MergingIterator<VariantContext> mergingIterator = new MergingIterator<>(comparator, iterators);
        SimpleInterval currentPosition = new SimpleInterval("FAKE", 1, 1);
        final List<VariantContext> variantsAtThisPosition = new ArrayList<>(20);
        while (mergingIterator.hasNext()) {
            final VariantContext vc = mergingIterator.next();
            if (!currentPosition.overlaps(vc)) {
                processVariantsAtSamePosition(variantsAtThisPosition, writer);
                variantsAtThisPosition.clear();
                currentPosition = new SimpleInterval(vc.getContig(), vc.getStart(), vc.getStart());
            }
            variantsAtThisPosition.add(vc);
        }
        processVariantsAtSamePosition(variantsAtThisPosition, writer);
        mergingIterator.close();
        writer.close();

        return "SUCCESS";
    }

    private void processVariantsAtSamePosition(final List<VariantContext> variants, final VariantContextWriter writer) {
        if (variants.size() >= minSampleCount){
            final VariantContext mergedVc = AssemblyBasedCallerUtils.makeMergedVariantContext(variants);
            final VariantContext outputVc = new VariantContextBuilder()
                    .source(mergedVc.getSource())
                    .loc(mergedVc.getContig(), mergedVc.getStart(), mergedVc.getEnd())
                    .alleles(mergedVc.getAlleles())
                    .make();
            writer.add(outputVc);
        }
    }
}
