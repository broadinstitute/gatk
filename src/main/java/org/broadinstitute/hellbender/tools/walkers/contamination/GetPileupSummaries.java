package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Predicate;

/**
 * <p>Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with {@link CalculateContamination}.</p>
 *
 * <p>
 * The tool requires a <i>common</i> germline variant sites VCF, e.g. derived from the gnomAD resource, with population allele frequencies (AF) in the INFO field.
 * This resource must contain only biallelic SNPs and can be an eight-column sites-only VCF.
 * The tool ignores the filter status of the variant calls in this germline resource.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 *     In particular, the mutect_resources.wdl script prepares a suitable resource from a larger dataset. An example excerpt is shown.
 * </p>
 *
 * <pre>
 * #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
 * chr6	29942512	.	G	C	2974860	VQSRTrancheSNP99.80to99.90	AF=0.063
 * chr6	29942517	.	C	A	2975860	VQSRTrancheSNP99.80to99.90	AF=0.062
 * chr6	29942525	.	G	C	2975600	VQSRTrancheSNP99.60to99.80	AF=0.063
 * chr6	29942547	rs114945359	G	C	15667700	PASS	AF=0.077
 * </pre>
 *
 * <h3>Usage examples</h3>
 *
 * <pre>
 * gatk GetPileupSummaries \
 *   -I tumor.bam \
 *   -V common_biallelic.vcf.gz \
 *   -L common_biallelic.vcf.gz \
 *   -O pileups.table
 * </pre>
 *
 * <pre>
 * gatk GetPileupSummaries \
 *   -I normal.bam \
 *   -V common_biallelic.vcf.gz \
 *   -L common_biallelic.vcf.gz \
 *   -O pileups.table
 * </pre>
 *
 * Although the sites (-L) and variants (-V) resources will often be identical, this need not be the case.  For example,
 * <pre>
 * gatk GetPileupSummaries \
 *   -I normal.bam \
 *   -V gnomad.vcf.gz \
 *   -L common_snps.interval_list \
 *   -O pileups.table
 * </pre>
 * attempts to get pileups at a list of common snps and emits output for those sites that are present in gnomAD, using the
 * allele frequencies from gnomAD.  Note that the sites may be a subset of the variants, the variants may be a subset of the sites,
 * or they may overlap partially.  In all cases pileup summaries are emitted for the overlap and nowhere else.  The most common use
 * case in which sites and variants differ is when the variants resources is a large file and the sites is an interval list subset from that file.
 *
 * <p>
 * GetPileupSummaries tabulates results into six columns as shown below.
 * The alt_count and allele_frequency correspond to the ALT allele in the germline resource.
 * </p>
 *
 * <pre>
 * contig	position	ref_count	alt_count	other_alt_count	allele_frequency
 * chr6	29942512	9	0	0	0.063
 * chr6	29942517	13	1	0	0.062
 * chr6	29942525	13	7	0	0.063
 * chr6	29942547	36	0	0	0.077
 * </pre>
 *
 * <p>
 * Note the default maximum population AF ({@code --maximum-population-allele-frequency} or {@code -max-af})
 * is set to 0.2, which limits the sites the tool considers to those in the variants resource file that have
 * AF of 0.2 or less. Likewise, the default minimum population AF ({@code --minimum-population-allele-frequency}
 * or {@code -min-af}) is set to 0.01, which limits the sites the tool considers to those in the variants resource
 * file that have AF of 0.01 or more.
 * </p>
 *
 */
@CommandLineProgramProperties(
        summary = "Tabulates pileup metrics for inferring contamination",
        oneLineSummary = "Tabulates pileup metrics for inferring contamination",
        programGroup = CoverageAnalysisProgramGroup.class)
@DocumentedFeature
public class GetPileupSummaries extends LocusWalker {

    public static final String MAX_SITE_AF_LONG_NAME = "maximum-population-allele-frequency";
    public static final String MIN_SITE_AF_LONG_NAME = "minimum-population-allele-frequency";
    public static final String MAX_SITE_AF_SHORT_NAME = "max-af";
    public static final String MIN_SITE_AF_SHORT_NAME = "min-af";
    public static final String MIN_MAPPING_QUALITY_LONG_NAME = "min-mapping-quality";
    public static final String MIN_MAPPING_QUALITY_SHORT_NAME = "mmq";

    private static final double DEFAULT_MIN_POPULATION_AF = 0.01;
    private static final double DEFAULT_MAX_POPULATION_AF = 0.2;
    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 50;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table", optional=false)
    private File outputTable;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants and allele frequencies")
    public FeatureInput<VariantContext> variants;

    @Argument(fullName = MIN_SITE_AF_LONG_NAME,
            shortName = MIN_SITE_AF_SHORT_NAME,
            doc = "Minimum population allele frequency of sites to consider.  A low value increases accuracy at the expense of speed.", optional = true)
    private double minPopulationAlleleFrequency = DEFAULT_MIN_POPULATION_AF;

    @Argument(fullName = MAX_SITE_AF_LONG_NAME,
            shortName = MAX_SITE_AF_SHORT_NAME,
            doc = "Maximum population allele frequency of sites to consider.", optional = true)
    private double maxPopulationAlleleFrequency = DEFAULT_MAX_POPULATION_AF;

    @Argument(fullName = MIN_MAPPING_QUALITY_LONG_NAME, shortName = MIN_MAPPING_QUALITY_SHORT_NAME, doc = "Minimum read mapping quality", optional = true)
    private int minMappingQuality = DEFAULT_MINIMUM_MAPPING_QUALITY;

    // Filter pileup elements that appear within this many bases from the ends of reads
    private static int MAX_REQUIRED_DISTANCE_FROM_END = 5;

    public static final String POST_FILTER_PILEUP_BAM = "filtered-pileup-bam";
    @Argument(fullName = POST_FILTER_PILEUP_BAM, optional = true)
    private GATKPath filteredPileupBam = null;
    GATKReadWriter bamWriter;

    public static final String LEGIT_READ_BAM_ARG_NAME = "legit-reads-bam";
    @Argument(fullName = LEGIT_READ_BAM_ARG_NAME, optional = true)
    private GATKPath legitReadsBam = null;
    GATKReadWriter legitBamWriter;

    public static final String BAD_READ_BAM_ARG_NAME = "bad-reads-bam";
    @Argument(fullName = BAD_READ_BAM_ARG_NAME, optional = true)
    private GATKPath badReadsBam = null;
    GATKReadWriter badBamWriter;

    @Argument(fullName = "counts", optional = true)
    private GATKPath countsFile = null;

    private final List<PileupSummary> pileupSummaries = new ArrayList<>();

    private boolean sawVariantsWithoutAlleleFrequency = false;
    private boolean sawVariantsWithAlleleFrequency = false;

    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine;
    private SampleList samplesList;

    private static int INITIAL_SET_SIZE = 10_000;

    private Set<String> readNames = new HashSet<>(INITIAL_SET_SIZE);
    private String currentContig = "";

    int numDeletions = 0;
    int numRefAllele = 0;
    int numShouldBeAltButRef = 0;
    int numShouldBeAltButDel = 0;

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    public boolean requiresReference() {
        return false;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public boolean requiresFeatures() {
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>();
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
        filters.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.PRIMARY_LINE);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE);
        filters.add(ReadFilterLibrary.GOOD_CIGAR);
        filters.add(new WellformedReadFilter());
        return filters;
    }

    @Override
    public void onTraversalStart() {
        final boolean alleleFrequencyInHeader = ((VCFHeader) getHeaderForFeatures(variants)).getInfoHeaderLines().stream()
                .anyMatch(line -> line.getID().equals(VCFConstants.ALLELE_FREQUENCY_KEY));
        if (!alleleFrequencyInHeader) {
            throw new UserException.BadInput("Population vcf does not have an allele frequency (AF) info field in its header.");
        }

        if (filteredPileupBam != null){
            bamWriter = createSAMWriter(filteredPileupBam, false);
        }

        if (legitReadsBam != null){
            legitBamWriter = createSAMWriter(legitReadsBam, false);
        }

        if (badReadsBam != null){
            badBamWriter = createSAMWriter(badReadsBam, false);
        }

        M2ArgumentCollection MTAC = new M2ArgumentCollection();
        MTAC.likelihoodArgs.pairHMM = PairHMM.Implementation.ORIGINAL;
        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(MTAC.likelihoodArgs, true);
        samplesList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(getHeaderForReads())));

    }

    @Override
    public void apply(AlignmentContext alignmentContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final List<VariantContext> vcs = featureContext.getValues(variants);
        if (vcs.isEmpty()) {
            return;
        }
        final VariantContext vc = vcs.get(0); // sato: this is the variant in exac

        if (vc.isBiallelic() && vc.isSNP() && alleleFrequencyInRange(vc)) {
            final ReadPileup rawPileup = alignmentContext.getBasePileup();
            final Predicate<PileupElement> mqFilter = pe -> pe.getRead().getMappingQuality() >= minMappingQuality;
            final Predicate<PileupElement> endsOfReadsFilter = pe -> Math.min(pe.getOffset(), pe.getRead().getLength() - pe.getOffset()) >=  MAX_REQUIRED_DISTANCE_FROM_END;

            ReadPileup pileup = alignmentContext.getBasePileup()
                    .makeFilteredPileup(mqFilter.and(endsOfReadsFilter));
            final int numReadsRemoved = rawPileup.size() - pileup.size();

            final byte altBase = vc.getAlternateAllele(0).getBases()[0];
            final byte refBase = vc.getReference().getBases()[0];
            for (PileupElement elem : pileup) {
                if (elem.getBase() == refBase || elem.isDeletion()) {
                    // Get a new pileup element, using the realigned read.
                    // Here I should either directly align to the alt haplotype. Or I should align to the alt (best) haplotype, then to ref, see if that alignment is better.
                    // Or---should I do pairHMM? Realign to the reference. Start small, baby.
                    final int padding = 5;
                    referenceContext.setWindow(padding, padding);
                    final byte[] refBases = Arrays.copyOf(referenceContext.getBases(), referenceContext.getBases().length);
                    final byte[] altBases = Arrays.copyOf(referenceContext.getBases(), referenceContext.getBases().length);
                    final byte[] readBases = Arrays.copyOfRange(elem.getRead().getBases(), Math.max(0, elem.getOffset() - padding), Math.min(elem.getOffset() + padding + 1, elem.getRead().getLength()));
                    final int variantOffsetWrtRef = padding;
                    if (refBases[variantOffsetWrtRef] != refBase) {
                        throw new UserException("refBase doesn't match: " + refBases[variantOffsetWrtRef] + ", " + refBase);
                    }

                    altBases[variantOffsetWrtRef] = altBase;
                    final int haplotypeSize = refBases.length;
                    final Cigar refCigar = new Cigar(Arrays.asList(new CigarElement(haplotypeSize, CigarOperator.M)));
                    final Haplotype altHaplotype = new Haplotype(altBases, false, 0, refCigar);
                    final Haplotype refHaplotype = new Haplotype(refBases, true, 0, refCigar);
                    final int pileupPosition = pileup.getLocation().getStart();
                    SimpleInterval haplotypeLocation = new SimpleInterval(pileup.getLocation().getContig(), pileupPosition - padding, pileupPosition + padding);
                    altHaplotype.setGenomeLocation(haplotypeLocation);
                    refHaplotype.setGenomeLocation(haplotypeLocation);

                    // PairHMM
                    int globalCoordinate = pileup.getLocation().getStart();
                    final GATKRead clippedRead = ReadClipper.hardClipToRegion(elem.getRead(), globalCoordinate - padding, globalCoordinate + padding);
                    final Map<String, List<GATKRead>> reads = new HashMap<>();
                    reads.put(samplesList.getSample(0), Arrays.asList(clippedRead));
                    AssemblyResultSet ars = new AssemblyResultSet();
                    ars.add(refHaplotype);
                    ars.add(altHaplotype);
                    final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(ars, samplesList, reads);
                    final double refLikelihood = readLikelihoods.sampleMatrix(0).get(0, 0);
                    final double altLikelihood = readLikelihoods.sampleMatrix(0).get(1, 0);

                    String readName = elem.getRead().getName();
                    if (elem.isDeletion()) {
                        numDeletions++;
                        if (refLikelihood < altLikelihood) {
                            numShouldBeAltButDel++;
                            if (badBamWriter != null){
                                badBamWriter.addRead(elem.getRead());
                            }
                            // This is a no-op...I should realign the read to the alt haplotype!
                            // final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.FASTEST_AVAILABLE);
                            // final GATKRead realignedRead = AlignmentUtils.createReadAlignedToRef(elem.getRead(), altHaplotype, refHaplotype, refHaplotype.getStartPosition(), true, aligner);
                        } else {
                            if (legitBamWriter != null){
                                legitBamWriter.addRead(elem.getRead());
                            }
                        }
                    } else if (elem.getBase() == refBase) {
                        numRefAllele++;
                        if (refLikelihood < altLikelihood) {
                            numShouldBeAltButRef++;
                            if (badBamWriter != null){
                                badBamWriter.addRead(elem.getRead());
                            }

                            int numElemBefore = pileup.getReads().size();
                            pileup = pileup.makeFilteredPileup(pe -> !pe.getRead().getName().equals(elem.getRead().getName()));
                            int numElemAfter = pileup.getReads().size();
                            logger.info("Removed " + elem.getRead().getName() + " at " + elem.getRead().getContig() + ": " + elem.getRead().getStart());
                        } else {
                            if (legitBamWriter != null){
                                legitBamWriter.addRead(elem.getRead());
                            }
                        }
                    }
                }
            }

            pileupSummaries.add(new PileupSummary(vc, pileup));

            final List<GATKRead> pileupReads = pileup.getReads();
            final List<GATKRead> sortedPileupReads = new ArrayList<>();
            pileup.sortedIterator().forEachRemaining(pe -> sortedPileupReads.add(pe.getRead()));

            if (!currentContig.equals(pileup.getLocation().getContig())) {
                currentContig = pileup.getLocation().getContig();
                readNames.clear();
            }

            if (bamWriter != null && pileup.size() > 0) {
                // We must check that the read hasn't already been added as an output---with a LocusIterator,
                // it's possible that we visit reads out of order.
                Iterator<PileupElement> sortedIterator = pileup.sortedIterator();
                while (sortedIterator.hasNext()) {
                    PileupElement e = sortedIterator.next();
                    if (readNames.contains(e.getRead().getName())) {
                        continue;
                    } else {
                        readNames.add(e.getRead().getName());
                        bamWriter.addRead(e.getRead());
                    }

                }
            }
        }
    }

    private void filterPileupElementBeta(final PileupElement pe){
        final GATKRead read = pe.getRead();
        // AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype();
        // AlignmentUtils.createReadAlignedToRef();

    }

    @Override
    public void closeTool(){
        likelihoodCalculationEngine.close();
        if (bamWriter != null){
            try {
                bamWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        if (legitBamWriter != null){
            try {
                legitBamWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        if (badBamWriter != null){
            try {
                badBamWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        if (sawVariantsWithoutAlleleFrequency && !sawVariantsWithAlleleFrequency) {
            throw new UserException.BadInput("No variants in population vcf had an allele frequency (AF) field.");
        }
        final String sampleName = ReadUtils.getSamplesFromHeader(getHeaderForReads()).stream().findFirst().get();
        PileupSummary.writeToFile(sampleName, pileupSummaries, outputTable);


        // For good measure
        System.out.println("Total_ref_count_in_hom_alt," + numRefAllele);
        System.out.println("Ref_should_be_alt," + numShouldBeAltButRef);
        System.out.println("Total_deletion_count_in_hom_alt," + numDeletions);
        System.out.println("Deletion_should_be_alt," + numShouldBeAltButDel);

        if (countsFile != null){
            try (PrintWriter memoWriter = new PrintWriter(countsFile.toString())){
                memoWriter.println("Total_ref_count_in_hom_alt," + numRefAllele);
                memoWriter.println("Ref_should_be_alt," + numShouldBeAltButRef);
                memoWriter.println("Total_deletion_count_in_hom_alt," + numDeletions);
                memoWriter.println("Deletion_should_be_alt," + numShouldBeAltButDel);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }


        return "SUCCESS";
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
}
