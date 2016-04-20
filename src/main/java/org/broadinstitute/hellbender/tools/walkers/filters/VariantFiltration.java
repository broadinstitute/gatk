package org.broadinstitute.hellbender.tools.walkers.filters;

import com.google.common.collect.Sets;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.VariantContextUtils.JexlVCMatchExp;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;

import java.util.*;
import java.util.stream.Collectors;

import static java.util.Collections.singleton;


/**
 * Filter variant calls based on INFO and FORMAT annotations
 *
 * <p>
 * This tool is designed for hard-filtering variant calls based on certain criteria.
 * Records are hard-filtered by changing the value in the FILTER field to something other than PASS. Filtered records
 * will be preserved in the output unless their removal is requested in the command line. </p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant set to filter.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A filtered VCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T VariantFiltration \
 *   -R reference.fasta \
 *   -o output.vcf \
 *   --variant input.vcf \
 *   --filterExpression "AB < 0.2 || MQ0 > 50" \
 *   --filterName "Nov09filters" \
 *   --mask mask.vcf \
 *   --maskName InDel
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Filter variant calls based on INFO and FORMAT annotations.",
        oneLineSummary = "Hard-filter variants VCF (mark them as FILTER)",
        programGroup = VariantProgramGroup.class
)
public final class VariantFiltration extends VariantWalker {

    /**
     * Any variant which overlaps entries from the provided mask file will be filtered. If the user wants logic to be reversed,
     * i.e. filter variants that do not overlap with provided mask, then argument -filterNotInMask can be used.
     * Note that it is up to the user to adapt the name of the mask to make it clear that the reverse logic was used
     * (e.g. if masking against Hapmap, use -maskName=hapmap for the normal masking and -maskName=not_hapmap for the reverse masking).
     */
    @Argument(fullName="mask", shortName="mask", doc="Input mask", optional=true)
    public FeatureInput<Feature> mask;

    @Argument(doc="File to which variants should be written", fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, optional = false)
    public String out = null;

    /**
     * VariantFiltration accepts any number of JEXL expressions (so you can have two named filters by using
     * --filterName One --filterExpression "X < 1" --filterName Two --filterExpression "X > 2").
     */
    @Argument(fullName="filterExpression", shortName="filter", doc="One or more expression used with INFO fields to filter", optional=true)
    public List<String> filterExpressions = new ArrayList<>();

    /**
     * This name is put in the FILTER field for variants that get filtered.  Note that there must be a 1-to-1 mapping between filter expressions and filter names.
     */
    @Argument(fullName="filterName", shortName="filterName", doc="Names to use for the list of filters", optional=true)
    public List<String> filterNames = new ArrayList<>();

    /**
     * Similar to the INFO field based expressions, but used on the FORMAT (genotype) fields instead.
     * VariantFiltration will add the sample-level FT tag to the FORMAT field of filtered samples (this does not affect the record's FILTER tag).
     * One can filter normally based on most fields (e.g. "GQ < 5.0"), but the GT (genotype) field is an exception. We have put in convenience
     * methods so that one can now filter out hets ("isHet == 1"), refs ("isHomRef == 1"), or homs ("isHomVar == 1"). Also available are
     * expressions isCalled, isNoCall, isMixed, and isAvailable, in accordance with the methods of the Genotype object.
     */
    @Argument(fullName="genotypeFilterExpression", shortName="G_filter", doc="One or more expression used with FORMAT (sample/genotype-level) fields to filter (see documentation guide for more info)", optional=true)
    public List<String> genotypeFilterExpressions = new ArrayList<>();

    /**
     * Similar to the INFO field based expressions, but used on the FORMAT (genotype) fields instead.
     */
    @Argument(fullName="genotypeFilterName", shortName="G_filterName", doc="Names to use for the list of sample/genotype filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered", optional=true)
    public List<String> genotypeFilterNames = new ArrayList<>();

    /**
     * Works together with the --clusterWindowSize argument.
     */
    @Argument(fullName="clusterSize", shortName="cluster", doc="The number of SNPs which make up a cluster. Must be at least 2", optional=true)
    public Integer clusterSize = 3;

    /**
     * Works together with the --clusterSize argument.  To disable the clustered SNP filter, set this value to less than 1.
     */
    @Argument(fullName="clusterWindowSize", shortName="window", doc="The window size (in bases) in which to evaluate clustered SNPs", optional=true)
    public Integer clusterWindow = 0;

    @Argument(fullName="maskExtension", shortName="maskExtend", doc="How many bases beyond records from a provided 'mask' should variants be filtered", optional=true)
    public Integer maskExtension = 0;

    /**
     * When using the -mask argument, the maskName will be annotated in the variant record.
     * Note that when using the -filterNotInMask argument to reverse the masking logic,
     * it is up to the user to adapt the name of the mask to make it clear that the reverse logic was used
     * (e.g. if masking against Hapmap, use -maskName=hapmap for the normal masking and -maskName=not_hapmap for the reverse masking).
     */
    @Argument(fullName="maskName", shortName="maskName", doc="The text to put in the FILTER field if a 'mask' is provided and overlaps with a variant call", optional=true)
    public String maskName = "Mask";

    /**
     * By default, if the -mask argument is used, any variant falling in a mask will be filtered.
     * If this argument is used, logic is reversed, and variants falling outside a given mask will be filtered.
     * Use case is, for example, if we have an interval list or BED file with "good" sites.
     * Note that it is up to the user to adapt the name of the mask to make it clear that the reverse logic was used
     * (e.g. if masking against Hapmap, use -maskName=hapmap for the normal masking and -maskName=not_hapmap for the reverse masking).
     */
    @Argument(fullName="filterNotInMask", shortName="filterNotInMask", doc="Filter records NOT in given input mask.", optional=true)
    public boolean filterRecordsNotInMask = false;

    /**
     * By default, if JEXL cannot evaluate your expression for a particular record because one of the annotations is not present, the whole expression evaluates as PASSing.
     * Use this argument to have it evaluate as failing filters instead for these cases.
     */
    @Argument(fullName="missingValuesInExpressionsShouldEvaluateAsFailing", doc="When evaluating the JEXL expressions, missing values should be considered failing the expression", optional=true)
    public Boolean failMissingValues = false;

    /**
     * Invalidate previous filters applied to the VariantContext, applying only the filters here
     */
    @Argument(fullName="invalidatePreviousFilters",doc="Remove previous filters applied to the VCF",optional=true)
    boolean invalidatePreviousFilters = false;

    /**
     * Invert the selection criteria for --filterExpression
     */
    @Argument(fullName="invertFilterExpression", shortName="invfilter", doc="Invert the selection criteria for --filterExpression", optional=true)
    public boolean invertFilterExpression = false;

    /**
     * Invert the selection criteria for --genotypeFilterExpression
     */
    @Argument(fullName="invertGenotypeFilterExpression", shortName="invG_filter", doc="Invert the selection criteria for --genotypeFilterExpression", optional=true)
    public boolean invertGenotypeFilterExpression = false;

    /**
     * If this argument is provided, set filtered genotypes to no-call (./.).
     */
    @Argument(fullName="setFilteredGtToNocall", optional=true, doc="Set filtered genotypes to no-call")
    public boolean setFilteredGenotypesToNocall = false;

    // JEXL expressions for the filters
    private List<JexlVCMatchExp> filterExps;
    private List<JexlVCMatchExp> genotypeFilterExps;

    public static final String CLUSTERED_SNP_FILTER_NAME = "SnpCluster";

    private VariantContextWriter writer;

    private final List<Allele> diploidNoCallAlleles = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    /**
     * Invert logic if specified
     *
     * @param logic boolean logical operation value
     * @param invert whether to invert logic
     * @return invert logic if invert flag is true, otherwise leave the logic
     */
    private static boolean invertLogic(final boolean logic, final boolean invert){
        return invert ? !logic : logic;
    }

    /**
     * Prepend inverse phrase to description if --invertFilterExpression
     *
     * @param description the description
     * @return the description with inverse prepended if --invert_filter_expression
     */
    private String possiblyInvertFilterExpression( final String description ){
        return invertFilterExpression ? "Inverse of: " + description : description;
    }

    private void initializeVcfWriter() {
        //TODO remove hardwiring to output VCFs
        writer = new VariantContextWriterBuilder().setOutputFile(out).setOutputFileType(VariantContextWriterBuilder.OutputType.VCF).unsetOption(Options.INDEX_ON_THE_FLY).build();

        // setup the header fields
        final Set<VCFHeaderLine> hInfo = new HashSet<>();
        hInfo.addAll(getHeaderForVariants().getMetaDataInInputOrder());

        if ( clusterWindow > 0 ) {
            hInfo.add(new VCFFilterHeaderLine(CLUSTERED_SNP_FILTER_NAME, "SNPs found in clusters"));
        }

        if ( !genotypeFilterExps.isEmpty() ) {
            hInfo.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
        }

        try {
            for ( final JexlVCMatchExp exp : filterExps ) {
                hInfo.add(new VCFFilterHeaderLine(exp.name, possiblyInvertFilterExpression(exp.exp.toString())));
            }
            for ( final JexlVCMatchExp exp : genotypeFilterExps ) {
                hInfo.add(new VCFFilterHeaderLine(exp.name, possiblyInvertFilterExpression(exp.exp.toString())));
            }

            if ( mask != null ) {
                if (filterRecordsNotInMask) {
                    hInfo.add(new VCFFilterHeaderLine(maskName, "Doesn't overlap a user-input mask"));
                } else {
                    hInfo.add(new VCFFilterHeaderLine(maskName, "Overlaps a user-input mask"));
                }
            }
        } catch (final IllegalArgumentException e) {
            throw new UserException.BadInput(e.getMessage());
        }

        writer.writeHeader(new VCFHeader(hInfo, getHeaderForVariants().getGenotypeSamples()));
    }

    @Override
    public void onTraversalStart() {
        if (clusterSize <= 1){
            throw new UserException.BadArgumentValue("clusterSize", "values lower than 2 are not allowed");
        }
        if ( maskExtension < 0 ) {
            throw new UserException.BadArgumentValue("maskExtension", "negative values are not allowed");
        }

        if (filterRecordsNotInMask && mask == null) {
            throw new UserException.BadArgumentValue("filterNotInMask", "argument not allowed if mask argument is not provided");
        }
        filterExps = VariantContextUtils.initializeMatchExps(filterNames, filterExpressions);
        genotypeFilterExps = VariantContextUtils.initializeMatchExps(genotypeFilterNames, genotypeFilterExpressions);

        VariantContextUtils.engine.get().setSilent(true);

        initializeVcfWriter();
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {

        final VariantContext vc1 = invalidatePreviousFilters ? (new VariantContextBuilder(variant)).unfiltered().make() : variant;
        final VariantContext vc = isMaskFilterPresent(vc1) ? vc1: addMaskIfCoversVariant(vc1, featureContext);

        filter(vc, featureContext);
    }

    /**
     * Add mask to variant context filters if it covers its location
     * @return VariantContext with the mask added if the VariantContext is within the extended mask area
     */
    private VariantContext addMaskIfCoversVariant(final VariantContext vc, final FeatureContext featureContext) {
        final List<Feature> maskVariants = featureContext.getValues(mask, maskExtension, maskExtension);

        final boolean variantsMasked = maskVariants.isEmpty() == filterRecordsNotInMask;
        if (variantsMasked) {
            final Set<String> oldFiltersPlusNewOne = Sets.union(vc.getFilters(), singleton(maskName));
            return new VariantContextBuilder(vc).filters(oldFiltersPlusNewOne).make();
        } else {
            return vc;
        }
    }

    private boolean isMaskFilterPresent(final VariantContext vc) {
        return vc.getFilters() != null && vc.getFilters().contains(maskName);
    }

    private void filter(final VariantContext vc, final FeatureContext featureContext) {
        final VariantContextBuilder builder = new VariantContextBuilder(vc);

        // make new Genotypes based on filters
        if ( !genotypeFilterExps.isEmpty() || setFilteredGenotypesToNocall ) {
            builder.genotypes(makeGenotypes(vc));
        }

        // make a new variant context based on filters
        final Set<String> filters = new LinkedHashSet<>(vc.getFilters());

        // test for clustered SNPs if requested
        if (areClusteredSNPs(featureContext, vc)){
            filters.add(CLUSTERED_SNP_FILTER_NAME);
        }

        for ( final JexlVCMatchExp exp : filterExps ) {
            try {
                if ( invertLogic(VariantContextUtils.match(vc, exp), invertFilterExpression) ) {
                    filters.add(exp.name);
                }
            } catch (final Exception e) {
                // do nothing unless specifically asked to; it just means that the expression isn't defined for this context
                if ( failMissingValues  ) {
                    filters.add(exp.name);
                }
            }
        }

        if ( filters.isEmpty() ) {
            builder.passFilters();
        } else {
            builder.filters(filters);
        }

        writer.add(builder.make());
    }

    private GenotypesContext makeGenotypes(final VariantContext vc) {
        final GenotypesContext genotypes = GenotypesContext.create(vc.getGenotypes().size());

        // for each genotype, check filters then create a new object
        for ( final Genotype g : vc.getGenotypes() ) {
            if ( g.isCalled() ) {
                final List<String> filters = new ArrayList<>();
                if ( g.isFiltered() ) {
                    filters.add(g.getFilters());
                }

                // Add if expression filters the variant context
                for ( final JexlVCMatchExp exp : genotypeFilterExps ) {
                    if ( invertLogic(VariantContextUtils.match(vc, g, exp), invertGenotypeFilterExpression) ) {
                        filters.add(exp.name);
                    }
                }

                // if sample is filtered and --setFilteredGtToNocall, set genotype to non-call
                if ( !filters.isEmpty() && setFilteredGenotypesToNocall ) {
                    genotypes.add(new GenotypeBuilder(g).filters(filters).alleles(diploidNoCallAlleles).make());
                } else {
                    genotypes.add(new GenotypeBuilder(g).filters(filters).make());
                }
            } else {
                genotypes.add(g);
            }
        }
        return genotypes;
    }

    /**
     * Return true if there is a window of size {@link clusterWindow} that contains as least {@link clusterSize} SNPs.
     */
    private boolean areClusteredSNPs(final FeatureContext featureContext, final VariantContext current) {
        if (clusterWindow < 1){ //as per argument doc, snpsInVicinity < 1 imply no clustering
            return false;
        }
        if (!current.isSNP()){
            return false;
        }
        //Need to fetch SNPs from left and right of the current SNPs.
        //Note: The "vicinity" here is the region around the current variant (up and down n bases, where n = clusterWindow).
        final List<VariantContext> snpsInVicinity = featureContext.getValues(getDrivingVariantsFeatureInput(), clusterWindow, clusterWindow)
                                                    .stream().filter(v -> v.isSNP()).collect(Collectors.toList());

        if (snpsInVicinity.size() < clusterSize){  //not enough variants - will never be a cluster no matter what.
            return false;
        }
        snpsInVicinity.sort(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR); //Note: by construction, they will be from the same contig.

        //Scan the list of variants and check if any stretch of clusterSize variants
        //that contains the 'current variant' has more variants than clusterWindow.
        final int currentStartPos = current.getStart();
        int firstIndex = 0;

        //start with the lowest pos
        int firstPos = snpsInVicinity.get(firstIndex).getStart();//there is at least 1 element so it's ok
        final int n = clusterSize - 1;
        while(firstPos <= currentStartPos && firstIndex + n < snpsInVicinity.size()) {
            final int firstPlusNPos = snpsInVicinity.get(firstIndex + n).getStart();

            if (firstPlusNPos - firstPos < clusterWindow){
                return true;
            }

            firstIndex++;
            firstPos = snpsInVicinity.get(firstIndex).getStart();//there is at least 1 element so it's ok
        }
        return false;
    }

    @Override
    public void closeTool() {
        if (writer != null){
            writer.close();
        }
    }

}
