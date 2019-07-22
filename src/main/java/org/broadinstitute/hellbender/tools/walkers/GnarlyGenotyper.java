package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.tools.evoquer.GnarlyGenotyperEngine;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBOptions;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_QualByDepth;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StrandBiasTest;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import org.reflections.Reflections;

import java.io.File;
import java.util.*;

/**
 * Perform "quick and dirty" joint genotyping on one or more samples pre-called with HaplotypeCaller
 *
 * <p>
 * This tool is designed to perform joint genotyping on multiple samples pre-called with HaplotypeCaller to produce a
 * multi-sample callset in a super extra highly scalable manner. In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller
 * with `-ERC GVCF` or `-ERC BP_RESOLUTION`.
 *
 * <h3>Input</h3>
 * <p>
 * A GenomicsDB containing the samples to joint-genotype.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A final VCF in which all samples have been jointly genotyped.
 * </p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Perform joint genotyping on a set of GVCFs stored in a GenomicsDB</h4>
 * <pre>
 * gatk --javaOptions "-Xmx4g" GnarlyGenotyper \
 *   -R reference.fasta \
 *   -V gendb://genomicsdb \
 *   --onlyOutputCallsStartingInIntervals \
 *   -O output.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <p><ul><li>Only GenomicsDB instances can be used as input for this tool.</li>
 * <li>To generate all the annotations necessary for VQSR, input variants must include the QUALapprox, VarDP and MQ_DP
 * </ul></p>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool assumes all diploid genotypes.</p>
 *
 */
@CommandLineProgramProperties(summary = "Perform \"quick and dirty\" joint genotyping on one or more samples pre-called with HaplotypeCaller",
        oneLineSummary = "Perform \"quick and dirty\" joint genotyping on one or more samples pre-called with HaplotypeCaller",
        programGroup = ShortVariantDiscoveryProgramGroup.class)
@DocumentedFeature
public final class GnarlyGenotyper extends VariantWalker {

    private static final OneShotLogger warning = new OneShotLogger(GnarlyGenotyper.class);

    private static final boolean SUMMARIZE_PLs = false;  //for very large numbers of samples, save on space and hail import time by summarizing PLs with genotype quality metrics

    public static final int PIPELINE_MAX_ALT_COUNT = 6;

    private static GnarlyGenotyperEngine genotyperEngine;

    private final RMSMappingQuality mqCalculator = RMSMappingQuality.getInstance();

    private final Set<Class<? extends InfoFieldAnnotation>> allASAnnotations = new HashSet<>();


    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written", optional=false)
    private File outputFile;

    @Argument(fullName = "output-database-name", shortName = "output-db",
            doc="File to which the sites-only annotation database derived from these input samples should be written", optional=true)
    private String outputDbName = null;

    @ArgumentCollection
    private GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    @Argument(fullName = "keep-all-sites", shortName = "keep-all",
            doc="Retain low quality and non-variant sites, applying appropriate filters", optional=true)
    private boolean keepAllSites = false;

    /**
     * This option can only be activated if intervals are specified.
     */
    @Advanced
    @Argument(fullName = GenotypeGVCFs.ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME,
            doc="Restrict variant output to sites that start within provided intervals",
            optional=true)
    private boolean onlyOutputCallsStartingInIntervals = false;

    @Argument(fullName = GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME,
            shortName = GenomicsDBImport.MERGE_INPUT_INTERVALS_LONG_NAME,
            doc = "Boolean flag to read in all data in between intervals.  Improves performance reading from GenomicsDB " +
                    "using large lists of +intervals, as in exome sequencing, especially if GVCF data only exists for " +
                    "specified intervals.")
    private boolean mergeInputIntervals = false;

    @Hidden
    @Argument(fullName = "strip-allele-specific-annotations", shortName = "strip-as", doc = "Remove raw AS values and don't calculate finalized values")
    private boolean stripASAnnotations = false;


    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set
     * when appropriate. Note that dbSNP is not used in any way for the genotyping calculations themselves.
     */
    @ArgumentCollection
    private final DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    private VariantContextWriter vcfWriter;
    private VariantContextWriter annotationDBwriter = null;

    /** these are used when {@link #onlyOutputCallsStartingInIntervals) is true */
    private List<SimpleInterval> intervals;

    @Override
    public boolean requiresReference() {
        return true;
    }

    /**
     * Get the largest interval per contig that contains the intervals specified on the command line.
     * @param getIntervals intervals to be transformed
     * @param sequenceDictionary used to validate intervals
     * @return a list of one interval per contig spanning the input intervals after processing and validation
     */
    @Override
    protected List<SimpleInterval> transformTraversalIntervals(final List<SimpleInterval> getIntervals, final SAMSequenceDictionary sequenceDictionary) {
        if (mergeInputIntervals) {
            return IntervalUtils.getSpanningIntervals(getIntervals, sequenceDictionary);
        } else {
            return getIntervals;
        }
    }

    @Override
    protected GenomicsDBOptions getGenomicsDBOptions() {
        if (genomicsDBOptions == null) {
            genomicsDBOptions = new GenomicsDBOptions(referenceArguments.getReferencePath(), true, PIPELINE_MAX_ALT_COUNT);
        }
        return genomicsDBOptions;
    }

    @Override
    public void onTraversalStart() {
        final VCFHeader inputVCFHeader = getHeaderForVariants();

        if(onlyOutputCallsStartingInIntervals) {
            if( !intervalArgumentCollection.intervalsSpecified()) {
                throw new CommandLineException.MissingArgument("-L or -XL", "Intervals are required if --" + GenotypeGVCFs.ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME + " was specified.");
            }
        }
        intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary()) :
                Collections.emptyList();

        final SampleList samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples());

        setupVCFWriter(inputVCFHeader, samples);

        genotyperEngine = new GnarlyGenotyperEngine(keepAllSites, genotypeArgs.MAX_ALTERNATE_ALLELES, SUMMARIZE_PLs, stripASAnnotations);

        Reflections reflections = new Reflections("org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific");
        //not InfoFieldAnnotation.class because we don't want AS_InbreedingCoeff
        allASAnnotations.addAll(reflections.getSubTypesOf(AS_StrandBiasTest.class));
        allASAnnotations.addAll(reflections.getSubTypesOf(AS_RankSumTest.class));
        allASAnnotations.add(AS_RMSMappingQuality.class);
        allASAnnotations.add(AS_QualByDepth.class);
    }

    private void setupVCFWriter(VCFHeader inputVCFHeader, SampleList samples) {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        // Remove GCVFBlocks
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCFWriter.GVCF_BLOCK));

        //add header for new filter
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));

        // add headers for annotations added by this tool
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AC_ADJUSTED_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.FISHER_STRAND_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.STRAND_ODDS_RATIO_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.SB_TABLE_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.EXCESS_HET_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.QUAL_BY_DEPTH_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.RMS_MAPPING_QUALITY_KEY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.AS_ALT_ALLELE_DEPTH_KEY));
        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        vcfWriter = createVCFWriter(outputFile);
        if (outputDbName != null) {
            annotationDBwriter = createVCFWriter(new File(outputDbName));
        }

        final Set<String> sampleNameSet = samples.asSetOfSamples();
        final VCFHeader dbHeader = new VCFHeader(headerLines);
        if (SUMMARIZE_PLs) {
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.GENOTYPE_QUALITY_BY_ALLELE_BALANCE));
            headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.GENOTYPE_QUALITY_BY_ALT_CONFIDENCE));
        }
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        vcfWriter.writeHeader(vcfHeader);
        if (outputDbName != null) {
            annotationDBwriter.writeHeader(dbHeader);
        }
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    @Override
    public void apply(VariantContext variant, ReadsContext reads, ReferenceContext ref, FeatureContext features) {
        SimpleInterval variantStart = new SimpleInterval(variant.getContig(), variant.getStart(), variant.getStart());
        //return early if there's no non-symbolic ALT since GDB already did the merging
        if ( !variant.isVariant() || !GenotypeGVCFs.isProperlyPolymorphic(variant)
                || variant.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) == 0
                || (onlyOutputCallsStartingInIntervals && !intervals.stream().anyMatch(interval -> interval.contains(variantStart)))) {
            if (keepAllSites) {
                VariantContextBuilder builder = new VariantContextBuilder(mqCalculator.finalizeRawMQ(variant));  //don't fill in QUAL here because there's no alt data
                builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
                builder.attribute(GATKVCFConstants.AC_ADJUSTED_KEY, 0);
                vcfWriter.add(builder.make());
            }
            return;
        }

        //return early if variant doesn't meet QUAL threshold
        if (!variant.hasAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY)) {
            warning.warn("Variant will not be output because it is missing the " + GATKVCFConstants.RAW_QUAL_APPROX_KEY + "key assigned by the ReblockGVCFs tool -- if the input did come from ReblockGVCFs, check the GenomicsDB vidmap.json annotation info");
            return;
        }

        final VariantContextBuilder annotationDBBuilder = new VariantContextBuilder(variant);
        final VariantContext finalizedVC = GnarlyGenotyperEngine.finalizeGenotype(variant, annotationDBBuilder);
        if (finalizedVC != null && (!onlyOutputCallsStartingInIntervals || intervals.stream().anyMatch(interval -> interval.contains(variantStart)))) {
            vcfWriter.add(finalizedVC);
        }
        if (annotationDBwriter != null) {
            annotationDBwriter.add(annotationDBBuilder.make());  //we don't seem to have a sites-only option anymore, so do it manually
        }
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null) {
            vcfWriter.close();
        }
        if (annotationDBwriter != null) {
            annotationDBwriter.close();
        }
    }
}