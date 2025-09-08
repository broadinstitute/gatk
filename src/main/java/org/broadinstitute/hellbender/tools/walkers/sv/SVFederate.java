package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SelectSVPairs;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngine;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.util.*;


/**
 * TODO: docs
 */
@CommandLineProgramProperties(
        summary = "Merges structural variants from two distinct callsets",
        oneLineSummary = "SV federation",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class SVFederate extends MultiVariantWalker {

    @ArgumentCollection
    protected SVFederationVariantsArgumentCollection inputArgs = new SVFederationVariantsArgumentCollection();

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected GATKPath outputFile;

    /**
     * Expected format is tab-delimited and contains columns VID_A, VID_B, SCORE.
     * First line must be a header with column names. Comment lines starting with
     * {@link TableUtils#COMMENT_PREFIX} are ignored.
     */
    public static final String SV_PAIR_FILE_LONG_NAME = "sv-pairs";
    @Argument(
            doc = "SV pair file (.tsv) containing the candidate SV pairs and matching scores",
            fullName = SV_PAIR_FILE_LONG_NAME
    )
    public GATKPath svPairFilePath;

    public static final String PREFIX_A_LONG_NAME = "prefix-A";
    @Argument(
            doc = "Prefix for fields relating to the records from VCF A",
            fullName = PREFIX_A_LONG_NAME
    )
    protected String prefixA = null;

    public static final String PREFIX_B_LONG_NAME = "prefix-B";
    @Argument(
            doc = "Prefix for fields relating to the records from VCF B",
            fullName = PREFIX_B_LONG_NAME
    )
    protected String prefixB = null;

    public static final String VARIANT_PREFIX_LONG_NAME = "variant-prefix";
    @Argument(
            doc = "If supplied, generate variant IDs with this prefix",
            fullName = VARIANT_PREFIX_LONG_NAME,
            optional = true
    )
    protected String variantPrefix = null;

    public static final String MAX_RECORDS_IN_RAM_LONG_NAME = "max-records-in-ram";
    @Argument(fullName = MAX_RECORDS_IN_RAM_LONG_NAME,
            doc = "When writing VCF files that need to be sorted, this will specify the number of records stored in " +
                    "RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a " +
                    "VCF file, and increases the amount of RAM needed.",
            optional=true)
    public int maxRecordsInRam = 10000;

    protected SAMSequenceDictionary dictionary;
    protected ReferenceSequenceFile reference;
    protected SortingCollection<VariantContext> sortingBuffer;
    protected VariantContextWriter writer;
    protected VCFHeader header;
    protected int numVariantsBuilt = 0;

    protected Map<String, Map<String, String>> sourceToPairMap;
    protected HashMap<String, String> vidAtoB;
    protected HashMap<String, String> vidBtoA;

    protected Map<String, List<Map<String, VariantContext>>> sourceToVariantMap;
    protected HashMap<String, VariantContext> vidToRecA;
    protected HashMap<String, VariantContext> vidToRecB;
    protected HashMap<String, String> sourceToPrefixMap;

    protected CanonicalSVCollapser collapser;
    protected CanonicalSVLinkage<SVCallRecord> linkage;


    @Override
    protected MultiVariantInputArgumentCollection getMultiVariantInputArgumentCollection() {
        return new MultiVariantInputArgumentCollection() {
            private static final long serialVersionUID = 1L;

            @Override
            public List<GATKPath> getDrivingVariantPaths() {
                // driving variants will be determined by initializeDrivingVariants();
                // directly overriding getMultiVariantInputArgumentCollection
                // to return an instance of SVFederationVariantsArgumentCollection() did not work
                return Collections.emptyList();
            }
        };
    }


    @Override
    protected void initializeDrivingVariants() {
        getDrivingVariantsFeatureInputs().addAll(inputArgs.getFeatureInputsForDrivingVariants());

        super.initializeDrivingVariants();
    }


    @Override
    public boolean requiresReference() {
        return true;
    }


    @Override
    public void onTraversalStart() {
        super.onTraversalStart();  // loads ploidy table, reference dictionary, initializes writer

        reference = ReferenceUtils.createReferenceReader(referenceArguments.getReferenceSpecifier());
        dictionary = reference.getSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        writer = createVCFWriter(outputFile);
        header = createHeader();
        writer.writeHeader(header);
        sortingBuffer = SortingCollection.newInstance(
                VariantContext.class,
                new VCFRecordCodec(header, true),
                header.getVCFRecordComparator(),
                maxRecordsInRam,
                tmpDir.toPath());

        // get map of A and B feature input names to their respective VID-to-VariantContext maps
        vidToRecA = new HashMap<>();
        vidToRecB = new HashMap<>();
        sourceToVariantMap = new HashMap<>();
        final String sourceA = getDrivingVariantsFeatureInputs().get(0).getName();
        final String sourceB = getDrivingVariantsFeatureInputs().get(1).getName();
        sourceToVariantMap.put(sourceA, Arrays.asList(vidToRecA, vidToRecB));
        sourceToVariantMap.put(sourceB, Arrays.asList(vidToRecB, vidToRecA));

        // get map of A and B feature input names to their respective prefixes
        sourceToPrefixMap = new HashMap<>();
        sourceToPrefixMap.put(sourceA, prefixA);
        sourceToPrefixMap.put(sourceB, prefixB);

        // load SV pairs
        final SelectSVPairs selector = new SelectSVPairs(svPairFilePath);
        vidAtoB = selector.getVidAToBMap();
        vidBtoA = selector.getVidBToAMap();

        // get map of A and B feature input names to their respective SV match maps
        // order is fixed in SVFederationArgumentCollection.getDrivingVariantPaths()
        sourceToPairMap = new HashMap<>();
        sourceToPairMap.put(getDrivingVariantsFeatureInputs().get(0).getName(), vidAtoB);
        sourceToPairMap.put(getDrivingVariantsFeatureInputs().get(1).getName(), vidBtoA);

        reference = ReferenceUtils.createReferenceReader(referenceArguments.getReferenceSpecifier());
        // TODO breakpoint summary strategy
        collapser = new CanonicalSVCollapser(reference,
                CanonicalSVCollapser.AltAlleleSummaryStrategy.MOST_SPECIFIC_SUBTYPE,
                CanonicalSVCollapser.BreakpointSummaryStrategy.REPRESENTATIVE,
                CanonicalSVCollapser.FlagFieldLogic.OR);
        linkage = new CanonicalSVLinkage<>(dictionary, true);
    }


    protected VCFHeader createHeader() {
        final VCFHeader header = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder());
        header.setSequenceDictionary(dictionary);

        header.addMetaDataLine(new VCFInfoHeaderLine(prefixA + "_VID", 1, VCFHeaderLineType.String,
                "Variant ID in " + prefixA));
        header.addMetaDataLine(new VCFInfoHeaderLine(prefixB + "_VID", 1, VCFHeaderLineType.String,
                "Variant ID in " + prefixB));
        header.addMetaDataLine(new VCFInfoHeaderLine(prefixA + "_AF", VCFHeaderLineCount.A, VCFHeaderLineType.Float,
                "Allele frequency in " + prefixA));
        header.addMetaDataLine(new VCFInfoHeaderLine(prefixB + "_AF", VCFHeaderLineCount.A, VCFHeaderLineType.Float,
                "Allele frequency in " + prefixB));
        header.addMetaDataLine(new VCFInfoHeaderLine(prefixA + "_AN", 1, VCFHeaderLineType.Integer,
                "Allele number in " + prefixA));
        header.addMetaDataLine(new VCFInfoHeaderLine(prefixB + "_AN", 1, VCFHeaderLineType.Integer,
                "Allele number in " + prefixB));
        header.addMetaDataLine(new VCFInfoHeaderLine(prefixA + "_AC", VCFHeaderLineCount.A, VCFHeaderLineType.Integer,
                "Allele count in " + prefixA));
        header.addMetaDataLine(new VCFInfoHeaderLine(prefixB + "_AC", VCFHeaderLineCount.A, VCFHeaderLineType.Integer,
                "Allele count in " + prefixB));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.RECIPROCAL_OVERLAP_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Reciprocal overlap between merged variants"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SIZE_SIMILARITY_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Size similarity between merged variants"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BREAKPOINT_DISTANCE_START_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Distance in bp between start coordinates of merged variants"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BREAKPOINT_DISTANCE_END_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Distance in bp between end coordinates of merged variants"));
        header.addMetaDataLine(new VCFInfoHeaderLine("LOG_AF_DIFFERENCE", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Absolute value of the difference of log allele frequencies between variant sources"));

        return header;
    }

    public static Double computeLogAlleleFrequencyDifference(final double alleleFrequencyA,
                                                             final double alleleFrequencyB) {
        if (Double.isNaN(alleleFrequencyA) || Double.isNaN(alleleFrequencyB)) {
            return Double.NaN;
        } else {
            return Math.abs(Math.log10(alleleFrequencyA) - Math.log10(alleleFrequencyB));
        }
    }

    protected SVCallRecord merge(final VariantContext thisVariant,
                                 final VariantContext thatVariant) {
        final SVCallRecord thisRecord = SVCallRecordUtils.create(thisVariant, true, false, dictionary);
        final SVCallRecord thatRecord = SVCallRecordUtils.create(thatVariant, true, false, dictionary);

        final String thisPrefix = sourceToPrefixMap.get(thisVariant.getSource());
        final String thatPrefix = sourceToPrefixMap.get(thatVariant.getSource());

        final SVClusterEngine.OutputCluster outputCluster =
                new SVClusterEngine.OutputCluster(List.of(thisRecord, thatRecord));
        final SVCallRecord merged = collapser.collapse(outputCluster);
        final Map<String, Object> attributes = merged.getAttributes();

        final CanonicalSVLinkage.CanonicalLinkageResult result = linkage.areClusterable(thisRecord, thatRecord);
        attributes.put(GATKSVVCFConstants.RECIPROCAL_OVERLAP_INFO, result.getReciprocalOverlap());
        attributes.put(GATKSVVCFConstants.SIZE_SIMILARITY_INFO, result.getSizeSimilarity());
        attributes.put(GATKSVVCFConstants.BREAKPOINT_DISTANCE_START_INFO, result.getBreakpointDistance1());
        attributes.put(GATKSVVCFConstants.BREAKPOINT_DISTANCE_END_INFO, result.getBreakpointDistance2());

        // store cohort AFs and calculate total AF
        // TODO: handle mCNVs
        final int thisAN = thisVariant.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY, 0);
        final int thatAN = thatVariant.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY, 0);
        final int totalAN = thisAN + thatAN;
        attributes.put(thisPrefix + "_AN", thisAN);
        attributes.put(thatPrefix + "_AN", thatAN);
        attributes.put(VCFConstants.ALLELE_NUMBER_KEY, totalAN);

        final int thisAC = thisVariant.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0);
        final int thatAC = thatVariant.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0);
        final int totalAC = thisAC + thatAC;
        attributes.put(thisPrefix + "_AC", thisAC);
        attributes.put(thatPrefix + "_AC", thatAC);
        attributes.put(VCFConstants.ALLELE_COUNT_KEY, totalAC);

        final double thisAF = thisVariant.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, Double.NaN);
        final double thatAF = thatVariant.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, Double.NaN);
        final double totalAF = (double) totalAC / totalAN;
        attributes.put(thisPrefix + "_AF", thisAF);
        attributes.put(thatPrefix + "_AF", thatAF);
        attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, totalAF);
        attributes.put("LOG_AF_DIFFERENCE", computeLogAlleleFrequencyDifference(thisAF, thatAF));  // TODO keep? if so move key to constants

        // store original VIDs
        attributes.put(thisPrefix + "_VID", thisRecord.getId());
        attributes.put(thatPrefix + "_VID", thatRecord.getId());

        return merged;
    }

    @Override
    public Object onTraversalSuccess() {
        for (final VariantContext variant : sortingBuffer) {
            writer.add(variant);
        }
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (sortingBuffer != null) {
            sortingBuffer.cleanup();
        }
        if (writer != null) {
            writer.close();
        }
    }

    protected void write(final SVCallRecord call) {
        sortingBuffer.add(buildVariantContext(call));
    }

    protected VariantContext buildVariantContext(final SVCallRecord call) {
        // Assign new variant ID
        final String newId = variantPrefix == null ? call.getId() : String.format("%s%08x", variantPrefix, numVariantsBuilt++);

        // Build new variant with empty genotypes
        final SVCallRecord finalCall = new SVCallRecord(newId, call.getContigA(), call.getPositionA(), call.getStrandA(),
                call.getContigB(), call.getPositionB(), call.getStrandB(), call.getType(), call.getComplexSubtype(),
                call.getComplexEventIntervals(), call.getLength(), call.getEvidence(), call.getAlgorithms(), call.getAlleles(), Collections.emptyList(),
                call.getAttributes(), call.getFilters(), call.getLog10PError(), dictionary);
        final VariantContextBuilder builder = SVCallRecordUtils.getVariantBuilder(finalCall);
        builder.rmAttribute(GATKSVVCFConstants.CLUSTER_MEMBER_IDS_KEY);
        return builder.make();
    }


    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext,
                      FeatureContext featureContext) {
        final String source = variant.getSource();  // which callset is this variant from
        final Map<String, String> pairMap = sourceToPairMap.get(source);  // A-->B or B-->A matching map
        // map VID to variant for this callset (to store this variant) and the other callset (to retrieve its match)
        final Map<String, VariantContext> thisVariantMap = sourceToVariantMap.get(source).get(0);
        final Map<String, VariantContext> thatVariantMap = sourceToVariantMap.get(source).get(1);
        final String vid = variant.getID();
        if (pairMap.containsKey(vid)) {
            // get VID of matching variant in other callset
            final String match = pairMap.get(vid);
            if (thatVariantMap.containsKey(match)) {
                // retrieve the matching variant which was saved in thatVariantMap, merge them, and write
                // TODO: external AF annotation only mode
                write(merge(variant, thatVariantMap.get(match)));  // write() handles conversion to VariantContext and sorting in buffer
                thatVariantMap.remove(match);  // all done with this variant - delete it to save memory
            } else {
                thisVariantMap.put(vid, variant);  // store this variant until its match is read
            }
        } else {
            // if the variant has no match, output to sorting buffer
            write(SVCallRecordUtils.create(variant, true, false, dictionary));
        }
    }
}
