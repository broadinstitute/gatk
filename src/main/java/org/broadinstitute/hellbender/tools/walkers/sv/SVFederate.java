package org.broadinstitute.hellbender.tools.walkers.sv;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.StructuralVariantAnnotationType;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SelectSVPairs;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngine;
import org.broadinstitute.hellbender.tools.sv.cluster.SVFederationCollapser;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFRecordCodec;


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

    public static final String AF_GROUPINGS_A_LONG_NAME = "af-groupings-A";
    @Argument(fullName = AF_GROUPINGS_A_LONG_NAME,
            doc = "Comma-separated list of groups to stratify AFs by that are annotated in VCF A.",
            optional=true)
    public String afGroupingsAInput = null;

    public static final String AF_GROUPINGS_B_LONG_NAME = "af-groupings-B";
    @Argument(fullName = AF_GROUPINGS_B_LONG_NAME,
            doc = "Comma-separated list of groups to stratify AFs by that are annotated in VCF B.",
            optional=true)
    public String afGroupingsBInput = null;

    public static final String XY_LONG_NAME = "XY-identifier";
    @Argument(fullName = XY_LONG_NAME,
            doc = "String used to annotate frequency information among XY individuals in input VCFs.",
            optional=true)
    public String xyIdentifier = null;

    public static final String XX_LONG_NAME = "XX-identifier";
    @Argument(fullName = XX_LONG_NAME,
            doc = "String used to annotate frequency information among XX individuals in input VCFs.",
            optional=true)
    public String xxIdentifier = null;

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
    protected Map<String, List<String>> sourceToAFGroupingsMap;

    protected List<String> sexes;
    protected List<String> afGroupingsA;
    protected List<String> afGroupingsB;
    protected List<String> afGroupingsAll;

    protected CanonicalSVCollapser collapser;
    protected CanonicalSVLinkage<SVCallRecord> linkage;


    private static final class VCFHeaderLineBuilder {
        private final String key;
        private final VCFHeaderLineCount count;
        private final int countInt;
        private final VCFHeaderLineType type;
        private final String description;

        private VCFHeaderLineBuilder(String key, VCFHeaderLineCount count, VCFHeaderLineType type, String description) {
            this.key = key;
            this.count = count;
            this.countInt = -1;
            this.type = type;
            this.description = description;
        }

        private VCFHeaderLineBuilder(String key, int count, VCFHeaderLineType type, String description) {
            this.key = key;
            this.countInt = count;
            this.count = null;
            this.type = type;
            this.description = description;
        }

        private String getKey() { return key; }
        private String getKeyWithPrefix(String prefix) { return prefix + "_" + key; }
        private String getKeyWithPrefixes(String prefix, String group1, String group2) {
            if (group2 == null || group2.isEmpty()) {
                return getKeyWithPrefix(prefix) + "_" + group1;
            }
            return getKeyWithPrefix(prefix) + "_" + group1 + "_" + group2;
        }

        private VCFInfoHeaderLine addPrefix(String prefix) {
            final String keyWithPrefix = getKeyWithPrefix(prefix);
            final String descriptionWithPrefix = description + " in " + prefix + ".";
            if (count != null) {
                return new VCFInfoHeaderLine(keyWithPrefix, count, type, descriptionWithPrefix);
            } else {
                return new VCFInfoHeaderLine(keyWithPrefix, countInt, type, descriptionWithPrefix);
            }
        }

        private VCFInfoHeaderLine addPrefixes(String prefix, String group1, String group2) {
            final String keyWithPrefix = getKeyWithPrefixes(prefix, group1, group2);
            final String descriptionWithPrefix = description + " in " + group1 + ((group2 == null || group2.isEmpty()) ? "" : " " + group2) + " individuals in " + prefix + ".";
            if (count != null) {
                return new VCFInfoHeaderLine(keyWithPrefix, count, type, descriptionWithPrefix);
            } else {
                return new VCFInfoHeaderLine(keyWithPrefix, countInt, type, descriptionWithPrefix);
            }
        }
    }
    /**
     * INFO fields to annotate from each source cohort
     */
    private static final List<VCFHeaderLineBuilder> COHORT_INFO_FIELDS = List.of(
            new VCFHeaderLineBuilder("VID", 1, VCFHeaderLineType.String, "Variant ID"),
            new VCFHeaderLineBuilder(GATKSVVCFConstants.SVTYPE, 1, VCFHeaderLineType.String, "SV type"),
            new VCFHeaderLineBuilder(GATKSVVCFConstants.CPX_TYPE, 1, VCFHeaderLineType.String, "Complex SV subtype")
    );

    private static final List<VCFHeaderLineBuilder> BIALLELIC_COHORT_INFO_FIELDS = List.of(
        new VCFHeaderLineBuilder(VCFConstants.ALLELE_COUNT_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Allele count"),
        new VCFHeaderLineBuilder(VCFConstants.ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele frequency"),
        new VCFHeaderLineBuilder(VCFConstants.ALLELE_NUMBER_KEY, 1, VCFHeaderLineType.Integer, "Allele number"),
        new VCFHeaderLineBuilder("N_HOMREF", 1, VCFHeaderLineType.Integer, "Number of samples with homozygous reference genotypes (biallelic sites only)"),
        new VCFHeaderLineBuilder("N_HET", 1, VCFHeaderLineType.Integer, "Number of samples with heterozygous genotypes (biallelic sites only)"),
        new VCFHeaderLineBuilder("N_HOMALT", 1, VCFHeaderLineType.Integer, "Number of samples with homozygous alternate genotypes (biallelic sites only)")
    );

    private static final List<VCFHeaderLineBuilder> BIALLELIC_XY_COHORT_INFO_FIELDS = List.of(
        new VCFHeaderLineBuilder("N_HEMIALT", 1, VCFHeaderLineType.Integer, "Number of XY samples with hemizygous alternate genotypes (biallelic sites only)")
    );

    private static final List<VCFHeaderLineBuilder> MULTIALLELIC_COHORT_INFO_FIELDS = List.of(
        new VCFHeaderLineBuilder("RD_CN_ESTIMATED_AF", 1, VCFHeaderLineType.Float, "Estimated AF from RD_CN"),
        new VCFHeaderLineBuilder("CN_NUMBER", 1, VCFHeaderLineType.Integer, "Total number of samples with estimated copy numbers (multiallelic CNVs only)"),
        new VCFHeaderLineBuilder("CN_COUNT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Number of samples observed at each copy state, starting from CN=0 (multiallelic CNVs only)"),
        new VCFHeaderLineBuilder("CN_NONREF_COUNT", 1, VCFHeaderLineType.Integer, "Number of samples with non-reference copy states (multiallelic CNVs only)"),
        new VCFHeaderLineBuilder("CN_NONREF_FREQ", 1, VCFHeaderLineType.Float, "Frequency of samples with non-reference copy states (multiallelic CNVs only)"),
        new VCFHeaderLineBuilder("CN_FREQ", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Frequency of samples observed at each copy state, starting from CN=0 (multiallelic CNVs only)")
    );

    private static final List<VCFHeaderLineBuilder> ALL_COHORT_INFO_FIELDS = Stream.of(
        COHORT_INFO_FIELDS, BIALLELIC_COHORT_INFO_FIELDS, MULTIALLELIC_COHORT_INFO_FIELDS)
        .flatMap(List::stream)
        .collect(Collectors.toList());

    private static final List<VCFHeaderLineBuilder> INFO_FIELDS_TO_STRATIFY = Stream.of(
        BIALLELIC_COHORT_INFO_FIELDS, MULTIALLELIC_COHORT_INFO_FIELDS)
        .flatMap(List::stream)
        .collect(Collectors.toList());

    private static final List<VCFHeaderLineBuilder> ALL_BIALLELIC_INFOS = Stream.of(
        BIALLELIC_COHORT_INFO_FIELDS, BIALLELIC_XY_COHORT_INFO_FIELDS)
        .flatMap(List::stream)
        .collect(Collectors.toList());

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

        // Parse comma-separated AF groupings inputs into lists
        if (afGroupingsAInput != null && !afGroupingsAInput.isEmpty()) {
            afGroupingsA = Arrays.stream(afGroupingsAInput.split(","))
                    .map(String::trim)
                    .filter(s -> !s.isEmpty())
                    .collect(Collectors.toList());
        }
        if (afGroupingsBInput != null && !afGroupingsBInput.isEmpty()) {
            afGroupingsB = Arrays.stream(afGroupingsBInput.split(","))
                    .map(String::trim)
                    .filter(s -> !s.isEmpty())
                    .collect(Collectors.toList());
        }

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

        // get map of A and B feature input names to their respective AF groupings
        sourceToAFGroupingsMap = new HashMap<>();
        sourceToAFGroupingsMap.put(sourceA, afGroupingsA);
        sourceToAFGroupingsMap.put(sourceB, afGroupingsB);

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
        collapser = new SVFederationCollapser(reference,
                CanonicalSVCollapser.AltAlleleSummaryStrategy.MOST_SPECIFIC_SUBTYPE,
                CanonicalSVCollapser.BreakpointSummaryStrategy.REPRESENTATIVE,
                CanonicalSVCollapser.FlagFieldLogic.OR);
        linkage = new CanonicalSVLinkage<>(dictionary, true);
    }


    protected VCFHeader createHeader() {
        final VCFHeader header = new VCFHeader(getHeaderForVariants().getMetaDataInInputOrder());
        header.setSequenceDictionary(dictionary);

        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.RECIPROCAL_OVERLAP_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Reciprocal overlap between merged variants"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.SIZE_SIMILARITY_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "Size similarity between merged variants"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BREAKPOINT_DISTANCE_START_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Distance in bp between start coordinates of merged variants"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.BREAKPOINT_DISTANCE_END_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Distance in bp between end coordinates of merged variants"));
        header.addMetaDataLine(new VCFInfoHeaderLine("LOG_AF_DIFFERENCE", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Absolute value of the difference of log allele frequencies between variant sources"));

        // Add prefixed cohort-specific INFO fields
        for (final VCFHeaderLineBuilder line : ALL_COHORT_INFO_FIELDS) {
            header.addMetaDataLine(line.addPrefix(prefixA));
            header.addMetaDataLine(line.addPrefix(prefixB));
        }

        sexes = new ArrayList<>();
        if (xyIdentifier != null) {
            sexes.add(xyIdentifier);
        }
        if (xxIdentifier != null) {
            sexes.add(xxIdentifier);
        }

        afGroupingsAll = Stream.concat(afGroupingsA.stream(), afGroupingsB.stream())
            .distinct()
            .collect(Collectors.toList());

        // add sex and group stratified frequency info fields per cohort
        for (final VCFHeaderLineBuilder line : INFO_FIELDS_TO_STRATIFY) {
            for (final String sex : sexes) {
                header.addMetaDataLine(line.addPrefixes(prefixA, sex, ""));
                header.addMetaDataLine(line.addPrefixes(prefixB, sex, ""));
            }
            for (final String grouping : afGroupingsA) {
                header.addMetaDataLine(line.addPrefixes(prefixA, grouping, ""));
                for (final String sex : sexes) {
                    header.addMetaDataLine(line.addPrefixes(prefixA, grouping, sex));
                }
            }
            for (final String grouping : afGroupingsB) {
                header.addMetaDataLine(line.addPrefixes(prefixB, grouping, ""));
                for (final String sex : sexes) {
                    header.addMetaDataLine(line.addPrefixes(prefixB, grouping, sex));
                }
            }
        }

        // add xy-specific frequency info fields per cohort stratified by groups
        for (final VCFHeaderLineBuilder line : BIALLELIC_XY_COHORT_INFO_FIELDS) {
            header.addMetaDataLine(line.addPrefixes(prefixA, xyIdentifier, ""));
            header.addMetaDataLine(line.addPrefixes(prefixB, xyIdentifier, ""));
            for (final String grouping : afGroupingsA) {
                header.addMetaDataLine(line.addPrefixes(prefixA, grouping, xyIdentifier));
            }
            for (final String grouping : afGroupingsB) {
                header.addMetaDataLine(line.addPrefixes(prefixB, grouping, xyIdentifier));
            }
        }

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

    protected void annotateCohortFrequencyInformation(final Map<String, Object> attributes,
                                                      final VariantContext variant,
                                                      final String prefix,
                                                      final List<String> afGroupings,
                                                      final List<VCFHeaderLineBuilder> frequencyInfoLines) {
        // if field (with optional groupings) is present in the cohort variant, annotate with cohort prefix
        // rely on check for presence to only annotate XY fields with XY identifier
        for (final VCFHeaderLineBuilder line : frequencyInfoLines) {
            if (variant.hasAttribute(line.getKey())) {
                attributes.put(line.getKeyWithPrefix(prefix), variant.getAttribute(line.getKey()));
            }
            for (final String sex : sexes) {
                final String key = line.getKey() + "_" + sex;
                if (variant.hasAttribute(key)) {
                    attributes.put(line.getKeyWithPrefixes(prefix, sex, ""), variant.getAttribute(key));
                }
            }
            for (final String grouping : afGroupings) {
                final String key = line.getKey() + "_" + grouping;
                if (variant.hasAttribute(key)) {
                    attributes.put(line.getKeyWithPrefixes(prefix, grouping, ""), variant.getAttribute(key));
                }
                for (final String sex : sexes) {
                    final String keyWithSex = key + "_" + sex;
                    if (variant.hasAttribute(keyWithSex)) {
                        attributes.put(line.getKeyWithPrefixes(prefix, grouping, sex), variant.getAttribute(keyWithSex));
                    }
                }
            }
        }
    }

    protected String getInfoKeyWithGroups(final String baseKey, final String group1, final String group2) {
        return baseKey + ((group1 == null || group1.isEmpty()) ? "" : "_" + group1) + ((group2 == null || group2.isEmpty()) ? "" : "_" + group2);
    }

    protected void biallelicFrequencyMergeHelper(final Map<String, Object> attributes,
                                                 final VariantContext thisVariant,
                                                 final VariantContext thatVariant,
                                                 final String group1, final String group2) {
        final String keyAN = getInfoKeyWithGroups(VCFConstants.ALLELE_NUMBER_KEY, group1, group2);
        final String keyAC = getInfoKeyWithGroups(VCFConstants.ALLELE_COUNT_KEY, group1, group2);
        final String keyAF = getInfoKeyWithGroups(VCFConstants.ALLELE_FREQUENCY_KEY, group1, group2);

        final int thisAN = thisVariant.getAttributeAsInt(keyAN, 0);
        final int thatAN = thatVariant.getAttributeAsInt(keyAN, 0);
        final int totalAN = thisAN + thatAN;
        attributes.put(keyAN, totalAN);

        final int thisAC = thisVariant.getAttributeAsInt(keyAC, 0);
        final int thatAC = thatVariant.getAttributeAsInt(keyAC, 0);
        final int totalAC = thisAC + thatAC;
        attributes.put(keyAC, totalAC);

        if (totalAN > 0) {
            final double totalAF = (double) totalAC / totalAN;
            attributes.put(keyAF, totalAF);
        }

        final String keyNHet = getInfoKeyWithGroups("N_HET", group1, group2);
        final String keyNHomRef = getInfoKeyWithGroups("N_HOMREF", group1, group2);
        final String keyNHomAlt = getInfoKeyWithGroups("N_HOMALT", group1, group2);
        final String keyNHemiAlt = getInfoKeyWithGroups("N_HEMIALT", group1, group2);

        final int thisNHet = thisVariant.getAttributeAsInt(keyNHet, 0);
        final int thatNHet = thatVariant.getAttributeAsInt(keyNHet, 0);
        attributes.put(keyNHet, thisNHet + thatNHet);

        final int thisNHomRef = thisVariant.getAttributeAsInt(keyNHomRef, 0);
        final int thatNHomRef = thatVariant.getAttributeAsInt(keyNHomRef, 0);
        attributes.put(keyNHomRef, thisNHomRef + thatNHomRef);

        final int thisNHomAlt = thisVariant.getAttributeAsInt(keyNHomAlt, 0);
        final int thatNHomAlt = thatVariant.getAttributeAsInt(keyNHomAlt, 0);
        attributes.put(keyNHomAlt, thisNHomAlt + thatNHomAlt);

        // annotate hemialt only if at least one of the variants has the key to safeguard nonsense combinations or autosomes
        if (thisVariant.hasAttribute(keyNHemiAlt) || thatVariant.hasAttribute(keyNHemiAlt)) {
            final int thisNHemiAlt = thisVariant.getAttributeAsInt(keyNHemiAlt, 0);
            final int thatNHemiAlt = thatVariant.getAttributeAsInt(keyNHemiAlt, 0);
            attributes.put(keyNHemiAlt, thisNHemiAlt + thatNHemiAlt);
        }
    }

    protected void annotateFederatedBiallelicFrequencyInformation(final Map<String, Object> attributes,
                                                                  final VariantContext thisVariant,
                                                                  final VariantContext thatVariant) {
        biallelicFrequencyMergeHelper(attributes, thisVariant, thatVariant, "", "");
        for (final String sex: sexes) {
            biallelicFrequencyMergeHelper(attributes, thisVariant, thatVariant, sex, "");
        }
        for (final String group: afGroupingsAll) {
            biallelicFrequencyMergeHelper(attributes, thisVariant, thatVariant, group, "");
            for (final String sex: sexes) {
                biallelicFrequencyMergeHelper(attributes, thisVariant, thatVariant, group, sex);
            }
        }
    }

    protected void multiallelicFrequencyMergeHelper(final Map<String, Object> attributes,
                                                    final VariantContext thisVariant,
                                                    final VariantContext thatVariant,
                                                    final String group1, final String group2) {
        final String keyCnNumber = getInfoKeyWithGroups("CN_NUMBER", group1, group2);
        final String keyCnCount = getInfoKeyWithGroups("CN_COUNT", group1, group2);
        final String keyCnFreq = getInfoKeyWithGroups("CN_FREQ", group1, group2);
        final String keyCnNonrefCount = getInfoKeyWithGroups("CN_NONREF_COUNT", group1, group2);
        final String keyCnNonrefFreq = getInfoKeyWithGroups("CN_NONREF_FREQ", group1, group2);

        final List<Integer> thisCnCounts = thisVariant.getAttributeAsIntList(keyCnCount, 0);
        final List<Integer> thatCnCounts = thatVariant.getAttributeAsIntList(keyCnCount, 0);
        final int maxCnCountSize = Math.max(thisCnCounts == null ? 0 : thisCnCounts.size(), thatCnCounts == null ? 0 : thatCnCounts.size());
        final int[] totalCnCounts = new int[maxCnCountSize];
        for (int i = 0; i < maxCnCountSize; i++) {
            final int thisCount = (thisCnCounts != null && i < thisCnCounts.size()) ? thisCnCounts.get(i) : 0;
            final int thatCount = (thatCnCounts != null && i < thatCnCounts.size()) ? thatCnCounts.get(i) : 0;
            totalCnCounts[i] = thisCount + thatCount;
        }
        attributes.put(keyCnCount, totalCnCounts);

        final int thisCnNumber = thisVariant.getAttributeAsInt(keyCnNumber, 0);
        final int thatCnNumber = thatVariant.getAttributeAsInt(keyCnNumber, 0);
        final int totalCnNumber = thisCnNumber + thatCnNumber;
        attributes.put(keyCnNumber, totalCnNumber);

        final int thisCnNonrefCount = thisVariant.getAttributeAsInt(keyCnNonrefCount, 0);
        final int thatCnNonrefCount = thatVariant.getAttributeAsInt(keyCnNonrefCount, 0);
        final int totalCnNonrefCount = thisCnNonrefCount + thatCnNonrefCount;
        attributes.put(keyCnNonrefCount, totalCnNonrefCount);

        if (totalCnNumber > 0) {
            final double[] totalCnFreqs = new double[maxCnCountSize];
            for (int i = 0; i < maxCnCountSize; i++) {
                totalCnFreqs[i] = (double) totalCnCounts[i] / totalCnNumber;
            }
            attributes.put(keyCnFreq, totalCnFreqs);

            final double totalCnNonrefFreq = (double) totalCnNonrefCount / totalCnNumber;
            attributes.put(keyCnNonrefFreq, totalCnNonrefFreq);
        }
    }


    protected void annotateFederatedMultiallelicFrequencyInformation(final Map<String, Object> attributes,
                                                                     final VariantContext thisVariant,
                                                                     final VariantContext thatVariant) {
        multiallelicFrequencyMergeHelper(attributes, thisVariant, thatVariant, "", "");
        for (final String sex: sexes) {
            multiallelicFrequencyMergeHelper(attributes, thisVariant, thatVariant, sex, "");
        }
        for (final String group: afGroupingsAll) {
            multiallelicFrequencyMergeHelper(attributes, thisVariant, thatVariant, group, "");
            for (final String sex: sexes) {
                multiallelicFrequencyMergeHelper(attributes, thisVariant, thatVariant, group, sex);
            }
        }

    }

    protected SVCallRecord merge(final VariantContext thisVariant,
                                 final VariantContext thatVariant) {
        final SVCallRecord thisRecord = SVCallRecordUtils.create(thisVariant, true, false, dictionary);
        final SVCallRecord thatRecord = SVCallRecordUtils.create(thatVariant, true, false, dictionary);

        final String thisPrefix = sourceToPrefixMap.get(thisVariant.getSource());
        final String thatPrefix = sourceToPrefixMap.get(thatVariant.getSource());

        final List<String> thisAFGroupings = sourceToAFGroupingsMap.get(thisVariant.getSource());
        final List<String> thatAFGroupings = sourceToAFGroupingsMap.get(thatVariant.getSource());

        final SVClusterEngine.OutputCluster outputCluster =
                new SVClusterEngine.OutputCluster(List.of(thisRecord, thatRecord));
        final SVCallRecord merged = collapser.collapse(outputCluster);
        final Map<String, Object> attributes = merged.getAttributes();

        attributes.put(thisPrefix + "_VID", thisRecord.getId());
        attributes.put(thatPrefix + "_VID", thatRecord.getId());

        final CanonicalSVLinkage.CanonicalLinkageResult result = linkage.areClusterable(thisRecord, thatRecord);
        attributes.put(GATKSVVCFConstants.RECIPROCAL_OVERLAP_INFO, result.getReciprocalOverlap());
        attributes.put(GATKSVVCFConstants.SIZE_SIMILARITY_INFO, result.getSizeSimilarity());
        attributes.put(GATKSVVCFConstants.BREAKPOINT_DISTANCE_START_INFO, result.getBreakpointDistance1());
        attributes.put(GATKSVVCFConstants.BREAKPOINT_DISTANCE_END_INFO, result.getBreakpointDistance2());

        // annotate per-cohort INFO fields
        for (final VCFHeaderLineBuilder line : COHORT_INFO_FIELDS) {
            attributes.put(line.getKeyWithPrefix(thisPrefix), thisVariant.hasAttribute(line.getKey()) ? thisVariant.getAttribute(line.getKey()) : VCFConstants.MISSING_VALUE_v4);
            attributes.put(line.getKeyWithPrefix(thatPrefix), thatVariant.hasAttribute(line.getKey()) ? thatVariant.getAttribute(line.getKey()) : VCFConstants.MISSING_VALUE_v4);
        }

        final boolean thisIsCnv = thisRecord.getType() == StructuralVariantAnnotationType.CNV;
        final boolean thatIsCnv = thatRecord.getType() == StructuralVariantAnnotationType.CNV;

        // if biallelic, annotate per-cohort biallelic info fields such as AF if they exist
        if (!thisIsCnv) {
            annotateCohortFrequencyInformation(attributes, thisVariant, thisPrefix, thisAFGroupings, ALL_BIALLELIC_INFOS);
        }

        if (!thatIsCnv) {
            annotateCohortFrequencyInformation(attributes, thatVariant, thatPrefix, thatAFGroupings, ALL_BIALLELIC_INFOS);
        }

        if (thisIsCnv || thatIsCnv) {
            // if either variant is multiallelic, set federated biallelic fields such as AF to missing
            for (final VCFHeaderLineBuilder line : BIALLELIC_COHORT_INFO_FIELDS) {
                attributes.put(line.getKey(), VCFConstants.MISSING_VALUE_v4);
            }
            attributes.put("LOG_AF_DIFFERENCE", VCFConstants.MISSING_VALUE_v4);

            // if either variant is multiallleic, annotate per-cohort multiallelic info fields such as CN_FREQ if they exist
            annotateCohortFrequencyInformation(attributes, thisVariant, thisPrefix, thisAFGroupings, MULTIALLELIC_COHORT_INFO_FIELDS);
            annotateCohortFrequencyInformation(attributes, thatVariant, thatPrefix, thatAFGroupings, MULTIALLELIC_COHORT_INFO_FIELDS);

            // compute federated CN statistics
            annotateFederatedMultiallelicFrequencyInformation(attributes, thisVariant, thatVariant);

        } else {
            // for biallelic variants, compute federated AF and related statistics
            annotateFederatedBiallelicFrequencyInformation(attributes, thisVariant, thatVariant);

            final double thisAF = thisVariant.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, Double.NaN);
            final double thatAF = thatVariant.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, Double.NaN);
            attributes.put("LOG_AF_DIFFERENCE", computeLogAlleleFrequencyDifference(thisAF, thatAF));
        }

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
