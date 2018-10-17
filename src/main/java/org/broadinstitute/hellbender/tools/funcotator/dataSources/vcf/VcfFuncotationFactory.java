package org.broadinstitute.hellbender.tools.funcotator.dataSources.vcf;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.nio.file.Path;
import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

/**
 * A class to create annotations from VCF feature sources.
 * Created by jonn on 3/23/18.
 */
public class VcfFuncotationFactory extends DataSourceFuncotationFactory {

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    /** Standard Logger.  */
    protected static final Logger logger = LogManager.getLogger(VcfFuncotationFactory.class);

    //==================================================================================================================
    // Private Members:

    /**
     * Name of this {@link VcfFuncotationFactory}.
     */
    private final String name;

    /**
     * The file on which this {@link VcfFuncotationFactory} is based.
     */
    private final Path sourceFilePath;

    /**
     * The field names that this {@link VcfFuncotationFactory} supports
     * and default values for each.
     */
    private final LinkedHashMap<String, Object> supportedFieldNamesAndDefaults;

    /**
     * A list of values to use when there are no annotations for an allele.
     */
    private final LinkedHashSet<String> supportedFieldNames;

    /**
     * Should contain metadata only for the fields in supportedFieldNames.
     */
    private final FuncotationMetadata supportedFieldMetadata;

    /**
     * Cache for speed.  Please note that the cache is done on the reference.
     */
    private final LRUCache<Triple<VariantContext, ReferenceContext, List<Feature>>, List<Funcotation>> cache = new LRUCache<>();

    /**
     * If the VCF has multiple lines with the same position, ref, and alt.
     */
    private final static String DUPLICATE_RECORD_DELIMITER = "|";

    @VisibleForTesting
    int cacheHits = 0;
    @VisibleForTesting
    int cacheMisses = 0;

    //==================================================================================================================
    // Constructors:

    /**
     * Create a {@link VcfFuncotationFactory}.
     * @param name A {@link String} containing the name of this {@link VcfFuncotationFactory}.
     * @param version  The version {@link String} of the backing data source from which {@link Funcotation}s will be made.
     * @param sourceFilePath {@link Path} to the VCF file from which {@link VariantContext}s will be read in and used as Features from which to create {@link Funcotation}s.
     * @param annotationOverridesMap A {@link LinkedHashMap<String,String>} containing user-specified overrides for specific {@link Funcotation}s.
     * @param mainSourceFileAsFeatureInput The backing {@link FeatureInput} for this {@link VcfFuncotationFactory}, from which all {@link Funcotation}s will be created.
     */
    public VcfFuncotationFactory(final String name,
                                 final String version,
                                 final Path sourceFilePath,
                                 final LinkedHashMap<String, String> annotationOverridesMap,
                                 final FeatureInput<? extends Feature> mainSourceFileAsFeatureInput) {

        super(mainSourceFileAsFeatureInput);

        this.name = name;
        this.version = version;
        this.sourceFilePath = sourceFilePath;

        // Handle the supported field names here:
        supportedFieldNamesAndDefaults = new LinkedHashMap<>();
        supportedFieldNames = new LinkedHashSet<>();
        populateSupportedFieldNamesFromVcfFile();

        // This step has to occur after supported field names and name have been populated.
        supportedFieldMetadata = createFuncotationMetadata(sourceFilePath);

        if (supportedFieldNames.size() == 0) {
            logger.warn("WARNING: VcfFuncotationFactory has nothing to annotate from VCF File: " + sourceFilePath.toUri().toString());
        }
        else {
            // Now check if we have any overrides to take care of:
            this.annotationOverrideMap = new LinkedHashMap<>();
            for ( final Map.Entry<String, String> entry : annotationOverridesMap.entrySet() ) {
                if ( supportedFieldNamesAndDefaults.containsKey(entry.getKey()) ) {
                    annotationOverrideMap.put(entry.getKey(), entry.getValue());
                }
            }
        }
    }

    private FuncotationMetadata createFuncotationMetadata(final Path sourceFilePath) {
        // Read the VCF to just get the header
        try ( final FeatureDataSource<VariantContext> vcfReader = new FeatureDataSource<>(sourceFilePath.toString()) ) {
            final Object header = vcfReader.getHeader();
            if ( ! (header instanceof VCFHeader) ) {
                throw new IllegalArgumentException(sourceFilePath + " does not have a valid VCF header");
            }
            final VCFHeader sourceVcfHeader = (VCFHeader) header;
            final List<VCFInfoHeaderLine> metadataVcfInfoHeaderLines = createFuncotationVcfInfoHeaderLines(sourceVcfHeader);
            return VcfFuncotationMetadata.create(metadataVcfInfoHeaderLines);
        }
    }

    @VisibleForTesting
    List<VCFInfoHeaderLine> createFuncotationVcfInfoHeaderLines(final VCFHeader vcfHeader) {
        final List<VCFInfoHeaderLine> supportedVcfInfoHeaderLines = vcfHeader.getInfoHeaderLines().stream()
                .filter(vcfInfoHeaderLine -> supportedFieldNames.contains(createFinalFieldName(name, vcfInfoHeaderLine.getID())))
                .map(vcfInfoHeaderLine -> copyWithRename(vcfInfoHeaderLine, name))
                .collect(Collectors.toList());

        // Add in the ID field to the meta data:
        final VCFInfoHeaderLine idHeaderLine = new VCFInfoHeaderLine(
                createFinalFieldName(name, "ID"),
                VCFHeaderLineCount.A,
                VCFHeaderLineType.String,
                "ID of the variant from the data source creating this annotation."
                );
        supportedVcfInfoHeaderLines.add( idHeaderLine );

        // Make sure to rename the input VCF field names to the output funcotation field names for this funcotation factory.
        return supportedVcfInfoHeaderLines;
    }

    private static VCFInfoHeaderLine copyWithRename(final VCFInfoHeaderLine vcfInfoHeaderLine, final String name) {
        if (vcfInfoHeaderLine.getCountType() == VCFHeaderLineCount.INTEGER) {
            return new VCFInfoHeaderLine(createFinalFieldName(name, vcfInfoHeaderLine.getID()),
                    vcfInfoHeaderLine.getCount(), vcfInfoHeaderLine.getType(), vcfInfoHeaderLine.getDescription());
        } else {
            return new VCFInfoHeaderLine(createFinalFieldName(name, vcfInfoHeaderLine.getID()),
                    vcfInfoHeaderLine.getCountType(), vcfInfoHeaderLine.getType(), vcfInfoHeaderLine.getDescription());
        }
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public Class<? extends Feature> getAnnotationFeatureClass() {
        return VariantContext.class;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public FuncotatorArgumentDefinitions.DataSourceType getType() {
        return FuncotatorArgumentDefinitions.DataSourceType.VCF;
    }

    @Override
    public LinkedHashSet<String> getSupportedFuncotationFields() {
        return supportedFieldNames;
    }

    @Override
    protected List<Funcotation> createDefaultFuncotationsOnVariant( final VariantContext variant, final ReferenceContext referenceContext ) {
        if ( supportedFieldNames.size() != 0 ) {
            return createDefaultFuncotationsOnVariantHelper(variant, referenceContext, Collections.emptySet());
        }
        else {
            return Collections.emptyList();
        }
    }

    @Override
    /**
     * {@inheritDoc}
     *
     * This is really the only entry point for the Vcf FuncotationFactory.
     *
     * {@link VcfFuncotationFactory} can be used with or without Gencode annotations.
     */
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {

        final List<Funcotation> outputFuncotations = new ArrayList<>();

        // TODO: Caching logic can be refactored and shared in other funcotation factories:  https://github.com/broadinstitute/gatk/issues/4974
        final Triple<VariantContext, ReferenceContext, List<Feature>> cacheKey = createCacheKey(variant, referenceContext, featureList);
        final List<Funcotation> cacheResult = cache.get(cacheKey);
        if (cacheResult != null) {
            cacheHits++;
            return cacheResult;
        }

        // Only create annotations if we have data to annotate:
        if ( supportedFieldNames.size() != 0 ) {

            // Get rid of any null features.
            // By this point we know the feature type is correct, so we cast it:
            final List<VariantContext> funcotationFactoryVariants = featureList.stream().filter(f -> f != null)
                    .map(f -> (VariantContext) f).collect(Collectors.toList());

            // Create a map that will keep the final outputs.  Default it to default funcotations for each alt allele in the
            //  query variant.
            final Map<Allele, Funcotation> outputOrderedMap = new LinkedHashMap<>();

            for ( final VariantContext funcotationFactoryVariant : funcotationFactoryVariants ) {

                // The funcotationFactoryVariants already overlap the query variant in position, now get which
                //  match in ref/alt as well.  And make sure to handle multiallelics in both the query variant and the
                //  funcotation factory variant.
                // matchIndices length will always be the same as the number of alt alleles in the variant (first parameter)
                //  Note that this is not the same length as the funcotationFactoryVariant.
                final int[] matchIndices = GATKVariantContextUtils.matchAllelesOnly(variant, funcotationFactoryVariant);

                for (int i = 0; i < matchIndices.length; i++) {
                    final int matchIndex = matchIndices[i];
                    final Allele queryAltAllele = variant.getAlternateAllele(i);
                    if (matchIndex != -1) {

                        final LinkedHashMap<String, Object> annotations = new LinkedHashMap<>(supportedFieldNamesAndDefaults);

                        for (final Map.Entry<String, Object> entry : funcotationFactoryVariant.getAttributes().entrySet()) {
                            populateAnnotationMap(funcotationFactoryVariant, variant, matchIndex, annotations, entry);
                        }

                        // Add the ID of the variant:
                        annotations.put(createFinalFieldName(name, "ID"), variant.getID());

                        final TableFuncotation newFuncotation = TableFuncotation.create(annotations, queryAltAllele, name, supportedFieldMetadata);
                        outputOrderedMap.merge(queryAltAllele, newFuncotation, VcfFuncotationFactory::mergeDuplicateFuncotationFactoryVariant);
                    }
                }
            }
            variant.getAlternateAlleles().forEach(a -> outputFuncotations.add(outputOrderedMap.computeIfAbsent(a, allele -> createDefaultFuncotation(allele))));
        }
        cacheMisses++;
        cache.put(cacheKey, outputFuncotations);

        // The output number of funcotations should equal to the variant.getAlternateAlleles().size()
        return outputFuncotations;
    }

    /**
     *
     * @param funcotation1 Must have same alt allele and datasource name as other funcotation.
     * @param funcotation2 Must have same alt allele and datasource name as other funcotation.
     * @return a funcotation that contains the merged fields.
     */
    private static Funcotation mergeDuplicateFuncotationFactoryVariant(final Funcotation funcotation1, final Funcotation funcotation2) {
        Utils.validateArg(funcotation1.getAltAllele().equals(funcotation2.getAltAllele()), "Merge called on funcotations that have differing alt alleles.");
        Utils.validateArg(funcotation1.getDataSourceName().equals(funcotation2.getDataSourceName()), "Merge called on funcotations that have differing datasource names.");

        final LinkedHashSet<String> allFieldNames = funcotation1.getFieldNames();
        allFieldNames.addAll(funcotation2.getFieldNames());

        final Map<String, Object> mergedFieldsMap = allFieldNames.stream()
                .collect(Collectors.toMap(f -> f, f -> mergeFuncotationValue(f, funcotation1, funcotation2, VcfFuncotationFactory::renderFieldConflicts)));
        return TableFuncotation.create(mergedFieldsMap, funcotation1.getAltAllele(), funcotation1.getDataSourceName(),
                merge(funcotation1.getMetadata(), funcotation2.getMetadata()));
    }

    /**
     * Given two FuncotationMetadata, create a new one with merged content.  This is basically a set union.
     *
     * No checking is done to make sure that you don't have metadata that are equivalent, but not exact.  All merging
     *  is done on exact match.
     */
    private static FuncotationMetadata merge(final FuncotationMetadata funcotationMetadata1, final FuncotationMetadata funcotationMetadata2) {
        final LinkedHashSet<VCFInfoHeaderLine> rawMetadata = new LinkedHashSet<>(funcotationMetadata1.retrieveAllHeaderInfo());
        rawMetadata.addAll(funcotationMetadata2.retrieveAllHeaderInfo());
        return VcfFuncotationMetadata.create(new ArrayList<>(rawMetadata));
    }

    /**
     *  Return a merged annotation value for the two regions and given annotation name.  Automatically solves conflicts.
     *
     * @param fieldName the annotation to determine.
     * @param funcotation1 first region to merge.
     * @param funcotation2 second region to merge.
     * @param conflictFunction the function to run to solve conflicts.
     * @return string with the new, merged value of the annotation.  Returns {@code null} if the annotation name
     * does not exist in either region.
     */
    private static String mergeFuncotationValue(final String fieldName, final Funcotation funcotation1,
                                                final Funcotation funcotation2, final BiFunction<String, String, String> conflictFunction) {
        final boolean doesRegion1ContainAnnotation = funcotation1.hasField(fieldName);
        final boolean doesRegion2ContainAnnotation = funcotation2.hasField(fieldName);

        if (doesRegion1ContainAnnotation && doesRegion2ContainAnnotation) {

            // Both regions contain an annotation and presumably these are of different values.
            return conflictFunction.apply(funcotation1.getField(fieldName),
                    funcotation2.getField(fieldName));
        } else if (doesRegion1ContainAnnotation) {
            return funcotation1.getField(fieldName);
        } else if (doesRegion2ContainAnnotation) {
            return funcotation2.getField(fieldName);
        }

        return null;
    }

    private static String renderFieldConflicts(final String value1, final String value2) {
        return value1 + DUPLICATE_RECORD_DELIMITER + value2;
    }

    private void populateAnnotationMap(final VariantContext funcotationFactoryVariant, final VariantContext queryVariant, final int funcotationFactoryAltAlleleIndex, final LinkedHashMap<String, Object> annotations, final Map.Entry<String, Object> attributeEntry) {
        final String valueString;
        final String attributeName = attributeEntry.getKey();

        // Handle collections a little differently:
        if (attributeEntry.getValue() instanceof Collection<?>) {
            @SuppressWarnings("unchecked") final Collection<Object> objectList = ((Collection<Object>) attributeEntry.getValue());
            final VCFHeaderLineCount countType = supportedFieldMetadata.retrieveHeaderInfo(createFinalFieldName(this.name, attributeName)).getCountType();

            if (funcotationFactoryVariant.isBiallelic() && queryVariant.isBiallelic()) {
                valueString = objectList.stream().map(Object::toString).collect(Collectors.joining(String.valueOf(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR_CHAR)));
            } else {
                valueString = determineValueStringFromMultiallelicAttributeList(funcotationFactoryAltAlleleIndex, objectList, countType);
            }
        } else {
            valueString = attributeEntry.getValue().toString();
        }

        annotations.put(createFinalFieldName(name, attributeName), valueString);
    }

    /**
     * Given a (possibly) multiallelic attribute value, convert it into a string that (if applicable) can be used for
     *  one alternate allele.
     *
     * @param funcotationFactoryAltAlleleIndex index of the alt allele in the variant from this funcotation factory.
     * @param attributeEntryValues the list of values for an attribute
     * @param countType the count type for the attribute.
     * @return string that can be used in a funcotation.  If the count type is R, then it will have two comma-separated
     *  values.  One for the ref and one for the alt allele being indexed.
     */
    private String determineValueStringFromMultiallelicAttributeList(final int funcotationFactoryAltAlleleIndex, final Collection<Object> attributeEntryValues, final VCFHeaderLineCount countType) {

        switch (countType) {
            case A:
                return attributeEntryValues.toArray()[funcotationFactoryAltAlleleIndex].toString();

            case R:
                final Object referenceAlleleValue = attributeEntryValues.toArray()[0];
                return referenceAlleleValue.toString() + VCFConstants.INFO_FIELD_ARRAY_SEPARATOR_CHAR + attributeEntryValues.toArray()[funcotationFactoryAltAlleleIndex+1].toString();

            default:
                return attributeEntryValues.stream().map(Object::toString).collect(Collectors.joining(String.valueOf(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR_CHAR)));
        }
    }

    @Override
    /**
     * {@inheritDoc}
     * {@link VcfFuncotationFactory} can be used with or without Gencode annotations.
     */
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList, final List<GencodeFuncotation> gencodeFuncotations) {
        return createFuncotationsOnVariant(variant, referenceContext, featureList);
    }


    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    private List<Funcotation> createDefaultFuncotationsOnVariantHelper( final VariantContext variant, final ReferenceContext referenceContext, final Set<Allele> annotatedAltAlleles  ) {

        final List<Funcotation> funcotationList = new ArrayList<>();

        if ( supportedFieldNames.size() != 0 ) {

            final List<Allele> alternateAlleles = variant.getAlternateAlleles();

            for ( final Allele altAllele : alternateAlleles ) {
                if ( !annotatedAltAlleles.contains(altAllele) ) {
                    funcotationList.add(createDefaultFuncotation(altAllele));
                }
            }
        }

        return funcotationList;
    }

    private TableFuncotation createDefaultFuncotation(final Allele altAllele) {
        return TableFuncotation.create(supportedFieldNamesAndDefaults, altAllele, name, supportedFieldMetadata);
    }

    /**
     * Populates {@link VcfFuncotationFactory#supportedFieldNames} and {@link VcfFuncotationFactory#supportedFieldNamesAndDefaults}.
     */
    private void populateSupportedFieldNamesFromVcfFile() {
        final VCFFileReader reader = new VCFFileReader(sourceFilePath.toFile());
        final VCFHeader header = reader.getFileHeader();

        final List<String> infoLineKeys = new ArrayList<>();
        final Map<String, Boolean> infoFieldFlagMap = new HashMap<>();

        // Get our list of keys and sort them:
        for ( final VCFInfoHeaderLine infoLine : header.getInfoHeaderLines() ) {
            infoLineKeys.add(infoLine.getID());
            infoFieldFlagMap.put(infoLine.getID(), infoLine.getType() == VCFHeaderLineType.Flag);
        }
        infoLineKeys.sort(Comparator.naturalOrder());

        // Add our sorted names to the supported list:
        for ( final String key : infoLineKeys ) {
            if ( infoFieldFlagMap.get(key) ) {
                supportedFieldNamesAndDefaults.put(createFinalFieldName(name, key), "false" );
            }
            else {
                supportedFieldNamesAndDefaults.put(createFinalFieldName(name, key), "" );
            }
            supportedFieldNames.add(createFinalFieldName(name, key));
        }

        // Add our ID to the supported fields:
        supportedFieldNamesAndDefaults.put(createFinalFieldName(name, "ID"), "" );
        supportedFieldNames.add(createFinalFieldName(name, "ID"));
    }

    @VisibleForTesting
    static String createFinalFieldName(final String funcotationFactoryName, final String fieldName) {
        return funcotationFactoryName + "_" + fieldName;
    }

    private Triple<VariantContext, ReferenceContext, List<Feature>> createCacheKey(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {
        return Triple.of(variant, referenceContext, featureList);
    }

    @Override
    public void close() {
        logger.info(getName() + " " + getVersion() + " cache hits/total: " + cacheHits + "/" + (cacheMisses + cacheHits));
    }

    //==================================================================================================================
    // Helper Data Types:

    // Modifed from https://docs.oracle.com/javase/7/docs/api/java/util/LinkedHashMap.html#removeEldestEntry(java.util.Map.Entry)
    class LRUCache<K, V> extends LinkedHashMap<K, V> {
        static final long serialVersionUID = 55337L;
        @VisibleForTesting
        static final int MAX_ENTRIES = 20;
        public LRUCache() {
            super(MAX_ENTRIES);
        }

        @Override
        protected boolean removeEldestEntry(Map.Entry<K, V> eldest) {
            return size() > MAX_ENTRIES;
        }
    }
}
