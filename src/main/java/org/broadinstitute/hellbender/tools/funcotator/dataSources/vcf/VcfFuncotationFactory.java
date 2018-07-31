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
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.FuncotationMetadata;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.nio.file.Path;
import java.util.*;
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

    @VisibleForTesting
    int cacheHits = 0;
    @VisibleForTesting
    int cacheMisses = 0;

    //==================================================================================================================
    // Constructors:

    public VcfFuncotationFactory(final String name, final String version, final Path sourceFilePath) {
        this(name, version, sourceFilePath, new LinkedHashMap<>());
    }

    public VcfFuncotationFactory(final String name, final String version, final Path sourceFilePath, final LinkedHashMap<String, String> annotationOverridesMap) {
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
                .collect(Collectors.toList());

        // Make sure to rename the input VCF field names to the output funcotation field names for this funcotation factory.
        return supportedVcfInfoHeaderLines.stream()
                .map(vcfInfoHeaderLine -> copyWithRename(vcfInfoHeaderLine, name))
                .collect(Collectors.toList());
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
    protected Class<? extends Feature> getAnnotationFeatureClass() {
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
            variant.getAlternateAlleles().forEach(a -> outputOrderedMap.put(a, createDefaultFuncotation(a)));

            // TODO: What happens if there is a duplicate pos,ref,alt in the datasource?  See (https://github.com/broadinstitute/gatk/issues/4972)
            for ( final VariantContext funcotationFactoryVariant : funcotationFactoryVariants ) {

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

                        outputOrderedMap.put(queryAltAllele, TableFuncotation.create(annotations, queryAltAllele, name, supportedFieldMetadata));
                    }
                }
            }
            outputOrderedMap.keySet().forEach(f -> outputFuncotations.add(outputOrderedMap.get(f)));
        }
        cacheMisses++;
        cache.put(cacheKey, outputFuncotations);

        // The output number of funcotations should equal to the variant.getAlternateAlleles().size()
        return outputFuncotations;
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
