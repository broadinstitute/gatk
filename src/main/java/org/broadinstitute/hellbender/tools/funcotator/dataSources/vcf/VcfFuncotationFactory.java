package org.broadinstitute.hellbender.tools.funcotator.dataSources.vcf;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.funcotator.DataSourceFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.FuncotatorArgumentDefinitions;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.gencode.GencodeFuncotation;

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
     * {@link VcfFuncotationFactory} can be used with or without Gencode annotations.
     */
    protected List<Funcotation> createFuncotationsOnVariant(final VariantContext variant, final ReferenceContext referenceContext, final List<Feature> featureList) {

        final List<Funcotation> outputFuncotations = new ArrayList<>();

        // Only create annotations if we have data to annotate:
        if ( supportedFieldNames.size() != 0 ) {

            final List<Allele> alternateAlleles = variant.getAlternateAlleles();

            // Create a set to put our annotated Alternate alleles in.
            // We'll use this to determine if the alt allele has been annotated.
            final Set<Allele> annotatedAltAlleles = new HashSet<>(alternateAlleles.size());

            if ( !featureList.isEmpty() ) {
                for ( final Feature feature : featureList ) {

                    if ( feature != null ) {

                        // By this point we know the feature type is correct, so we cast it:
                        final VariantContext variantFeature = (VariantContext) feature;

                        // Now we create one funcotation for each Alternate allele:
                        for ( final Allele altAllele : alternateAlleles ) {
                            if (!(variantFeature.hasAlternateAllele(altAllele) && variantFeature.getReference().equals(variant.getReference()))) {
                                continue;
                            }
                            // Add all Info keys/values to a copy of our default map:
                            final LinkedHashMap<String, Object> annotations = new LinkedHashMap<>(supportedFieldNamesAndDefaults);
                            for ( final Map.Entry<String, Object> entry : variantFeature.getAttributes().entrySet() ) {

                                final String valueString;

                                // Handle collections a little differently:
                                if (entry.getValue() instanceof Collection<?>) {
                                    @SuppressWarnings("unchecked") final Collection<Object> objectList = ((Collection<Object>) entry.getValue());
                                    valueString = objectList.stream().map(Object::toString).collect(Collectors.joining(","));
                                } else {
                                    valueString = entry.getValue().toString();
                                }

                                annotations.put(name + "_" + entry.getKey(), valueString);
                            }

                            // Add our funcotation to the funcotation list:
                            outputFuncotations.add(new TableFuncotation(annotations, altAllele, name));
                            annotatedAltAlleles.add(altAllele);
                        }
                    }
                }
            }

            // If we didn't add funcotations for an allele, we should add in blank funcotations to that allele for each field that can be produced
            // by this VcfFuncotationFactory:
            if ( annotatedAltAlleles.size() != alternateAlleles.size() ) {
                outputFuncotations.addAll( createDefaultFuncotationsOnVariantHelper(variant, referenceContext, annotatedAltAlleles) );
            }
        }

        return outputFuncotations;
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
                    funcotationList.add(new TableFuncotation(supportedFieldNamesAndDefaults, altAllele, name));
                }
            }
        }

        return funcotationList;
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
                supportedFieldNamesAndDefaults.put( name + "_" + key, "false" );
            }
            else {
                supportedFieldNamesAndDefaults.put( name + "_" + key, "" );
            }
            supportedFieldNames.add( name + "_" + key );
        }
    }

    //==================================================================================================================
    // Helper Data Types:

}
