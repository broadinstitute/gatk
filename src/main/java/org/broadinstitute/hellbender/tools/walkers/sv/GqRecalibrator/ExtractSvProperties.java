package org.broadinstitute.hellbender.tools.walkers.sv.GqRecalibrator;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.*;


/**
 * - Extract properties for use in training and/or filtering SV genotypes
 *    - take properties from input VCF
 *    - also compute properties relating to overlaps between variants and genome tracks
 * - Save properties to gzipped TSV files
 *     - properties are either arrays (numVariants x 1) or matrices (numVariants x numSamples)
 *     - ordinal-encode properties that are String or Set<String> (except variant IDs which won't compress)
 *     - directly write out numbers
 * - Also save JSON with summary for each property (name, type, size, and details of any ordinal encoding)
 */






@CommandLineProgramProperties(
        summary = "Extract properties for use in training and/or filtering SV genotypes: 1) take properties from input"+
                "VCF 2) also compute properties relating to overlaps between variants and genome tracks.\n" +
                "Save properties to gzipped TSV files: properties are either arrays (numVariants x 1) or matrices " +
                " (numVariants x numSamples). Ordinal-encode properties that are String or Set<String> (except variant"+
                " IDs which won't compress); directly write out numbers.\n" +
                "Also save JSON with summary for each property (name, type, size, and details of any ordinal encoding)",
        oneLineSummary = "Extract data for training and/or filtering SV genotypes",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ExtractSvProperties extends VariantWalker {
    @Argument(fullName=StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              doc="Output folder that will store gzipped TSVs with extracted properties", optional=true)
    public GATKPath outputFolder = null;

    @Argument(fullName="use-copy-number-calls", doc="If true, attempt to use copy number info if genotype is NO CALL. "+
                       "This also entails keeping track of copy number call quality (RD_GQ) and reconciling that with GQ",
                       optional=true)
    public boolean useCopyNumberCalls = false;

    static final String minSamplesToEstimateAlleleFrequencyKey = "min-samples-to-estimate-allele-frequency";
    @Argument(fullName=minSamplesToEstimateAlleleFrequencyKey, shortName="ms", optional=true,
              doc="If the VCF does not have allele frequency, estimate it from the sample population if there are at least this many samples. Otherwise throw an exception.")
    public int minSamplesToEstimateAlleleFrequency = 100;

    @Argument(fullName="genome-track", shortName="gt", optional=true)
    final List<String> genomeTrackFiles = new ArrayList<>();

    List<TrackOverlapDetector> trackOverlapDetectors = null;

    protected List<String> sampleIds = null; // numSamples list of IDs for samples in VCF, used to keep properties in order

    final Map<String, GzippedTsvWriter> tsvWritersMap = new HashMap<>();

    private int numVariants;
    private int numSamples;
    private Set<String> vcfHeaderIds = null;
    private static final String ID_KEY = "ID";
    private static final String NO_CALL_COUNTS_KEY = "NO_CALL_COUNTS";
    private static final String SVLEN_KEY = "SVLEN";
    private static final String EV_KEY = "EV";
    private static final String CONC_ST_KEY = "CONC_ST";
    private static final String FILTER_KEY = "FILTER";
    private static final String ALGORITHMS_KEY = "ALGORITHMS";
    private static final String EVIDENCE_KEY = "EVIDENCE";
    private static final String STATUS_KEY = "STATUS";
    private static final String NON_REF_GENOTYPE_CONCORDANCE_KEY = "NON_REF_GENOTYPE_CONCORDANCE";
    private static final String VAR_PPV_KEY = "VAR_PPV";
    private static final String VAR_SENSITIVITY_KEY = "VAR_SENSITIVITY";
    private static final String TRUTH_AF_KEY = "TRUTH_AF";
    private static final String MULTIALLELIC_FILTER = "MULTIALLELIC";
    private static final String CALL_QUALITY_KEY = "CALL_QUALITY";
    private static final String IS_COPY_NUMBER_CALL_KEY = "IS_CN_CALL";

    private static final String PE_GQ_KEY = "PE_GQ";
    private static final String RD_CN_KEY = "RD_CN";
    private static final String RD_GQ_KEY = "RD_GQ";
    private static final String SR_GQ_KEY = "SR_GQ";
    private static final short MISSING_GQ_VAL = -1;

    int getNumVariants() { return numVariants; }
    int getNumSamples() { return numSamples; }
    int getNumProperties() { return tsvWritersMap.size(); }

    // stats on tool actions, input/output VCFs
    private int numVarGenotypes;
    private int numRefGenotypes;
    private int numNoCallGenotypes;

    @Override
    public void onTraversalStart() {
        ensureOutputFolderExists();
        // get sample IDs to take properties in order
        sampleIds = getHeaderForVariants().getGenotypeSamples();
        numSamples = sampleIds.size();
        // get valid property IDs from header
        vcfHeaderIds = Stream.concat(
            getHeaderForVariants().getInfoHeaderLines().stream().map(VCFCompoundHeaderLine::getID),
            getHeaderForVariants().getFormatHeaderLines().stream().map(VCFCompoundHeaderLine::getID)
        ).collect(Collectors.toSet());
        // load genome tracks
        trackOverlapDetectors = genomeTrackFiles.stream().map(TrackOverlapDetector::new).collect(Collectors.toList());
        // initialize counts
        numVariants = 0;
        numVarGenotypes = 0;
        numRefGenotypes = 0;
        numNoCallGenotypes = 0;
        // initialize array variables
        sampleAlleleCounts = new byte[numSamples];
        sampleNoCallCounts = new byte[numSamples];
        isCopyNumberCall = useCopyNumberCalls ? new boolean[numSamples] : null;
        callQuality = useCopyNumberCalls ? new short[numSamples] : null;
    }

    void ensureOutputFolderExists() {
        final File outputFolderFile = new File(outputFolder.toPath().toUri());
        if(!outputFolderFile.exists()) {
            if(!outputFolderFile.mkdirs()) {
                throw new RuntimeException("Unable to create outputFolder: " + outputFolderFile);
            }
        }
    }

    private void getTrackProperties(final TrackOverlapDetector trackOverlapDetector,
                                    final VariantContext variantContext) {
        final List<SimpleInterval> genomeTrackOverlapLocations =
            SimpleSvInterval.streamFrom(variantContext)
                    .flatMap(SimpleSvInterval::streamGenomeTrackOverlapLocations)
                    .distinct()
                    .collect(Collectors.toList());

        final double[] overlaps;
        try {
            if (trackOverlapDetector.hasOther()) { // This track has paired intervals
                overlaps = genomeTrackOverlapLocations.stream()
                    .mapToDouble(location -> trackOverlapDetector.getPrimaryOverlapFraction(location) +
                                             trackOverlapDetector.getOtherOverlapfraction(location))
                    .toArray();
                // check for spanning by streaming over all possible pairs
                final boolean spans = IntStream.range(0, genomeTrackOverlapLocations.size() - 1)
                    .anyMatch(
                        i -> IntStream.range(i + 1, genomeTrackOverlapLocations.size()).anyMatch(
                            j -> trackOverlapDetector.spansPrimaryAndOther(
                                genomeTrackOverlapLocations.get(i),
                                genomeTrackOverlapLocations.get(j)
                            )
                        )
                    );
                getTsvWriter(trackOverlapDetector.getName() + "_spans").append(spans);
            } else {  // this track has simple intervals
                overlaps = genomeTrackOverlapLocations.stream()
                        .mapToDouble(trackOverlapDetector::getPrimaryOverlapFraction)
                        .toArray();
            }
            getTsvWriter(trackOverlapDetector.getName() + "_min").append(
                Arrays.stream(overlaps).min().orElse(0.0)
            );
            getTsvWriter(trackOverlapDetector.getName() + "_max").append(
                    Arrays.stream(overlaps).max().orElse(0.0)
            );
        } catch(RuntimeException runtimeException) {
            throw new RuntimeException(
                "Error getting overlap for " + variantContext.getID() + " and " + trackOverlapDetector.getName(),
                    runtimeException
            );
        }
    }

    private static short getGenotypeAttributeAsShort(final Genotype genotype, final String key, Short defaultValue) {
        if(key.equals(VCFConstants.GENOTYPE_QUALITY_KEY)) {
            return (short)genotype.getGQ();
        }
        Object x = genotype.getExtendedAttribute(key);
        if (x == null ||
            x == VCFConstants.MISSING_VALUE_v4 ||
            (x instanceof Character && (char)x == VCFConstants.NO_CALL_ALLELE)) {
            if(defaultValue == null) {
                throw new IllegalArgumentException("Genotype is missing value of " + key);
            } else {
                return defaultValue;
            }
        }
        if ( x instanceof Short ) return (Short)x;
        try {
            //noinspection ConstantConditions
            return Short.parseShort((String) x); // throws an exception if this isn't a string
        } catch(ClassCastException classCastException) {
            throw new GATKException("Unable to extract value " + x + " for key " + key);
        }
    }

    @SuppressWarnings("SameParameterValue")
    private short[] getGenotypeAttributeAsShort(final Iterable<Genotype> sampleGenotypes, final String attributeKey,
                                                final Short missingAttributeValue) {
        final short[] values = getShortPropertyArray(attributeKey);
        int index = 0;
        for(final Genotype genotype : sampleGenotypes) {
            values[index] = getGenotypeAttributeAsShort(genotype, attributeKey, missingAttributeValue);
            ++index;
        }
        return values;
    }

    private static Set<String> getVcfPropertyAsStringSet(Object rawProperty) {
        if ( rawProperty == null ||
                rawProperty == VCFConstants.EMPTY_ALLELE ||
                (rawProperty instanceof Character && (char)rawProperty == VCFConstants.NO_CALL_ALLELE) ) {
            return Collections.emptySet();
        } else {
            // Throws an exception if this isn't a String
            if (rawProperty instanceof String) {
                // Remove any weird brackets, and split by commas
                return Arrays.stream(((String)rawProperty)
                        .replaceAll("[\\[\\] ]", "")
                        .split(",")
                ).collect(Collectors.toSet());
            } else {
                throw new IllegalArgumentException("Value was not a String");
            }
        }
    }

    private static boolean isMissing(final List<String> values) {
        return values == null || values.isEmpty() || (
                values.size() == 1 && (values.get(0) == null || values.get(0).equals(VCFConstants.EMPTY_ALLELE))
        );
    }

    private static Set<String> getInfoAttributeAsStringSet(final VariantContext variantContext, final String key) {
        try {
            final List<String> values = variantContext.getAttributeAsStringList(key, null);
            return isMissing(values) ? Collections.emptySet() : new HashSet<>(values);
        } catch(IllegalArgumentException illegalArgumentException) {
            throw new IllegalArgumentException(
                "Getting " + key + " from INFO field: " + illegalArgumentException.getMessage(),
                illegalArgumentException
            );
        }
    }

    private static Set<String> getGenotypeAttributeAsStringSet(final Genotype genotype, final String key) {
        Object rawProperty = genotype.getExtendedAttribute(key);
        try {
            return getVcfPropertyAsStringSet(rawProperty);
        } catch(IllegalArgumentException illegalArgumentException) {
            throw new IllegalArgumentException(
                "Getting " + key + " for sample " + genotype.getSampleName() + ": " + illegalArgumentException.getMessage(),
                illegalArgumentException
           );
        }
    }

    @SuppressWarnings("SameParameterValue")
    private Set<String>[] getGenotypeAttributeAsStringSet(final Iterable<Genotype> sampleGenotypes,
                                                          final String key) {
        final Set<String>[] values = getStringSetPropertyArray(key);

        int index = 0;
        for(final Genotype genotype : sampleGenotypes) {
            values[index] = getGenotypeAttributeAsStringSet(genotype, key);
            ++index;
        }

        return values;
    }

    GzippedTsvWriter getTsvWriter(final String propertyName) {
        GzippedTsvWriter tsvWriter = tsvWritersMap.getOrDefault(propertyName, null);
        if(tsvWriter == null) {
            tsvWriter = new GzippedTsvWriter(propertyName, outputFolder.toPath());
            tsvWritersMap.put(propertyName, tsvWriter);
        }
        return tsvWriter;
    }

    private class AlleleCountsGetter {
        final byte nonRefCounts;
        final byte numCalledCounts;
        final byte noCallCounts;
        final short callQuality;
        final boolean isCopyNumberCall;

        AlleleCountsGetter(final Genotype genotype) {
            byte numAlleles = 0;
            byte noCallCounts = 0;
            byte nonRefCounts = 0;
            short callQuality = MISSING_GQ_VAL;
            for(final Allele allele : genotype.getAlleles()) {
                // Do normal counting of alleles from genotype field
                ++numAlleles;
                if(allele.isNoCall()) {
                    ++noCallCounts;
                } else if(!allele.isReference()) {
                    ++nonRefCounts;
                }
            }

            boolean isCopyNumberCall = false;
            if(useCopyNumberCalls && genotype.isNoCall()) {
                // Check if maybe this is a copy number call
                // Get the value from this key, if it's available
                final short rd_cn = getGenotypeAttributeAsShort(genotype, RD_CN_KEY, (short)-1);
                if(rd_cn != -1) { // got copy number info, use it
                    nonRefCounts = (byte)(rd_cn == 2 ? 0 : 1);
                    noCallCounts = 0;
                    numAlleles = 2;
                    isCopyNumberCall = true;
                    callQuality = getGenotypeAttributeAsShort(genotype, RD_GQ_KEY, MISSING_GQ_VAL);
                }
            }
            if(!isCopyNumberCall) {
                callQuality = getGenotypeAttributeAsShort(genotype, VCFConstants.GENOTYPE_QUALITY_KEY, MISSING_GQ_VAL);
            }
            this.nonRefCounts = nonRefCounts;
            this.numCalledCounts = (byte)(numAlleles - noCallCounts);
            this.noCallCounts = noCallCounts;
            this.callQuality = callQuality;
            this.isCopyNumberCall = isCopyNumberCall;
        }
    }

    /**
     * Wrap actual code of apply with a try/catch call so that problematic variants can be identified by exception
     */
    @Override
    public void apply(VariantContext variantContext, ReadsContext readsContext, ReferenceContext ref, FeatureContext featureContext) {
        try {
            wrapped_apply(variantContext);
        } catch(Exception exception) {
            throw new IllegalArgumentException(
                "Error processing variant " + variantContext.getID() + ": " + exception.getMessage(),
                exception
            );
        }
    }

    byte[] sampleAlleleCounts = null;
    byte[] sampleNoCallCounts = null;
    boolean[] isCopyNumberCall = null;
    short[] callQuality = null;
    final Map<String, short[]> shortPropertiesArrayMap = new HashMap<>();
    final Map<String, Set<String>[]> stringSetPropertiesArrayMap = new HashMap<>();

    short[] getShortPropertyArray(final String propertyName) {
        final short[] allocatedArray = shortPropertiesArrayMap.getOrDefault(propertyName, null);
        if(allocatedArray == null) {
            final short[] newArray = new short[numSamples];
            shortPropertiesArrayMap.put(propertyName, newArray);
            return newArray;
        } else {
            return allocatedArray;
        }
    }

    Set<String>[] getStringSetPropertyArray(final String propertyName) {
        final Set<String>[] allocatedArray = stringSetPropertiesArrayMap.getOrDefault(propertyName, null);
        if(allocatedArray == null) {
            @SuppressWarnings("unchecked")
            final Set<String>[] newArray = (Set<String>[]) new Set<?>[numSamples];
            stringSetPropertiesArrayMap.put(propertyName, newArray);
            return newArray;
        } else {
            return allocatedArray;
        }
    }


    /**
     * Accumulate properties for variant matrix, and allele counts, genotype quality for trio tensors
     */
    private void wrapped_apply(VariantContext variantContext) throws IOException {
        /////////////////////////////////////// Process allele counts and quality //////////////////////////////////////
        int numCalledAlleles = 0;
        int numNonRefAlleles = 0;
        int numVariantInputVar = 0;
        int numVariantInputNoCall = 0;
        int numVariantInputRef = 0;
        int sampleIndex = 0;
        for(final Genotype genotype : variantContext.getGenotypesOrderedBy(sampleIds)) {
            final AlleleCountsGetter alleleCountsGetter = new AlleleCountsGetter(genotype);
            sampleAlleleCounts[sampleIndex] = alleleCountsGetter.nonRefCounts;
            sampleNoCallCounts[sampleIndex] = alleleCountsGetter.noCallCounts;
            if(useCopyNumberCalls) {
                isCopyNumberCall[sampleIndex] = alleleCountsGetter.isCopyNumberCall;
                callQuality[sampleIndex] = alleleCountsGetter.callQuality;
            }
            ++sampleIndex;

            if(alleleCountsGetter.nonRefCounts > 0) {
                ++numVariantInputVar;
            } else if(alleleCountsGetter.noCallCounts > 0) {
                ++numVariantInputNoCall;
            } else {
                ++numVariantInputRef;
            }
            numNonRefAlleles += alleleCountsGetter.nonRefCounts;
            numCalledAlleles += alleleCountsGetter.numCalledCounts;
        }
        numVarGenotypes += numVariantInputVar;
        numNoCallGenotypes += numVariantInputNoCall;
        numRefGenotypes += numVariantInputRef;

        ++numVariants;
        /////////////////////////////////////////// Get values from FORMAT /////////////////////////////////////////////
        getTsvWriter(VCFConstants.ALLELE_COUNT_KEY).append(sampleAlleleCounts);
        getTsvWriter(NO_CALL_COUNTS_KEY).append(sampleNoCallCounts);
        if(useCopyNumberCalls) {
            // if using copy number calls, add the reconciled call quality, and whether this is a copy number call
            getTsvWriter(IS_COPY_NUMBER_CALL_KEY).append(isCopyNumberCall);
            getTsvWriter(CALL_QUALITY_KEY).append(callQuality);
        }
        final Iterable<Genotype> sampleGenotypes = variantContext.getGenotypesOrderedBy(sampleIds);
        // Append various qualities
        Stream.of(VCFConstants.GENOTYPE_QUALITY_KEY, PE_GQ_KEY, SR_GQ_KEY, RD_GQ_KEY)
            .filter(vcfHeaderIds::contains)
            .forEach(
                genotypeAttribute -> getTsvWriter(genotypeAttribute).append(
                    getGenotypeAttributeAsShort(sampleGenotypes, genotypeAttribute, MISSING_GQ_VAL)
                )
            );
        Stream.of(EV_KEY, CONC_ST_KEY)
            .filter(vcfHeaderIds::contains)
            .forEach(
                genotypeAttribute -> getTsvWriter(genotypeAttribute).append(
                    getGenotypeAttributeAsStringSet(sampleGenotypes, genotypeAttribute)
                )
            );
        //////////////////////////////////////// Get or estimate allele frequency //////////////////////////////////////
        float alleleFrequency = (float)variantContext.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, -1.0);
        if(alleleFrequency <= 0) {
            if(variantContext.getNSamples() <= minSamplesToEstimateAlleleFrequency) {
                throw new GATKException("VCF does not have " + VCFConstants.ALLELE_FREQUENCY_KEY + " annotated or enough samples to estimate it ("
                                        + minSamplesToEstimateAlleleFrequencyKey + "=" + minSamplesToEstimateAlleleFrequency + " but there are "
                                        + variantContext.getNSamples() + " samples)");
            }
            // VCF not annotated with allele frequency, guess it from allele counts. If somehow we have a variant with
            // no called alleles, just set alleleFrequency to 0. There's nothing to do with it anyway, so it hardly
            // matters
            alleleFrequency = numCalledAlleles > 0 ? numNonRefAlleles / (float) numCalledAlleles : 0F;
        }
        getTsvWriter(VCFConstants.ALLELE_FREQUENCY_KEY).append(alleleFrequency);
        ///////////////////////////////////////////////// Get FILTER ///////////////////////////////////////////////////
        getTsvWriter(FILTER_KEY).append(variantContext.getFilters());
        ////////////////////////////////////////////// Get INFO fields /////////////////////////////////////////////////
        // SVTYPE is mandatory
        final String svType = variantContext.getAttributeAsString(VCFConstants.SVTYPE, null);
        if(svType == null) {
            throw new GATKException("Missing " + VCFConstants.SVTYPE + " for variant " + variantContext.getID());
        }
        getTsvWriter(VCFConstants.SVTYPE).append(svType);
        // Some variants don't have annotated SVLEN, in which case infer from location
        getTsvWriter(SVLEN_KEY).append(SimpleSvInterval.getOrInferSvLen(variantContext));
        // get INFO fields that are sets of strings:
        Stream.of(ALGORITHMS_KEY, EVIDENCE_KEY, STATUS_KEY)
            .filter(vcfHeaderIds::contains)
            .forEach(
                infoAttribute -> getTsvWriter(infoAttribute).append(
                    getInfoAttributeAsStringSet(variantContext, infoAttribute)
                )
            );
        // add INFO fields that are floats
        Stream.of(NON_REF_GENOTYPE_CONCORDANCE_KEY, VAR_PPV_KEY, VAR_SENSITIVITY_KEY, TRUTH_AF_KEY)
            .filter(vcfHeaderIds::contains)
            .forEach(
                infoAttribute -> getTsvWriter(infoAttribute).append(
                    (float)variantContext.getAttributeAsDouble(infoAttribute, Double.NaN)
                )
            );
        /////////////////////////////////////////  Get genome track overlaps ///////////////////////////////////////////
        for(final TrackOverlapDetector trackOverlapDetector : trackOverlapDetectors) {
            getTrackProperties(trackOverlapDetector, variantContext);
        }

        getTsvWriter(ID_KEY).append(variantContext.getID(), false);
    }

    @Override
    public Object onTraversalSuccess() {
        if(getNumVariants() == 0) {
            throw new GATKException("No variants contained in vcf: " + drivingVariantFile);
        }

        printPropertiesDebugInfo();
        for(final GzippedTsvWriter tsvWritersMap : tsvWritersMap.values()) {
            tsvWritersMap.close();
        }
        final Path summaryPath = outputFolder.toPath().resolve("properties_summary.json");
        GzippedTsvWriter.savePropertiesSummaryJson(tsvWritersMap.values(), summaryPath);
        final Path sampleIdsPath = outputFolder.toPath().resolve("sample_ids.list");
        saveSampleIds(sampleIdsPath);
        return null;
    }

    void printPropertiesDebugInfo() {
        System.out.println("########################################");
        System.out.println("numVariants: " + getNumVariants());
        System.out.println("numSamples: " + getNumSamples());
        System.out.println("numProperties: " + getNumProperties());
        System.out.format("Input VCF had %.1f variants per sample.\n", numVarGenotypes / (double)numSamples);
        System.out.format("Genotype break-down: %d ref, %d non-ref, %d no-call\n", numRefGenotypes, numVarGenotypes, numNoCallGenotypes);
        printPropertiesValuesSummary();

        for(final GzippedTsvWriter tsvWriter : tsvWritersMap.values()) {
            if(tsvWriter.hasLabels()) {
                System.out.format("%s:\n", tsvWriter.getName());
                int idx = 0;
                for(final String label : tsvWriter.getAllLabels()) {
                    System.out.format("%d\t%s\n", idx, label);
                    ++idx;
                }
            }
        }
        System.out.println("########################################");
    }

    void printPropertiesValuesSummary() {
        final int nameWidth = FastMath.max(
                "propertyName".length(),
                tsvWritersMap.keySet().stream().mapToInt(String::length).max().orElse(0)
        );
        final int typeWidth = FastMath.max(
                "type".length(),
                tsvWritersMap.values().stream().mapToInt(writer -> writer.getTypeName().length()).max().orElse(0)
        );
        final String nameFormat = "%" + nameWidth + "s";
        final String typeFormat = "%" + typeWidth + "s";
        final String headerFormat = String.join(
                "\t", "%5s", nameFormat, typeFormat, "%8s", "%8s"
        ) + "\n";
        System.out.format(headerFormat, "index", "propertyName", "type", "nRows", "nColumns");
        final String dataFormat = String.join(
                "\t", "%5d", nameFormat, typeFormat, "%8d", "%8d"
        ) + "\n";
        int idx = 0;
        for(final GzippedTsvWriter writer : tsvWritersMap.values()) {
            System.out.format(dataFormat,
                    idx, writer.getName(), writer.getTypeName(), writer.getNumRows(), writer.getNumColumns());
            ++idx;
        }
    }

    void saveSampleIds(final Path sampleIdsPath) {

        try(final OutputStream outputStream = Files.newOutputStream(sampleIdsPath)) {
            for(final String sampleId : sampleIds) {
                outputStream.write((sampleId + "\n").getBytes(StandardCharsets.UTF_8));
            }
        } catch(IOException ioException) {
            throw new RuntimeException("Error writing GzippedTsvWriter properties summary", ioException);
        }
    }
}
