package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections4.map.DefaultedMap;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculatorImplementation;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.Serializable;
import java.util.Collections;
import java.util.Map;

/**
 * This is pulled out so that every caller isn't exposed to the arguments from every other caller.
 */
public class StandardCallerArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Copies the values from other into this StandardCallerArgumentCollection
     *
     * @param other StandardCallerArgumentCollection from which to copy values
     */
    public void copyStandardCallerArgsFrom( final StandardCallerArgumentCollection other ) {
        Utils.nonNull(other);

        this.genotypeArgs = new GenotypeCalculationArgumentCollection(other.genotypeArgs);
        this.genotypingOutputMode = other.genotypingOutputMode;
        this.alleles = other.alleles; // FeatureInputs are immutable outside of the engine, so this shallow copy is safe
        this.CONTAMINATION_FRACTION = other.CONTAMINATION_FRACTION;
        this.CONTAMINATION_FRACTION_FILE = other.CONTAMINATION_FRACTION_FILE != null ? new File(other.CONTAMINATION_FRACTION_FILE.getAbsolutePath()) : null;
        if ( other.sampleContamination != null ) {
            setSampleContamination(other.sampleContamination);
        }
        this.requestedAlleleFrequencyCalculationModel = other.requestedAlleleFrequencyCalculationModel;
        this.exactCallsLog = other.exactCallsLog != null ? new File(other.exactCallsLog.getAbsolutePath()) : null;
        this.outputMode = other.outputMode;
        this.annotateAllSitesWithPLs = other.annotateAllSitesWithPLs;
    }

    @ArgumentCollection
    public GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    @Argument(fullName = "genotyping-mode", doc = "Specifies how to determine the alternate alleles to use for genotyping", optional=true)
    public GenotypingOutputMode genotypingOutputMode = GenotypingOutputMode.DISCOVERY;

    /**
     * When the caller is put into GENOTYPE_GIVEN_ALLELES mode it will genotype the samples using only the alleles provide in this rod binding
     */
    @Argument(fullName="alleles", doc="The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES", optional=true)
    public FeatureInput<VariantContext> alleles;

    /**
     * When set to true an when in GENOTYPE_GIVEN_ALLELES mode all given alleles, even filtered ones, are genotyped
     */
    @Advanced
    @Argument(fullName = "genotype-filtered-alleles", doc = "Whether to genotype all given alleles, even filtered ones, --genotyping_mode is GENOTYPE_GIVEN_ALLELES", optional = true)
    public boolean genotypeFilteredAlleles = false;

    /**
     * If this fraction is greater is than zero, the caller will aggressively attempt to remove contamination through biased down-sampling of reads.
     * Basically, it will ignore the contamination fraction of reads for each alternate allele.  So if the pileup contains N total bases, then we
     * will try to remove (N * contamination fraction) bases for each alternate allele.
     */
    @Argument(fullName = "contamination-fraction-to-filter", shortName = "contamination", doc = "Fraction of contamination in sequencing data (for all samples) to aggressively remove", optional=true)
    public double CONTAMINATION_FRACTION = DEFAULT_CONTAMINATION_FRACTION;
    public static final double DEFAULT_CONTAMINATION_FRACTION = 0.0;

    /**
     *  This argument specifies a file with two columns "sample" and "contamination" specifying the contamination level for those samples.
     *  Samples that do not appear in this file will be processed with CONTAMINATION_FRACTION.
     **/
    @Advanced
    @Argument(fullName = "contamination-fraction-per-sample-file", shortName = "contamination-file", doc = "Tab-separated File containing fraction of contamination in sequencing data (per sample) to aggressively remove. Format should be \"<SampleID><TAB><Contamination>\" (Contamination is double) per line; No header.", optional = true)
    public File CONTAMINATION_FRACTION_FILE = null;

    private DefaultedMap<String,Double> sampleContamination;
    private boolean mapHasContaminationSet = false;

    /**
     * Returns true if there is some sample contamination present, false otherwise.
     * @return {@code true} iff there is some sample contamination
     */
    public boolean isSampleContaminationPresent() {
        return contaminationFractionIsSet(CONTAMINATION_FRACTION) || mapHasContaminationSet;
    }
    
    /**
     * Returns an unmodifiable view of the map of SampleId -> contamination.
     * 
     * The returned map will return a default value equal to the configured
     * {@link #CONTAMINATION_FRACTION} for samples whose contamination is not
     * explicitly set.
     */
    public Map<String,Double> getSampleContamination() {
        if (sampleContamination == null) {
            setSampleContamination(Collections.emptyMap()); // default to empty map
        }
        return Collections.unmodifiableMap(sampleContamination);
    }

    /**
     * Set the sample contamination map using the provided map. The resulting map will have
     * its default value for unknown keys set equal to {@link #CONTAMINATION_FRACTION}, regardless
     * of any default value set in the provided map (if it's a DefaultedMap).
     *
     * @param sampleContamination Map of sample to contamination fraction with which to initialize our
     *                            sample contamination map. Replaces any existing values in our map.
     *                            The resulting map will have {@link #CONTAMINATION_FRACTION} as the default
     *                            value for unknown keys, regardless of any default set in the provided map.
     */
    public void setSampleContamination(final Map<String, Double> sampleContamination) {
        this.sampleContamination = new DefaultedMap<>(CONTAMINATION_FRACTION);  //NOTE: a bit weird because it ignores the default from the argument and uses ours
        this.sampleContamination.putAll(sampleContamination);                   //make a copy to be safe

        this.mapHasContaminationSet = contaminationIsPresentInMap(this.sampleContamination);
    }

    /**
     * @param fraction double value to test
     * @return True if fraction represents non-zero contamination, otherwise false
     */
    private boolean contaminationFractionIsSet(final double fraction) {
        return ! Double.isNaN(fraction) && fraction > 0.0;
    }

    /**
     * Given a map of sample to contamination fraction, determines whether any samples have
     * a non-zero contamination fraction set.
     *
     * @param contaminationMap sample -> contamination fraction map to test
     * @return true if at least one sample has a non-zero contamination fraction, otherwise false
     */
    private boolean contaminationIsPresentInMap(final Map<String, Double> contaminationMap) {
        if ( contaminationMap == null ) {
            return false;
        }

        for ( final Map.Entry<String,Double> mapEntry : contaminationMap.entrySet() ) {
            if ( contaminationFractionIsSet(mapEntry.getValue()) ) {
                return true;
            }
        }

        return false;
    }

    /**
     * Controls the model used to calculate the probability that a site is variant plus the various sample genotypes in the data at a given locus.
     */
    @Hidden
    @Argument(fullName = "p-nonref-model", doc = "Non-reference probability calculation model to employ", optional = true)
    public AFCalculatorImplementation requestedAlleleFrequencyCalculationModel;

    @Hidden
    @Argument(shortName = "log-exact-calls", optional=true)
    public File exactCallsLog = null;

    @Argument(fullName = "output-mode", doc = "Specifies which type of calls we should output", optional = true)
    public OutputMode outputMode = OutputMode.EMIT_VARIANTS_ONLY;

    /**
     * Advanced, experimental argument: if SNP likelihood model is specified, and if EMIT_ALL_SITES output mode is set, when we set this argument then we will also emit PLs at all sites.
     * This will give a measure of reference confidence and a measure of which alt alleles are more plausible (if any).
     * WARNINGS:
     * - This feature will inflate VCF file size considerably.
     * - All SNP ALT alleles will be emitted with corresponding 10 PL values.
     * - An error will be emitted if EMIT_ALL_SITES is not set, or if anything other than diploid SNP model is used
     */
    @Advanced
    @Argument(fullName = "all-site-pls", doc = "Annotate all sites with PLs", optional = true)
    public boolean annotateAllSitesWithPLs = false;
}
