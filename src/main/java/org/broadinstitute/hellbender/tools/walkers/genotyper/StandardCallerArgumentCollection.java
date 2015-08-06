package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.collections4.map.DefaultedMap;
import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.Hidden;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculatorImplementation;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.Collections;
import java.util.Map;

/**
 * This is pulled out so that every caller isn't exposed to the arguments from every other caller.
 */
public class StandardCallerArgumentCollection {

    @ArgumentCollection
    public GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    @Argument(fullName = "genotyping_mode", shortName = "gt_mode", doc = "Specifies how to determine the alternate alleles to use for genotyping", optional=true)
    public GenotypingOutputMode genotypingOutputMode = GenotypingOutputMode.DISCOVERY;

    /**
     * When the UnifiedGenotyper is put into GENOTYPE_GIVEN_ALLELES mode it will genotype the samples using only the alleles provide in this rod binding
     */
    @Argument(fullName="alleles", shortName = "alleles", doc="The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES", optional=true)
    public FeatureInput<VariantContext> alleles;

    /**
     * If this fraction is greater is than zero, the caller will aggressively attempt to remove contamination through biased down-sampling of reads.
     * Basically, it will ignore the contamination fraction of reads for each alternate allele.  So if the pileup contains N total bases, then we
     * will try to remove (N * contamination fraction) bases for each alternate allele.
     */
    @Argument(fullName = "contamination_fraction_to_filter", shortName = "contamination", doc = "Fraction of contamination in sequencing data (for all samples) to aggressively remove", optional=true)
    public double CONTAMINATION_FRACTION = DEFAULT_CONTAMINATION_FRACTION;
    public static final double DEFAULT_CONTAMINATION_FRACTION = 0.0;

    /**
     *  This argument specifies a file with two columns "sample" and "contamination" specifying the contamination level for those samples.
     *  Samples that do not appear in this file will be processed with CONTAMINATION_FRACTION.
     **/
    @Advanced
    @Argument(fullName = "contamination_fraction_per_sample_file", shortName = "contaminationFile", doc = "Tab-separated File containing fraction of contamination in sequencing data (per sample) to aggressively remove. Format should be \"<SampleID><TAB><Contamination>\" (Contamination is double) per line; No header.", optional = true)
    public File CONTAMINATION_FRACTION_FILE = null;

    private DefaultedMap<String,Double> sampleContamination;

    /**
     * Returns the sample contamination or CONTAMINATION_FRACTION if no contamination level was specified for this sample.
     */
    public Double getSampleContamination(final String sampleId){
        Utils.nonNull(sampleId);
        if (sampleContamination == null){
            setSampleContamination(new DefaultedMap<>(CONTAMINATION_FRACTION));//default to empty map
        }
        return sampleContamination.get(sampleId);
    }

    public void setSampleContamination(final DefaultedMap<String, Double> sampleContamination) {
        this.sampleContamination = new DefaultedMap<>(CONTAMINATION_FRACTION);  //NOTE: a bit weird because it ignores the default from the argument and uses ours
        this.sampleContamination.putAll(sampleContamination);                   //make a copy to be safe
    }

    /**
     * Controls the model used to calculate the probability that a site is variant plus the various sample genotypes in the data at a given locus.
     */
    @Hidden
    @Argument(fullName = "p_nonref_model", shortName = "pnrm", doc = "Non-reference probability calculation model to employ", optional = true)
    public AFCalculatorImplementation requestedAlleleFrequencyCalculationModel;

    @Hidden
    @Argument(shortName = "logExactCalls", doc="x", optional=true)
    public File exactCallsLog = null;

    @Argument(fullName = "output_mode", shortName = "out_mode", doc = "Specifies which type of calls we should output", optional = true)
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
    @Argument(fullName = "allSitePLs", shortName = "allSitePLs", doc = "Annotate all sites with PLs", optional = true)
    public boolean annotateAllSitesWithPLs = false;

    /**
     * Creates a Standard caller argument collection with default values.
     */
    public StandardCallerArgumentCollection() { }

    /**
     * Holds a modifiers mask that identifies those fields that cannot be copied between
     * StandardCallerArgumentCollections.
     */
    private final int UNCOPYABLE_MODIFIER_MASK = Modifier.PRIVATE | Modifier.STATIC | Modifier.FINAL;
}
