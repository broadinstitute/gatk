package org.broadinstitute.hellbender.tools.walkers.fasta;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import picard.cmdline.programgroups.ReferenceProgramGroup;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;


/**
 * Generate an alternative reference sequence over the specified interval
 *
 * <p>Given a variant callset, this tool replaces the reference bases at variation sites with the bases supplied in the
 * corresponding callset records. Additionally, it allows for one or more "snp-mask" VCFs to set overlapping bases to 'N'.</p>
 *
 * <p>The output format can be partially controlled using the provided command-line arguments.
 * Specify intervals with the usual -L argument to output only the reference bases within your intervals.
 * Overlapping intervals are automatically merged; reference bases for each disjoint interval will be output as a
 * separate fasta sequence (named numerically in order).</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>If there are multiple variants that start at a site, it chooses one of them randomly.</li>
 *     <li>When there are overlapping indels (but with different start positions) only the first will be chosen.</li>
 *     <li>This tool works only for SNPs and for simple indels (but not for things like complex substitutions).</li>
 * </ul>

 * <h3>Input</h3>
 * <p>
 * The reference, requested intervals, and any number of variant files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A FASTA file representing the requested intervals.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk FastaAlternateReferenceMaker \
 *   -R reference.fasta \
 *   -O output.fasta \
 *   -L input.intervals \
 *   -V input.vcf \
 *
 *   [--snp-mask mask.vcf]
 * </pre>
 *
 */
@DocumentedFeature
@CommandLineProgramProperties(summary = "Create an alternative fasta reference by inserting variants from a vcf into an existing reference sequence.",
        oneLineSummary = "Create an alternative reference by combining a fasta with a vcf.",
        programGroup = ReferenceProgramGroup.class)
public class FastaAlternateReferenceMaker extends FastaReferenceMaker {

    //Pre-allocated as a lame optimization
    private static final byte[] NO_BASES = {};
    private static final byte[] N_BYTES = {'N'};
    private static final byte[] n_BYTES = {'n'};
    private static final byte[] A_BYTES = {'A'};
    private static final byte[] a_BYTES = {'a'};
    private static final byte[] G_BYTES = {'G'};
    private static final byte[] g_BYTES = {'g'};
    private static final byte[] C_BYTES = {'C'};
    private static final byte[] c_BYTES = {'c'};
    private static final byte[] T_BYTES = {'T'};
    private static final byte[] t_BYTES = {'t'};

    private static final String EMPTY_BASE = " ";

    /**
     * Variants from this input file are used by this tool to construct an alternate reference.
     */
    @Argument(fullName= StandardArgumentDefinitions.VARIANT_LONG_NAME,
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
            doc="A source of variants to merge with the reference sequence.")
    protected FeatureInput<VariantContext> variants;

    public static final String SNP_MASK_LONG_NAME = "snp-mask";

    /**
     * SNPs from this file are used as a mask (inserting N's in the sequence) when constructing the alternate reference
     */
    @Argument(fullName= SNP_MASK_LONG_NAME, doc="SNP mask VCF file", optional=true)
    protected FeatureInput<VariantContext> snpmask;

    public static final String SNP_MASK_PRIORITY_LONG_NAME = "snp-mask-priority";
    /**
     * Gives priority to a SNP mask over an input VCF for a site. Only has an effect if the --snp-mask argument is used.
     */
    @Argument(fullName= SNP_MASK_PRIORITY_LONG_NAME, doc="Give the snp mask priority over the input VCF.", optional=true)
    protected boolean snpmaskPriority = false;

    public static final String USE_IUPAC_SAMPLE_LONG_NAME = "use-iupac-sample";
    /**
     * This option will generate an error if the specified sample does not exist in the VCF.
     * Non-diploid (or non-called) genotypes are ignored.
     */
    @Argument(fullName= USE_IUPAC_SAMPLE_LONG_NAME, doc = "If specified, heterozygous SNP sites will be output using IUPAC ambiguity codes given the genotypes for this sample", optional=true)
    private String iupacSample = null;

    private int deletionBasesRemaining = 0;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        if (snpmaskPriority && snpmask == null){
            throw new CommandLineException("Cannot specify --" + SNP_MASK_PRIORITY_LONG_NAME + " without " + " --" + SNP_MASK_LONG_NAME);
        }

        if ( iupacSample != null ) {
            final VCFHeader variantHeader = (VCFHeader) getHeaderForFeatures(variants);
            final ArrayList<String> samples = variantHeader.getSampleNamesInOrder();
            if( samples == null || !samples.contains(iupacSample)) {
                throw new UserException.BadInput("the IUPAC sample specified is not present in the provided VCF file");
            }
        }
    }

    @Override
    public void apply(final ReferenceContext ref, final ReadsContext reads, final FeatureContext features) {
        final byte[] bases = handlePosition(ref.getInterval(), ref.getBase(), features);
        final SimpleInterval interval = ref.getInterval();
        if( bases.length == 0){
            advancePosition(interval);
        } else {
            for( byte base : bases) {
                addToReference(interval, base);
            }
        }
    }

    private byte[] handlePosition(SimpleInterval interval, byte base, FeatureContext features) {
        if (deletionBasesRemaining > 0) {
            deletionBasesRemaining--;
            return NO_BASES;
        }

        // If we have a mask at this site, use it
        if ( snpmaskPriority ){
            if (isMasked(features) )
                return N_BYTES;
        }

        // Check to see if we have a called snp
        for ( final VariantContext vc : features.getValues(variants) ) {
            if ( vc.isFiltered() || vc.getStart() != interval.getStart()  )
                continue;

            if ( vc.isSimpleDeletion()) {
                deletionBasesRemaining = vc.getReference().length() - 1;
                // delete the next n bases, not this one
                return baseToByteArray(base);
            } else if ( vc.isSimpleInsertion() || vc.isSNP() ) {
                // Get the first alt allele that is not a spanning deletion. If none present, use the empty allele
                final Optional<Allele> optionalAllele = getFirstConcreteAltAllele(vc.getAlternateAlleles());
                final Allele allele = optionalAllele.orElseGet(() -> Allele.create(EMPTY_BASE, false));
                if ( vc.isSimpleInsertion() ) {
                    return allele.getBases();
                } else {
                    final String iupacBase = (iupacSample != null) ? getIUPACBase(vc.getGenotype(iupacSample)) : allele.toString();
                    return iupacBase.getBytes();
                }
            }
        }

        if ( !snpmaskPriority ){
            if ( isMasked(features)) {
                return N_BYTES;
            }
        }

        // if we got here then we're just ref
        return baseToByteArray(base);
    }

    private byte[] baseToByteArray(byte base) {
        switch(base) {
            case 'A': return A_BYTES;
            case 'a': return a_BYTES;
            case 'C': return C_BYTES;
            case 'c': return c_BYTES;
            case 'G': return G_BYTES;
            case 'g': return g_BYTES;
            case 'T': return T_BYTES;
            case 't': return t_BYTES;
            case 'N': return N_BYTES;
            case 'n': return n_BYTES;
            default: return new byte[]{base};
        }
    }

    /**
     * Get the first non-symbolic, non-spanning deletion, called allele
     * @param altAlleles the alternate alleles
     * @return the first non spanning deletion allele or null
     */
    private Optional<Allele> getFirstConcreteAltAllele(final List<Allele> altAlleles ) {
        return altAlleles.stream()
                .filter(allele -> !allele.isSymbolic())
                .filter(allele -> !GATKVCFConstants.isSpanningDeletion(allele))
                .filter(Allele::isCalled)
                .findFirst();
    }

    /**
     * Mask a SNP (inserting N's in the sequence)
     *
     * @param features the Reference Metadata available at a particular site in the genome
     * @return mask at the locus or null if no SNP at that locus
     */
    private boolean isMasked(final FeatureContext features){
        return features.getValues(snpmask).stream().anyMatch(VariantContext::isSNP);
    }

    /**
     * Returns the IUPAC encoding for the given genotype or the reference base if not possible
     *
     * @param genotype  the genotype to encode
     * @return non-null, non-empty String of bases
     */
    private String getIUPACBase(final Genotype genotype) {
        Utils.nonNull(genotype, () -> "The genotype is null for sample " + iupacSample);

        // If we have a spanning deletion, if both alleles are spanning deletions, use the empty allele. Otherwise, use the allele that is not a
        // spanning deletion.
        if ( genotype.getAlleles().contains(Allele.SPAN_DEL) ) {
            if ( genotype.isHomVar() ) {
                return EMPTY_BASE;
            } else {
                return genotype.getAllele(0).equals(Allele.SPAN_DEL) ? genotype.getAllele(1).getBaseString() : genotype.getAllele(0).getBaseString();
            }
        }

        if ( !genotype.isHet() ) {
            return genotype.getAllele(0).getBaseString();
        }

        final byte allele1 = genotype.getAllele(0).getBases()[0];
        final byte allele2 = genotype.getAllele(1).getBases()[0];
        return new String(new byte[] {BaseUtils.basesToIUPAC(allele1, allele2)});
    }
}