package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Set of arguments for the {@link HaplotypeCallerEngine}
 */
public class HaplotypeCallerArgumentCollection extends AssemblyBasedCallerArgumentCollection implements Serializable{
    private static final long serialVersionUID = 1L;

    /**
     * Which annotations to add to the output VCF file. The single value 'none' removes the default annotations.
     * See the VariantAnnotator -list argument to view available annotations.
     */
    @Advanced
    @Argument(fullName = "annotation", shortName = "A", doc = "One or more specific annotations to apply to variant calls", optional = true)
    public List<String> annotationsToUse = new ArrayList<>();

    /**
     * Which groups of annotations to add to the output VCF file. The single value 'none' removes the default group. See
     * the VariantAnnotator -list argument to view available groups. Note that this usage is not recommended because
     * it obscures the specific requirements of individual annotations. Any requirements that are not met (e.g. failing
     * to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Argument(fullName = "group", shortName = "G", doc = "One or more classes/groups of annotations to apply to variant calls", optional = true)
    public List<String> annotationGroupsToUse = new ArrayList<>(Arrays.asList(new String[]{"StandardAnnotation", "StandardHCAnnotation"}));

    /**
     * You can use this argument to specify that HC should process a single sample out of a multisample BAM file. This
     * is especially useful if your samples are all in the same file but you need to run them individually through HC
     * in -ERC GVC mode (which is the recommended usage). Note that the name is case-sensitive.
     */
    @Argument(fullName = "sample_name", shortName = "sn", doc = "Name of single sample to use from a multi-sample bam", optional = true)
    public String sampleNameToUse = null;

    // -----------------------------------------------------------------------------------------------
    // general advanced arguments to control haplotype caller behavior
    // -----------------------------------------------------------------------------------------------

    /**
     * When HC is run in reference confidence mode with banding compression enabled (-ERC GVCF), homozygous-reference
     * sites are compressed into bands of similar genotype quality (GQ) that are emitted as a single VCF record. See
     * the FAQ documentation for more details about the GVCF format.
     *
     * This argument allows you to set the GQ bands. HC expects a list of strictly increasing GQ values
     * that will act as exclusive upper bounds for the GQ bands. To pass multiple values,
     * you provide them one by one with the argument, as in `-GQB 10 -GQB 20 -GQB 30` and so on
     * (this would set the GQ bands to be `[0, 10), [10, 20), [20, 30)` and so on, for example).
     * Note that GQ values are capped at 99 in the GATK, so values must be integers in [1, 100].
     * If the last value is strictly less than 100, the last GQ band will start at that value (inclusive)
     * and end at 100 (exclusive).
     */
    @Advanced
    @Argument(fullName = "GVCFGQBands", shortName = "GQB", doc= "Exclusive upper bounds for reference confidence GQ bands " +
            "(must be in [1, 100] and specified in increasing order)", optional = true)
    public List<Integer> GVCFGQBands = new ArrayList<>(70);
    {
            for (int i=1; i<=60; ++i) {
                GVCFGQBands.add(i);
            }
            GVCFGQBands.add(70); GVCFGQBands.add(80); GVCFGQBands.add(90); GVCFGQBands.add(99);
    };

    /**
     * This parameter determines the maximum size of an indel considered as potentially segregating in the
     * reference model.  It is used to eliminate reads from being indel informative at a site, and determines
     * by that mechanism the certainty in the reference base.  Conceptually, setting this parameter to
     * X means that each informative read is consistent with any indel of size < X being present at a specific
     * position in the genome, given its alignment to the reference.
     */
    @Advanced
    @Argument(fullName = "indelSizeToEliminateInRefModel", shortName = "ERCIS", doc = "The size of an indel to check for in the reference model", optional = true)
    public int indelSizeToEliminateInRefModel = 10;


    @Advanced
    @Argument(fullName = "useAllelesTrigger", shortName = "allelesTrigger", doc = "Use additional trigger on variants found in an external alleles file", optional = true)
    public boolean USE_ALLELES_TRIGGER = false;
}
