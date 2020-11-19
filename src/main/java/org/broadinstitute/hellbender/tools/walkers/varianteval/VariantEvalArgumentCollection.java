package org.broadinstitute.hellbender.tools.walkers.varianteval;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;

/**
 * The collection of arguments for VariantEval
 */
public class VariantEvalArgumentCollection {

    // constants
    public final static String ALL_SAMPLE_NAME = "all";
    public final static String ALL_FAMILY_NAME = "all";

    /**
     * The variant file(s) to evaluate.
     */
    @Argument(fullName="eval", shortName = "eval", doc="Input evaluation file(s)", optional=false)
    public List<FeatureInput<VariantContext>> evals;

    /**
     * The variant file(s) to compare against.
     */
    @Argument(fullName = StandardArgumentDefinitions.COMPARISON_LONG_NAME, shortName = StandardArgumentDefinitions.COMPARISON_SHORT_NAME, doc="Input comparison file(s)", optional=true)
    public List<FeatureInput<VariantContext>> compsProvided = new ArrayList<>();
    public List<FeatureInput<VariantContext>> comps = new ArrayList<>();

    /**
     * dbSNP comparison VCF.  By default, the dbSNP file is used to specify the set of "known" variants.
     * Other sets can be specified with the -known-name (--known_names) argument.
     */
    @ArgumentCollection
    public DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    /**
     * Some analyses want to count overlap not with dbSNP (which is in general very open) but
     * actually want to itemize their overlap specifically with a set of gold standard sites
     * such as HapMap, OMNI, or the gold standard indels.  This argument provides a mechanism
     * for communicating which file to use
     */
    @Argument(fullName="gold-standard", shortName = "gold", doc="Evaluations that count calls at sites of true variation (e.g., indel calls) will use this argument as their gold standard for comparison", optional=true)
    public FeatureInput<VariantContext> goldStandard = null;

    // Partitioning the data arguments
    @Argument(shortName="select", doc="One or more stratifications to use when evaluating the data", optional=true)
    public ArrayList<String> selectExps = new ArrayList<String>();

    @Argument(shortName="select-name", doc="Names to use for the list of stratifications (must be a 1-to-1 mapping)", optional=true)
    public ArrayList<String> selectNames = new ArrayList<String>();

    @Argument(fullName="sample", shortName="sn", doc="Derive eval and comp contexts using only these sample genotypes, when genotypes are available in the original context", optional=true)
    public Set<String> sampleExpressions = new TreeSet<>();

    @Argument(fullName = StandardArgumentDefinitions.PEDIGREE_FILE_LONG_NAME, shortName = StandardArgumentDefinitions.PEDIGREE_FILE_SHORT_NAME, doc="Pedigree file for determining the population \"founders\"", optional=true)
    public GATKPath pedigreeFile;

    /**
     * List of feature tracks to be used for specifying "known" variants other than dbSNP.
     */
    @Argument(shortName="known-name", doc="Name of feature bindings containing variant sites that should be treated as known when splitting eval features into known and novel subsets", optional=true)
    public Set<String> knownNames = new HashSet<String>();

    // Stratification arguments
    @Argument(fullName="stratification-module", shortName="ST", doc="One or more specific stratification modules to apply to the eval track(s) (in addition to the standard stratifications, unless -noS is specified)", optional=true)
    public List<String> stratificationsToUse = new ArrayList<>();

    @Argument(fullName="do-not-use-all-standard-stratifications", shortName="no-st", doc="Do not use the standard stratification modules by default (instead, only those that are specified with the -S option)", optional=true)
    public Boolean noStandardStratifications = false;

    /**
     * See the -list argument to view available modules.
     */
    @Argument(fullName="eval-module", shortName="EV", doc="One or more specific eval modules to apply to the eval track(s) (in addition to the standard modules, unless -no-ev is specified)", optional=true)
    public List<String> modulesToUse = new ArrayList<>();

    @Argument(fullName="do-not-use-all-standard-modules", shortName="no-ev", doc="Do not use the standard modules by default (instead, only those that are specified with the -EV option)", optional=true)
    public Boolean noStandardModules = false;

    @Argument(fullName="min-phase-quality", shortName="mpq", doc="Minimum phasing quality", optional=true)
    public double minPhaseQuality = 10.0;

    @Argument(shortName="mvq", fullName="mendelian-violation-qual-threshold", doc="Minimum genotype QUAL score for each trio member required to accept a site as a violation. Default is 50.", optional=true)
    public double mendelianViolationQualThreshold = 50;

    @Argument(shortName="ploidy", fullName="sample-ploidy", doc="Per-sample ploidy (number of chromosomes per sample)", optional=true)
    public int ploidy = GATKVariantContextUtils.DEFAULT_PLOIDY;

    @Argument(fullName="ancestral-alignments", shortName="aa", doc="Fasta file with ancestral alleles", optional=true)
    public File ancestralAlignmentsFile = null;

    @Argument(fullName="require-strict-allele-match", shortName="strict", doc="If provided only comp and eval tracks with exactly matching reference and alternate alleles will be counted as overlapping", optional=true)
    public boolean requireStrictAlleleMatch = false;

    @Argument(fullName="keep-ac0", shortName="keep-ac0", doc="If provided, modules that track polymorphic sites will not require that a site have AC > 0 when the input eval has genotypes", optional=true)
    public boolean keepSitesWithAC0 = false;

    /**
     * If true, VariantEval will treat -eval 1 -eval 2 as separate tracks from the same underlying
     * variant set, and evaluate the union of the results.  Useful when you want to do -eval chr1.vcf -eval chr2.vcf etc.
     */
    @Argument(fullName="merge-evals", shortName="merge-evals", doc="If provided, all -eval tracks will be merged into a single eval track", optional=true)
    public boolean mergeEvals = false;

    /**
     * File containing tribble-readable features for the IntervalStratificiation
     */
    @Argument(fullName="strat-intervals", shortName="strat-intervals", doc="File containing tribble-readable features for the IntervalStratificiation", optional=true)
    public FeatureInput<Feature> intervalsFile = null;

    /**
     * File containing tribble-readable features containing known CNVs.  For use with VariantSummary table.
     */
    @Argument(fullName="known-cnvs", shortName="known-cnvs", doc="File containing tribble-readable features describing a known list of copy number variants", optional=true)
    public FeatureInput<Feature> knownCNVsFile = null;

    @Hidden
    @Argument(fullName="num-samples", doc="If provided, modules that track polymorphic sites will not require that a site have AC > 0 when the input eval has genotypes", optional=true)
    public int numSamplesFromArgument = 0;

    public List<FeatureInput<VariantContext>> getFeatureInputsForDrivingVariants() {
        List<FeatureInput<VariantContext>> ret = new ArrayList<>(evals);

        if (dbsnp.dbsnp != null) {
            ret.add(dbsnp.dbsnp);
        }

        ret.addAll(compsProvided);

        return ret;
    }

    public FeatureInput<Feature> getKnownCNVsFile() {
        return knownCNVsFile;
    }

    public List<FeatureInput<VariantContext>> getEvals() {
        return Collections.unmodifiableList(evals);
    }

    public List<FeatureInput<VariantContext>> getComps() {
        return Collections.unmodifiableList(comps);
    }

    public boolean ignoreAC0Sites() {
        return ! keepSitesWithAC0;
    }

    public FeatureInput<VariantContext> getGoldStand() {
        return goldStandard;
    }

    public List<FeatureInput<VariantContext>> getCompsProvided() {
        return Collections.unmodifiableList(compsProvided);
    }

    public boolean isMergeEvals() {
        return mergeEvals;
    }

    public int getPloidy() {
        return ploidy;
    }

    public FeatureInput<Feature> getIntervalsFile() {
        return intervalsFile;
    }

    public double getMendelianViolationQualThreshold() {
        return mendelianViolationQualThreshold;
    }
}
