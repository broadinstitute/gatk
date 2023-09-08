package org.broadinstitute.hellbender.tools.walkers.varianteval;

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.Utils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.File;
import java.util.Collections;
import java.util.List;

/**
 *
 *
 * <p>
 * Given a variant callset, it is common to calculate various quality control metrics. These metrics include the number of
 * raw or filtered SNP counts; ratio of transition mutations to transversions; concordance of a particular sample's calls
 * to a genotyping chip; number of s per sample; etc. Furthermore, it is often useful to stratify these metrics
 * by various criteria like functional class (missense, nonsense, silent), whether the site is CpG site, the amino acid
 * degeneracy of the site, etc. VariantEval facilitates these calculations in two ways: by providing several built-in
 * evaluation and stratification modules, and by providing a framework that permits the easy development of new evaluation
 * and stratification modules.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more variant sets to evaluate plus any number of comparison sets.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Evaluation tables detailing the results of the eval modules which were applied.
 * For example:
 * </p>
 * <pre>
 * output.eval.grp:
 * ##:GATKReport.v0.1 CountVariants : Counts different classes of variants in the sample
 * CountVariants  CompFeatureInput  CpG      EvalFeatureInput  JexlExpression  Novelty  nProcessedLoci  nCalledLoci  nRefLoci  nVariantLoci  variantRate ...
 * CountVariants  dbsnp             CpG      eval              none            all      65900028        135770       0         135770        0.00206024  ...
 * CountVariants  dbsnp             CpG      eval              none            known    65900028        47068        0         47068         0.00071423  ...
 * CountVariants  dbsnp             CpG      eval              none            novel    65900028        88702        0         88702         0.00134601  ...
 * CountVariants  dbsnp             all      eval              none            all      65900028        330818       0         330818        0.00502000  ...
 * CountVariants  dbsnp             all      eval              none            known    65900028        120685       0         120685        0.00183133  ...
 * CountVariants  dbsnp             all      eval              none            novel    65900028        210133       0         210133        0.00318866  ...
 * CountVariants  dbsnp             non_CpG  eval              none            all      65900028        195048       0         195048        0.00295976  ...
 * CountVariants  dbsnp             non_CpG  eval              none            known    65900028        73617        0         73617         0.00111710  ...
 * CountVariants  dbsnp             non_CpG  eval              none            novel    65900028        121431       0         121431        0.00184265  ...
 * ...
 * </pre>
 * </p>
 *
 * <h3>Usage examples</h3>
 * <pre>
 * gatk VariantEval \
 *   -R reference.fasta \
 *   -O output.eval.grp \
 *   --eval set1:set1.vcf \
 *   --eval set2:set2.vcf \
 *   [--comp comp.vcf]
 * </pre>
 *
 * Count Mendelian violations for each family in a callset with multiple families (and provided pedigree)
 * <pre>
 * gatk VariantEval \
 *   -R reference.fasta \
 *   -O output.MVs.byFamily.table \
 *   --eval multiFamilyCallset.vcf \
 *   -no-ev -noST \
 *   -ST Family \
 *   -EV MendelianViolationEvaluator
 * </pre>
 *
 * <h3>Caveat</h3>
 *
 * <p>Some stratifications and evaluators are incompatible with each other due to their respective memory requirements,
 * such as AlleleCount and VariantSummary, or Sample and VariantSummary. If you specify such a combination, the program
 * will output an error message and ask you to disable one of these options. We do not currently provide an exhaustive
 * list of incompatible combinations, so we recommend trying out combinations that you are interested in on a dummy
 * command line, to rapidly ascertain whether it will work or not.</p>
 *
 */
@CommandLineProgramProperties(
        summary = "Given a variant callset, it is common to calculate various quality control metrics. These metrics include the number of " +
                "raw or filtered SNP counts; ratio of transition mutations to transversions; concordance of a particular sample's calls " +
                "to a genotyping chip; number of singletons per sample; etc. Furthermore, it is often useful to stratify these metrics " +
                "by various criteria like functional class (missense, nonsense, silent), whether the site is CpG site, the amino acid " +
                "degeneracy of the site, etc. VariantEval facilitates these calculations in two ways: by providing several built-in " +
                "evaluation and stratification modules, and by providing a framework that permits the easy development of new evaluation " +
                "and stratification modules.",
        oneLineSummary = "General-purpose tool for variant evaluation (% in dbSNP, genotype concordance, Ti/Tv ratios, and a lot more)",
        programGroup = VariantEvaluationProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class VariantEval extends MultiVariantWalkerGroupedOnStart {
    protected VariantEvalEngine engine;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    protected File outFile;

    @ArgumentCollection
    protected VariantEvalArgumentCollection variantEvalArgs = new VariantEvalArgumentCollection();

    /**
     * Note that the --list argument requires a fully resolved and correct command-line to work.
     */
    @Argument(fullName="list", shortName="ls", doc="List the available eval modules and exit", optional=true)
    protected Boolean LIST = false;

    @Override
    protected MultiVariantInputArgumentCollection getMultiVariantInputArgumentCollection() {
        return new MultiVariantInputArgumentCollection() {
            private static final long serialVersionUID = 1L;

            @Override
            public List<GATKPath> getDrivingVariantPaths() {
                //driving variants will be determined by initializeDrivingVariants()
                return Collections.emptyList();
            }
        };
    }

    @Override
    protected void initializeDrivingVariants() {
        getDrivingVariantsFeatureInputs().addAll(variantEvalArgs.getFeatureInputsForDrivingVariants());

        super.initializeDrivingVariants();
    }

    @Override
    public void onTraversalStart() {
        // Just list the modules, and exit quickly.
        if (LIST) { listModulesAndExit(); }

        Utils.nonNull(outFile);
        IOUtil.assertFileIsWritable(outFile);

        engine = new VariantEvalEngine(variantEvalArgs, this.features, getTraversalIntervals(), getSequenceDictionaryForDrivingVariants(), getSamplesForVariants());
    }

    /**
     * List all of the available evaluation modules, then exit successfully
     */
    public void listModulesAndExit() {
        logger.info("Available stratification modules:");
        logger.info("(Standard modules are starred)");

        for (String name: VariantEvalEngine.getStratifierClasses().keySet()) {

            logger.info("\t" + name + (VariantEvalEngine.getRequiredStratificationNames().contains(name) || VariantEvalEngine.getStandardStratificationNames().contains(name) ? "*" : ""));
        }
        logger.info("");

        logger.info("Available evaluation modules:");
        logger.info("(Standard modules are starred)");
        for (String veName : VariantEvalEngine.getEvaluatorClasses().keySet()) {
            logger.info("\t" + veName + (VariantEvalEngine.getStandardEvaluatorNames().contains(veName) ? "*" : ""));
        }
        logger.info("");

        System.exit(0);
    }

    @Override
    public void apply(final List<VariantContext> variantContexts, final ReferenceContext referenceContext, final List<ReadsContext> readsContexts) {
        engine.apply(variantContexts, referenceContext);
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Finalizing variant report");

        engine.finalizeReport(outFile);

        return null;
    }

    @Override
    public List<? extends CommandLinePluginDescriptor<?>> getPluginDescriptors() {
        // Reason this override is needed:
        // Because getDefaultReadFilters() is nearly always included (see GATKTool), this pulls in SampleReadFilter,
        // which has an argument named "sample", which conflicts with the local argument.  we get an exception like:
        // org.broadinstitute.barclay.argparser.CommandLineException$CommandLineParserInternalException: [sample, sample] has already been used.
        // at org.broadinstitute.barclay.argparser.CommandLineArgumentParser.handleArgumentAnnotation(CommandLineArgumentParser.java:1002)
        // therefore override the method and return an empty list.
        return Collections.emptyList();
    }
}
