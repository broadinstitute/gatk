package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStatificationEngine;
import org.broadinstitute.hellbender.tools.sv.stratify.SVStratificationEngineArgumentsCollection;
import org.broadinstitute.hellbender.utils.*;

import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>Stratifies structural variants into separate groups according to the following customizable criteria:
 * <ul>
 *     <li>SV type</li>
 *     <li>Size range</li>
 *     <li>Reference context</li>
 * </ul>
 * Users should provide a stratification configuration .tsv file (tab-delimited table) with the following column
 * header on the first line:
 * <ol>
 *     <li>NAME</li>
 *     <li>SVTYPE</li>
 *     <li>MIN_SIZE</li>
 *     <li>MAX_SIZE</li>
 *     <li>CONTEXT</li>
 * </ol>
 * </p>
 * <p>For example:</p>
 * <table border="1">
 *   <tr>
 *      <td>NAME</td><td>SVTYPE</td><td>MIN_SIZE</td><td>MAX_SIZE</td><td>CONTEXT</td>
 *   </tr>
 *   <tr>
 *      <td>DEL_large_SD</td><td>DEL</td><td>5000</td><td>-1</td><td>SD</td>
 *   </tr>
 *   <tr>
 *      <td>DUP_large_SD</td><td>DUP</td><td>5000</td><td>-1</td><td>SD</td>
 *   </tr>
 *   <tr>
 *      <td>DEL_small_SR_RM</td><td>DEL</td><td>-1</td><td>5000</td><td>SR,RM</td>
 *   </tr>
 *   <tr>
 *      <td>DUP_small_SR_RM</td><td>DUP</td><td>-1</td><td>5000</td><td>SR,RM</td>
 *   </tr>
 *   <tr>
 *      <td>INS_SR</td><td>INS</td><td>-1</td><td>-1</td><td>SR</td>
 *   </tr>
 * </table>
 * <p>
 * The "NAME" column is an arbitrary identifier, "SVTYPE" is the class of variant (DEL, DUP, INS, etc.), MIN_SIZE in an
 * inclusive size lower-bound, MAX_SIZE is an exclusive size upper-bound, and CONTEXT is a comma-delimited list of
 * reference contexts defined using the {@link SVStratificationEngineArgumentsCollection#contextFileList} and
 * {@link SVStratificationEngineArgumentsCollection#contextNameList} parameters. For example,
 * </p>
 * <pre>
 *     gatk GroupedSVCluster \
 *       --context-name RM \
 *       --context-intervals repeatmasker.bed \
 *       --context-name SD \
 *       --context-intervals segmental_duplications.bed \
 *       --context-name SR \
 *       --context-intervals simple_repeats.bed \
 *       ...
 * </pre>
 * <p>
 * The MIN_SIZE, MAX_SIZE, and CONTEXT columns may contain null values {"-1", "", "NULL", "NA"}. Null MIN_SIZE and
 * MAX_SIZE correspond to negative and positive infinity, respectively, and a null CONTEXT value means that variants
 * will not be matched based on context. Variants with undefined SVLEN will only match if both MIN_SIZE and MAX_SIZE
 * are null.
 * </p>
 *
 * <p>The {@link SVStratificationEngineArgumentsCollection#overlapFraction},
 * {@link SVStratificationEngineArgumentsCollection#numBreakpointOverlaps}, and
 * {@link SVStratificationEngineArgumentsCollection#numBreakpointOverlapsInterchrom} can be used to modify the overlap
 * criteria for assigning variants to each group based on overlap with the given reference context intervals. By
 * default, only one endpoint of the variant needs to lie in a context interval in order to match. INS variants are
 * treated as single points and only {@link SVStratificationEngineArgumentsCollection#numBreakpointOverlaps} is used,
 * ignoring {@link SVStratificationEngineArgumentsCollection#overlapFraction}. Similarly, CTX and BND variant
 * overlap is only defined by {@link SVStratificationEngineArgumentsCollection#numBreakpointOverlapsInterchrom}.
 * </p>
 *
 * <p>By default, each stratification group must be mutually exclusive, meaning that any given SV can only belong to
 * one group. An error is thrown if the tool encounters a variant that meets the criteria for more than one group.
 * This restriction can be overridden with the {@link SVStratify#ALLOW_MULTIPLE_MATCHES_LONG_NAME} argument.
 * Furthermore, SVs that do not fall into any of the groups will be assigned to the
 * {@link SVStratify#DEFAULT_STRATUM} group.</p>
 *
 * <p>The tool generates a set of VCFs as output, with each VCF containing the records of each group. In addition,
 * records are annotated with their respective strata names in the {@link GATKSVVCFConstants#STRATUM_INFO_KEY} INFO
 * field.</p>
 *
 * <p>This tool accepts multiple VCF inputs with no restrictions on site or sample overlap.</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         One or more SV VCFs
 *     </li>
 *     <li>
 *         Stratification configuration TSV file
 *     </li>
 *     <li>
 *         Reference dictionary
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Set of stratified VCFs
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk GroupedSVCluster \
 *       -V variants.vcf.gz \
 *       -O out \
 *       --sequence-dictionary reference.dict \
 *       --context-name RM \
 *       --context-intervals repeatmasker.bed \
 *       --stratify-config strata.tsv
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Splits VCFs by SV type, size, and reference context",
        oneLineSummary = "Splits VCFs by SV type, size, and reference context",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class SVStratify extends MultiVariantWalker {

    public static final String ALLOW_MULTIPLE_MATCHES_LONG_NAME = "allow-multiple-matches";

    // Default output group name for unmatched records
    public static final String DEFAULT_STRATUM = "default";

    @Argument(
            doc = "Output directory",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputPath;

    @Argument(
            doc = "Prefix for output filenames",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME
    )
    private String outputPrefix;

    @ArgumentCollection
    private final SVStratificationEngineArgumentsCollection stratArgs = new SVStratificationEngineArgumentsCollection();

    @Argument(
            doc = "Do not enforce mutual exclusivity for each stratification group",
            fullName = ALLOW_MULTIPLE_MATCHES_LONG_NAME
    )
    private boolean allowMultipleMatches = false;

    protected SAMSequenceDictionary dictionary;
    protected Map<String, VariantContextWriter> writers;
    protected SVStatificationEngine engine;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        dictionary = getMasterSequenceDictionary();
        Utils.validateArg(dictionary != null, "Reference dictionary is required; please specify with --" +
                StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME);
        engine = loadStratificationConfig(stratArgs, dictionary);
        logger.debug("Loaded stratification groups:");
        for (final SVStatificationEngine.Stratum s : engine.getStrata()) {
            logger.debug(s);
        }
        initializeWriters();
    }

    protected void createGroupWriter(final String name) {
        final String filename = outputPrefix + "." + name + ".vcf.gz";
        final Path path = outputPath.toPath().resolve(filename);
        final VariantContextWriter writer = createVCFWriter(path);
        final VCFHeader header = new VCFHeader(getHeaderForVariants());
        addStratifyMetadata(header);
        writer.writeHeader(header);
        if (writers.containsKey(name)) {
            throw new GATKException.ShouldNeverReachHereException("Stratification name already exists: " + name);
        }
        writers.put(name, writer);
    }

    public static void addStratifyMetadata(final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.STRATUM_INFO_KEY, 1,
                VCFHeaderLineType.String, "Stratum ID"));
    }

    protected void initializeWriters() {
        writers = new HashMap<>();
        createGroupWriter(DEFAULT_STRATUM);
        for (final SVStatificationEngine.Stratum s : engine.getStrata()) {
            createGroupWriter(s.getName());
        }
    }

    /**
     * Reusable method for loading the stratification configuration table. See tool doc for the expected format.
     */
    public static SVStatificationEngine loadStratificationConfig(final SVStratificationEngineArgumentsCollection args,
                                                                 final SAMSequenceDictionary dictionary) {
        Utils.validateArg(args.contextNameList.size() == args.contextFileList.size(), "Arguments --" +
                SVStratificationEngineArgumentsCollection.CONTEXT_NAME_FILE_LONG_NAME + " and --" + SVStratificationEngineArgumentsCollection.CONTEXT_INTERVAL_FILE_LONG_NAME +
                " must be specified the same number of times.");
        final Map<String, List<Locatable>> map = new HashMap<>();
        final Iterator<String> nameIterator = args.contextNameList.iterator();
        final Iterator<GATKPath> pathIterator = args.contextFileList.iterator();
        final GenomeLocParser genomeLocParser = new GenomeLocParser(dictionary);
        while (nameIterator.hasNext() && pathIterator.hasNext()) {
            final String name = nameIterator.next();
            final GATKPath path = pathIterator.next();
            final GenomeLocSortedSet genomeLocs = IntervalUtils.loadIntervals(Collections.singletonList(path.toString()), IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, genomeLocParser);
            final List<Locatable> intervals = Collections.unmodifiableList(genomeLocs.toList());
            if (map.containsKey(name)) {
                throw new UserException.BadInput("Duplicate context name was specified: " + name);
            }
            map.put(name, intervals);
        }
        final SVStatificationEngine engine = SVStatificationEngine.create(map, args.configFile, dictionary);
        if (engine.getStrata().stream().anyMatch(s -> s.getName().equals(DEFAULT_STRATUM))) {
            throw new UserException.BadInput("Stratification configuration contains entry with reserved " +
                    "ID \"" + DEFAULT_STRATUM + "\"");
        }
        return engine;
    }

    @Override
    public void closeTool() {
        for (final VariantContextWriter writer : writers.values()) {
            writer.close();
        }
        super.closeTool();
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        // Save a ton of compute by not copying genotypes into the new record
        final VariantContext variantNoGenotypes = new VariantContextBuilder(variant).genotypes(Collections.emptyList()).make();
        final SVCallRecord record = SVCallRecordUtils.create(variantNoGenotypes, dictionary);
        final Collection<SVStatificationEngine.Stratum> stratifications = engine.getMatches(record,
                stratArgs.overlapFraction, stratArgs.numBreakpointOverlaps, stratArgs.numBreakpointOverlapsInterchrom);
        final VariantContextBuilder builder = new VariantContextBuilder(variant);
        if (stratifications.isEmpty()) {
            writers.get(DEFAULT_STRATUM).add(builder.attribute(GATKSVVCFConstants.STRATUM_INFO_KEY, DEFAULT_STRATUM).make());
        } else {
            if (!allowMultipleMatches && stratifications.size() > 1) {
                final String matchesString = String.join(", ", stratifications.stream().map(SVStatificationEngine.Stratum::getName).collect(Collectors.toList()));
                throw new GATKException("Record " + record.getId() + " matched multiple groups: " + matchesString + ". Bypass this error using the --" + ALLOW_MULTIPLE_MATCHES_LONG_NAME + " argument");
            }
            for (final SVStatificationEngine.Stratum stratum : stratifications) {
                writers.get(stratum.getName()).add(builder.attribute(GATKSVVCFConstants.STRATUM_INFO_KEY, stratum.getName()).make());
            }
        }
    }
}
