package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
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
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SVLocatable;
import org.broadinstitute.hellbender.tools.sv.cluster.*;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>Clusters structural variants based on coordinates, event type, and supporting algorithms. Primary use cases include:</p>
 * <ul>
 *     <li>
 *         Clustering SVs produced by multiple callers, based on interval overlap, breakpoint proximity, and sample overlap.
 *     </li>
 *     <li>
 *         Merging multiple SV VCFs with disjoint sets of samples and/or variants.
 *     </li>
 *     <li>
 *         Defragmentation of copy number variants produced with depth-based callers.
 *     </li>
 * </ul>
 *
 * <p>The tool generates a new VCF with clusters collapsed into single representative records. By default, a MEMBERS
 * field is generated that lists the input variant IDs contained in that record's cluster.</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         One or more SV VCFs
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Clustered VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk CrossReferenceSVGenotypes \
 *       -V variants1.vcf.gz \
 *       -V variants2.vcf.gz \
 *       --test-vcf test.vcf.gz \
 *       -O annotated.vcf.gz
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Cross references structural variant genotypes",
        oneLineSummary = "Cross references structural variant genotypes",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public final class CrossReferenceSVGenotypes extends MultiVariantWalker {

    public static final String TEST_VCF_LONG_NAME = "test-vcf";

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFile;

    @Argument(
            doc = "Test VCF to annotate using cross-referenced genotypes provided with one or more -"
                    + StandardArgumentDefinitions.VARIANT_SHORT_NAME + " arguments",
            fullName = TEST_VCF_LONG_NAME
    )
    private GATKPath testVcf;

    @ArgumentCollection
    private final SVClusterEngineArgumentsCollection clusterParameterArgs = new SVClusterEngineArgumentsCollection();

    private SAMSequenceDictionary dictionary;
    private VariantContextWriter writer;
    private SortingClusterCollapser<CrossRefSVCallRecord, CrossRefOutputCluster<CrossRefSVCallRecord>> buffer;
    private CrossRefLinkage linkage;
    private Iterator<VariantContext> testVariantSource;
    private VariantContext currentCallsetVariant;
    private Comparator<VariantContext> variantComparator;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }

        linkage = new CrossRefLinkage(dictionary);
        linkage.setDepthOnlyParams(clusterParameterArgs.getDepthParameters());
        linkage.setMixedParams(clusterParameterArgs.getMixedParameters());
        linkage.setEvidenceParams(clusterParameterArgs.getPESRParameters());

        final CrossRefSVClusterEngine engine = new CrossRefSVClusterEngine<>(linkage, dictionary);
        final CrossRefClusterCollapser collapser = new CrossRefClusterCollapser();
        buffer = new SortingClusterCollapser<>(engine, collapser, dictionary);

        variantComparator = IntervalUtils.getDictionaryOrderComparator(dictionary);
        final FeatureDataSource<VariantContext> source = new FeatureDataSource<>(testVcf.toString());
        final Object headerObject = source.getHeader();
        if ( ! (headerObject instanceof VCFHeader) ) {
            throw new GATKException("Header for " + source.getName() + " is not in VCF header format");
        }
        final VCFHeader header = (VCFHeader) headerObject;
        testVariantSource = new FeatureDataSource<VariantContext>(testVcf.toString()).iterator();

        writer = createVCFWriter(outputFile);
        writer.writeHeader(createHeader(header));
    }

    @Override
    public Object onTraversalSuccess() {
        flush(true);
        return super.onTraversalSuccess();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        if (writer != null) {
            writer.close();
        }
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext,
                      final ReferenceContext referenceContext, final FeatureContext featureContext) {
        while (currentCallsetVariant == null && testVariantSource.hasNext()) {
            currentCallsetVariant = testVariantSource.next();
            if (variantComparator.compare(currentCallsetVariant, variant) <= 0) {
                buffer.add(new CrossRefSVCallRecord(SVCallRecordUtils.create(currentCallsetVariant), true));
                currentCallsetVariant = null;
            }
        }
        buffer.add(new CrossRefSVCallRecord(SVCallRecordUtils.create(variant), false));
        flush(false);
    }

    private void flush(final boolean force) {
        buffer.flush(force).stream().map(SVCallRecordUtils::getVariantBuilder).map(VariantContextBuilder::make).forEach(writer::add);
    }

    private VCFHeader createHeader(final VCFHeader header) {
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.GENOTYPE_SUPPORT_IDS_KEY, 1, VCFHeaderLineType.Integer, "Number of cross-referenced non-ref genotypes"));
        return header;
    }

    public static final class CrossRefClusterCollapser implements SVCollapser<CrossRefSVCallRecord, CrossRefOutputCluster<CrossRefSVCallRecord>> {
        @Override
        public CrossRefSVCallRecord collapse(final CrossRefOutputCluster<CrossRefSVCallRecord> cluster) {
            final SVCallRecord testRecord = cluster.getTestItem();
            final Map<String, Long> sampleToCarrierCountsMap = cluster.getMembers().stream()
                    .map(CrossRefSVCallRecord::getCarrierGenotypeList)
                    .flatMap(Collection::stream)
                    .collect(Collectors.groupingBy(Genotype::getSampleName, Collectors.counting()));
            final GenotypesContext oldGenotypes = testRecord.getGenotypes();
            final ArrayList<Genotype> newGenotypes = new ArrayList<>(oldGenotypes.size());
            for (final Genotype g : oldGenotypes) {
                final GenotypeBuilder builder = new GenotypeBuilder(g);
                builder.attribute(GATKSVVCFConstants.GENOTYPE_SUPPORT_IDS_KEY, sampleToCarrierCountsMap.getOrDefault(g.getSampleName(), 0L));
                newGenotypes.add(builder.make());
            }
            return new CrossRefSVCallRecord(SVCallRecordUtils.copyCallWithNewGenotypes(testRecord, GenotypesContext.create(newGenotypes)), true);
        }
    }

    public static final class CrossRefLinkage extends CanonicalSVLinkage<CrossRefSVCallRecord> {

        public CrossRefLinkage(final SAMSequenceDictionary dictionary) {
            super(dictionary, false);
        }

        @Override
        public boolean areClusterable(final SVCallRecord a, final SVCallRecord b) {
            if (!algorithmsMatch(a, b)) {
                return false;
            }
            return super.areClusterable(a, b);
        }

        @Override
        protected boolean typesMatch(final SVCallRecord a, final SVCallRecord b) {
            final StructuralVariantType aType = a.getType();
            final StructuralVariantType bType = b.getType();
            if (aType == bType) {
                return true;
            }
            if (a.getType() == StructuralVariantType.BND) {
                return checkBnd(a, b);
            }
            if (b.getType() == StructuralVariantType.BND) {
                return checkBnd(b, a);
            }
            return super.typesMatch(a, b);
        }

        private boolean checkBnd(final SVCallRecord a, final SVCallRecord b) {
            if (b.isDepthOnly()) {
                return false;
            }
            return a.getStrandA() == b.getStrandA() && a.getStrandB() == b.getStrandB();
        }

        protected boolean algorithmsMatch(final SVCallRecord a, final SVCallRecord b) {
            final List<String> algorithmsA = a.getAlgorithms();
            final List<String> algorithmsB = b.getAlgorithms();
            if (algorithmsA.isEmpty() && algorithmsB.isEmpty()) {
                return true;
            }
            if (algorithmsA.size() == 1 &&  algorithmsB.size() == 1 && algorithmsA.get(0).equals(algorithmsB.get(0))) {
                return true;
            }
            final Set<String> setA = new HashSet<>(algorithmsA);
            if (setA.containsAll(algorithmsB)) {
                return true;
            }
            final Set<String> setB = new HashSet<>(algorithmsB);
            return setB.containsAll(setA);
        }
    }

    public static final class CrossRefSVCallRecord extends SVCallRecord {

        private final boolean isTestVariant;

        public CrossRefSVCallRecord(final SVCallRecord record, final boolean isTestVariant) {
            super(record);
            this.isTestVariant = isTestVariant;
        }

        public static CrossRefSVCallRecord create(final SVCallRecord record, final boolean isTestVariant) {
            return new CrossRefSVCallRecord(record, isTestVariant);
        }

        public boolean isTestVariant() {
            return isTestVariant;
        }
    }

    public static class CrossRefOutputCluster<T extends SVLocatable> extends BasicOutputCluster<T> {

        final Long testItemId;
        final T testItem;

        public CrossRefOutputCluster(final Map<Long, T> members, final Long testItemId, final T testItem) {
            super(members);
            this.testItemId = Utils.nonNull(testItemId);
            this.testItem = Utils.nonNull(testItem);
        }

        public T getTestItem() {
            return testItem;
        }

        @Override
        public Set<Long> getAllIds() {
            final Set<Long> ids = new HashSet<>(super.getAllIds());
            ids.add(testItemId);
            return ids;
        }
    }

    public static class CrossRefCluster extends BasicCluster {

        final Long testItem;

        public CrossRefCluster(final List<Long> members,
                               final String contig,
                               final int maxClusterableStart,
                               final Long testItem) {
            super(members, contig, maxClusterableStart);
            this.testItem = Utils.nonNull(testItem);
        }

        public Long getTestItem() {
            return testItem;
        }

        @Override
        public Set<Long> getAllIds() {
            final Set<Long> ids = new HashSet<>(super.getAllIds());
            ids.add(testItem);
            return ids;
        }
    }
}
