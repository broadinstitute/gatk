package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
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
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVLinkage;
import org.broadinstitute.hellbender.tools.sv.cluster.CrossRefSVClusterEngine;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngineArgumentsCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;

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
    private CrossRefLinkage linkage;
    private Iterator<VariantContext> testVariantSource;
    private VariantContext currentTestVariant;
    private Comparator<VariantContext> variantComparator;
    private CrossRefSVClusterEngine engine;
    private Long nextItemId = 0L;
    private String currentContig = null;

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
        engine = new CrossRefSVClusterEngine(linkage, this::collapse, dictionary);

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
        flushTestVariants(null);
        flushClusters(true);
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
        flushTestVariants(variant);
        if (currentContig == null) {
            currentContig = variant.getContig();
        } else if (!currentContig.equals(variant.getContig())) {
            currentContig = variant.getContig();
            flushClusters(true);
        }
        engine.add(nextItemId++, SVCallRecordUtils.create(variant), true);
        flushClusters(false);
    }

    private void flushTestVariants(final VariantContext variant) {
        while (true) {
            if (currentTestVariant == null && testVariantSource.hasNext()) {
                currentTestVariant = testVariantSource.next();
            }
            if (currentTestVariant != null) {
                if (currentContig == null) {
                    currentContig = currentTestVariant.getContig();
                } else if (!currentTestVariant.getContig().equals(currentContig)) {
                    flushClusters(true);
                    currentContig = currentTestVariant.getContig();
                }
            }
            if (currentTestVariant != null && (variant == null || variantComparator.compare(currentTestVariant, variant) <= 0)) {
                engine.add(nextItemId++, SVCallRecordUtils.create(currentTestVariant), false);
                currentTestVariant = null;
            } else {
                break;
            }
        }
    }

    private void flushClusters(final boolean force) {
        engine.flush(force).stream()
                .map(SVCallRecordUtils::getVariantBuilder)
                .map(VariantContextBuilder::make)
                .forEach(writer::add);
    }

    private VCFHeader createHeader(final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CARRIER_COUNT_INFO, 1, VCFHeaderLineType.Integer, "Number of non-ref samples"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CARRIER_FREQ_INFO, 1, VCFHeaderLineType.Float, "Fraction of non-ref samples"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CROSS_REFERENCE_SITE_SUPPORT_INFO, 1, VCFHeaderLineType.Integer, "Number of cross-referenced sites matching this variant"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CROSS_REFERENCE_CARRIER_COUNT_INFO, 1, VCFHeaderLineType.Integer, "Number of cross-referenced non-ref samples"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CROSS_REFERENCE_CARRIER_FREQ_INFO, 1, VCFHeaderLineType.Float, "Fraction of cross-referenced non-ref samples"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GATKSVVCFConstants.CROSS_REFERENCE_MEMBERS_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Cross-referenced variant IDs"));
        header.addMetaDataLine(new VCFFormatHeaderLine(GATKSVVCFConstants.CROSS_REFERENCE_GENOTYPE_SUPPORT_FORMAT, 1, VCFHeaderLineType.Integer, "Number of cross-referenced non-ref genotypes"));
        return header;
    }

    public SVCallRecord collapse(final CrossRefSVClusterEngine.CrossRefOutputCluster cluster) {
        final SVCallRecord testRecord = cluster.getTestItem();
        final Map<String, Long> sampleToCarrierCountsMap = cluster.getItems().stream()
                .map(SVCallRecord::getCarrierGenotypeList)
                .flatMap(Collection::stream)
                .collect(Collectors.groupingBy(Genotype::getSampleName, Collectors.counting()));
        final int xrCarrierCount = sampleToCarrierCountsMap.size();
        final int sampleCount = testRecord.getGenotypes().size();
        final int carrierCount = cluster.getTestItem().getCarrierGenotypeList().size();
        final GenotypesContext oldGenotypes = testRecord.getGenotypes();
        final ArrayList<Genotype> newGenotypes = new ArrayList<>(oldGenotypes.size());
        for (final Genotype g : oldGenotypes) {
            final GenotypeBuilder builder = new GenotypeBuilder(g);
            final Long supportCount = sampleToCarrierCountsMap.getOrDefault(g.getSampleName(), 0L);
            builder.attribute(GATKSVVCFConstants.CROSS_REFERENCE_GENOTYPE_SUPPORT_FORMAT, supportCount);
            newGenotypes.add(builder.make());
        }
        final SVCallRecord recordWithGenotypes = SVCallRecordUtils.copyCallWithNewGenotypes(testRecord, GenotypesContext.create(newGenotypes));
        final Map<String, Object> attributes = new HashMap<>(recordWithGenotypes.getAttributes());
        final List<String> members = cluster.getItems().stream().map(SVCallRecord::getId).collect(Collectors.toList());
        attributes.put(GATKSVVCFConstants.CROSS_REFERENCE_SITE_SUPPORT_INFO, cluster.getItems().size());
        attributes.put(GATKSVVCFConstants.CROSS_REFERENCE_CARRIER_COUNT_INFO, xrCarrierCount);
        attributes.put(GATKSVVCFConstants.CROSS_REFERENCE_CARRIER_FREQ_INFO, xrCarrierCount / (double) sampleCount);
        attributes.put(GATKSVVCFConstants.CARRIER_COUNT_INFO, carrierCount);
        attributes.put(GATKSVVCFConstants.CARRIER_FREQ_INFO, carrierCount / (double) sampleCount);
        attributes.put(GATKSVVCFConstants.CROSS_REFERENCE_MEMBERS_INFO, String.join(",", members));
        return SVCallRecordUtils.copyCallWithNewAttributes(recordWithGenotypes, attributes);
    }

    public static final class CrossRefLinkage extends CanonicalSVLinkage<SVCallRecord> {

        public CrossRefLinkage(final SAMSequenceDictionary dictionary) {
            super(dictionary, false);
        }

        @Override
        public boolean typesMatch(final SVCallRecord a, final SVCallRecord b) {
            final StructuralVariantType aType = a.getType();
            final StructuralVariantType bType = b.getType();
            if (aType == bType) {
                return true;
            }
            if (aType == StructuralVariantType.BND || bType == StructuralVariantType.BND) {
                return strandsMatch(a, b) && a.isIntrachromosomal() == b.isIntrachromosomal();
            }
            if (a.isSimpleCNV() && b.isSimpleCNV()) {
                if (aType == StructuralVariantType.CNV || bType == StructuralVariantType.CNV) {
                    return true;
                }
            }
            return false;
        }
    }

}
