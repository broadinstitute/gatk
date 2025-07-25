package org.broadinstitute.hellbender.tools.walkers.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MultiVariantInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.tools.sv.SVCallRecordUtils;
import org.broadinstitute.hellbender.tools.sv.SelectSVPairs;
import org.broadinstitute.hellbender.tools.sv.cluster.CanonicalSVCollapser;
import org.broadinstitute.hellbender.tools.sv.cluster.SVClusterEngine;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;


/**
 * TODO: docs
 */
@CommandLineProgramProperties(
        summary = "Merges structural variants from two distinct callsets",
        oneLineSummary = "SV federation",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@BetaFeature
@DocumentedFeature
public class SVFederate extends SVMergingWalker {

    @ArgumentCollection
    protected SVFederationArgumentCollection inputArgs = new SVFederationArgumentCollection();

    protected VariantContextWriter writer;
    protected Map<String, Map<String, String>> sourceToPairMap;
    protected HashMap<String, String> vidAtoB;
    protected HashMap<String, String> vidBtoA;

    protected Map<String, List<Map<String, VariantContext>>> sourceToVariantMap;
    protected HashMap<String, VariantContext> vidToRecA;
    protected HashMap<String, VariantContext> vidToRecB;

    protected SAMSequenceDictionary dictionary;
    protected ReferenceSequenceFile reference;
    protected CanonicalSVCollapser collapser;


    @Override
    protected MultiVariantInputArgumentCollection getMultiVariantInputArgumentCollection() {
        return new MultiVariantInputArgumentCollection() {
            private static final long serialVersionUID = 1L;

            @Override
            public List<GATKPath> getDrivingVariantPaths() {
                // driving variants will be determined by initializeDrivingVariants()
                // directly overriding getMultiVariantInputArgumentCollection
                // to return an instance of SVFederationArgumentCollection() did not work
                return Collections.emptyList();
            }
        };
    }

    @Override
    protected void initializeDrivingVariants() {
        getDrivingVariantsFeatureInputs().addAll(inputArgs.getFeatureInputsForDrivingVariants());

        super.initializeDrivingVariants();
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();  // loads ploidy table, reference dictionary, initializes writer

        // get map of A and B feature input names to their respective SV match maps
        // order is fixed in SVFederationArgumentCollection.getDrivingVariantPaths()
        sourceToPairMap = new HashMap<>();
        sourceToPairMap.put(getDrivingVariantsFeatureInputs().get(0).getName(), vidAtoB);
        sourceToPairMap.put(getDrivingVariantsFeatureInputs().get(1).getName(), vidBtoA);

        // get map of A and B feature input names to their respective VID-to-VariantContext maps
        sourceToVariantMap = new HashMap<>();
        sourceToVariantMap.put(getDrivingVariantsFeatureInputs().get(0).getName(),
                Arrays.asList(vidToRecA, vidToRecB));
        sourceToVariantMap.put(getDrivingVariantsFeatureInputs().get(1).getName(),
                Arrays.asList(vidToRecB, vidToRecA));

        // load SV pairs
        final SelectSVPairs selector = new SelectSVPairs(inputArgs.getSVPairFilePath());
        vidAtoB = selector.getVidAToBMap();
        vidBtoA = selector.getVidBToAMap();

        reference = ReferenceUtils.createReferenceReader(referenceArguments.getReferenceSpecifier());
        // TODO breakpoint summary strategy
        collapser = new CanonicalSVCollapser(reference,
                CanonicalSVCollapser.AltAlleleSummaryStrategy.MOST_SPECIFIC_SUBTYPE,
                CanonicalSVCollapser.BreakpointSummaryStrategy.REPRESENTATIVE,
                CanonicalSVCollapser.FlagFieldLogic.OR);

    }


    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext,
                      FeatureContext featureContext) {
        final String source = variant.getSource();  // which callset is this variant from
        final Map<String, String> pairMap = sourceToPairMap.get(source);  // A-->B or B-->A matching map
        // map VID to variant for this callset (to store this variant) and the other callset (to retrieve its match)
        final Map<String, VariantContext> thisVariantMap = sourceToVariantMap.get(source).get(0);
        final Map<String, VariantContext> thatVariantMap = sourceToVariantMap.get(source).get(1);
        final String vid = variant.getID();
        if (pairMap.containsKey(vid)) {
            // get VID of matching variant in other callset
            final String match = pairMap.get(vid);
            if (thatVariantMap.containsKey(match)) {
                // retrieve the matching variant which was saved in thatVariantMap, merge them, and write
                final SVClusterEngine.OutputCluster outputCluster =
                        new SVClusterEngine.OutputCluster(Stream.of(variant, thatVariantMap.get(match))
                        .map(var -> SVCallRecordUtils.create(var, dictionary))
                        .collect(Collectors.toList()));
                // TODO: annotate cohort AFs
                // TODO: external AF annotation only mode
                final SVCallRecord merged = collapser.collapse(outputCluster);
                write(merged);  // handles conversion to VariantContext and sorting in buffer
                thatVariantMap.remove(match);  // all done with this variant - delete it to save memory
            } else {
                thisVariantMap.put(vid, variant);  // store this variant until its match is read
            }
        } else {
            // if the variant has no match, fill missing genotypes and output to sorting buffer
            write(SVCallRecordUtils.create(variant, dictionary));  // handles filling missing genotypes
        }
    }
}
