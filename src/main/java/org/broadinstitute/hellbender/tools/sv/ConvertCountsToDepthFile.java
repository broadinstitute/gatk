package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.codecs.*;
import org.broadinstitute.hellbender.utils.io.FeatureOutputStream;


import java.util.Collections;

@CommandLineProgramProperties(
        summary = "...",
        oneLineSummary = "...",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
public class ConvertCountsToDepthFile extends FeatureWalker<SimpleCount> {
    public static final String COMPRESSION_LEVEL_NAME = "compression-level";

    @Argument(
            doc = "...",
            fullName = "counts-filename", //StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.FEATURE_SHORT_NAME
    )
    private FeatureInput<SimpleCount> countsPath;

    @Argument(
            doc = "...",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private GATKPath outputFilePath;

    @Argument(
            doc = "...",
            fullName = "sample-name",
            shortName = "x"
    )
    private String sampleName;

    @Argument(
            doc = "Output compression level",
            fullName = COMPRESSION_LEVEL_NAME,
            minValue = 0, maxValue = 9, optional = true
    )
    private int compressionLevel = 4;

    private SAMSequenceDictionary dictionary;

    private FeatureOutputStream<DepthEvidence> outputSink;

    @Override
    protected boolean isAcceptableFeatureType(Class<? extends Feature> featureType) {
        return SimpleCount.class.isAssignableFrom(featureType);
    }

    @Override
    public void onTraversalStart() {
        final SampleLocatableMetadata header = (SampleLocatableMetadata)getDrivingFeaturesHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();

        final FeatureOutputCodec<? extends Feature, ? extends FeatureSink<? extends Feature>> codec =
                FeatureOutputCodecFinder.find(outputFilePath);

        if (!DepthEvidence.class.isAssignableFrom(codec.getFeatureType()))
            throw new UserException.BadInput(
                    outputFilePath.getURIString() + " is not a valid format for " +
                    DepthEvidence.class.getSimpleName());


        @SuppressWarnings("unchecked")
        final FeatureOutputCodec<DepthEvidence, FeatureOutputStream<DepthEvidence>> depthEvidenceCodec =
                (FeatureOutputCodec<DepthEvidence, FeatureOutputStream<DepthEvidence>>) codec;

        outputSink = depthEvidenceCodec.makeSink(
                outputFilePath,
                sequenceDictionary,
                Collections.singletonList(sampleName),
                compressionLevel
        );


        //Object sink = codec.makeSink(outputFilePath, sequenceDictionary, new ArrayList<>(Collections.singletonList(sampleName)), compressionLevel);
        //FeatureOutputStream<DepthEvidence> x = codec.makeSink(outputFilePath, sequenceDictionary, new ArrayList<>(Collections.singletonList(sampleName)), compressionLevel);

        /*
        outputSink = (FeatureSink<SimpleCount>)codec.makeSortMerger(outputFilePath,
                getDictionary(), new ArrayList<>(sampleNames), compressionLevel);*/

        /*
        final String headerString = (String) depthEvidenceCodec.getHeader(
                Collections.singletonList(sampleName),
                sequenceDictionary
        );*/

        //outputSink.writeHeader(headerString);

        //outputSink.writeHeader();
    }

    @Override
    public void apply(SimpleCount feature, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        final DepthEvidence depthEvidence = new DepthEvidence(
                feature.getContig(),
                feature.getStart(),
                feature.getEnd(),
                new int[] {feature.getCount()});

        outputSink.write(depthEvidence);
    }

    @Override
    public GATKPath getDrivingFeaturePath() {
        // TODO: not sure if this is intended to be implemented as the following or what purpose it serves
        return new GATKPath(countsPath.getFeaturePath());
    }

    @Override
    public Object onTraversalSuccess() {
        outputSink.close();
        return super.onTraversalSuccess();
    }
}
