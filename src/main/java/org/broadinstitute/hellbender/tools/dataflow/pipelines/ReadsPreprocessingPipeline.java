package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.transforms.View;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OpticalDuplicatesArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.datasources.ReadsDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantsDataflowSource;
import org.broadinstitute.hellbender.engine.dataflow.transforms.composite.AddContextDataToRead;
import org.broadinstitute.hellbender.engine.datasources.ReferenceAPISource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.ApplyBQSRTransform;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalOutput;
import org.broadinstitute.hellbender.tools.dataflow.transforms.bqsr.BaseRecalibratorTransform;
import org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates.MarkDuplicates;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dataflow.SmallBamWriter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.markduplicates.OpticalDuplicateFinder;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;

import java.util.List;


@CommandLineProgramProperties(
        summary = "Takes aligned reads (likely from BWA) and runs MarkDuplicates and BQSR. The final result is analysis-ready reads.",
        oneLineSummary = "Takes aligned reads (likely from BWA) and runs MarkDuplicates and BQSR. The final result is analysis-ready reads.",
        usageExample = "Hellbender ReadsPreprocessingPipeline -I single.bam -R referenceAPIName -BQSRKnownVariants variants.vcf -O output.bam",
        programGroup = DataFlowProgramGroup.class
)

/**
 * ReadsPreprocessingPipeline is our standard pipeline that takes aligned reads (likely from BWA) and runs MarkDuplicates
 * and BQSR. The final result is analysis-ready reads.
 */
public class ReadsPreprocessingPipeline extends DataflowCommandLineProgram {
    private static final long serialVersionUID = 1L;

    @Argument(doc = "uri for the input bam, either a local file path or a gs:// bucket path",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false)
    protected String bam;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String output;

    @Argument(doc = "the reference URL, starting with " + ReferenceAPISource.URL_PREFIX, shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME,
            fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, optional = false)
    protected String referenceURL;

    @Argument(doc = "the known variants", shortName = "BQSRKnownVariants", fullName = "baseRecalibrationKnownVariants", optional = true)
    protected List<String> baseRecalibrationKnownVariants;

    @ArgumentCollection
    protected OpticalDuplicatesArgumentCollection opticalDuplicatesArgumentCollection = new OpticalDuplicatesArgumentCollection();

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Override
    protected void setupPipeline( Pipeline pipeline ) {
        // Load the reads.
        final ReadsDataflowSource readsDataflowSource = new ReadsDataflowSource(bam, pipeline);
        final SAMFileHeader readsHeader = readsDataflowSource.getHeader();
        final List<SimpleInterval> intervals = intervalArgumentCollection.getSpecifiedOrAllIntervals(readsHeader.getSequenceDictionary());

        final PCollectionView<SAMFileHeader> headerSingleton = ReadsDataflowSource.getHeaderView(pipeline, readsHeader);
        final PCollection<GATKRead> initialReads = readsDataflowSource.getReadPCollection(intervals);

        final OpticalDuplicateFinder finder = opticalDuplicatesArgumentCollection.READ_NAME_REGEX != null ?
                new OpticalDuplicateFinder(opticalDuplicatesArgumentCollection.READ_NAME_REGEX, opticalDuplicatesArgumentCollection.OPTICAL_DUPLICATE_PIXEL_DISTANCE, null) : null;
        final PCollectionView<OpticalDuplicateFinder> finderPcolView = pipeline.apply(Create.of(finder)).apply(View.<OpticalDuplicateFinder>asSingleton());

        // Apply MarkDuplicates to produce updated GATKReads.
        final PCollection<GATKRead> markedReads = initialReads.apply(new MarkDuplicates(headerSingleton, finderPcolView));

        // Load the Variants and the Reference and join them to reads.
        final VariantsDataflowSource variantsDataflowSource = new VariantsDataflowSource(baseRecalibrationKnownVariants, pipeline);

        // Use the BQSR_REFERENCE_WINDOW_FUNCTION so that the reference bases required by BQSR for each read are fetched
        final ReferenceMultiSource referenceDataflowSource = new ReferenceMultiSource(pipeline.getOptions(), referenceURL, BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION);

        final PCollection<KV<GATKRead, ReadContextData>> readsWithContext = AddContextDataToRead.add(markedReads, referenceDataflowSource, variantsDataflowSource);

        // BQSR.
        // default arguments are best practice.
        RecalibrationArgumentCollection recalArgs = new RecalibrationArgumentCollection();
        final SAMSequenceDictionary readsDictionary = readsHeader.getSequenceDictionary();
        final SAMSequenceDictionary refDictionary = referenceDataflowSource.getReferenceSequenceDictionary(readsDictionary);
        checkSequenceDictionaries(refDictionary, readsDictionary);
        PCollectionView<SAMSequenceDictionary> refDictionaryView = pipeline.apply(Create.of(refDictionary)).setName("refDictionary").apply(View.asSingleton());
        BaseRecalibratorTransform baseRecalibrator = new BaseRecalibratorTransform(headerSingleton, refDictionaryView, recalArgs);
        final PCollection<BaseRecalOutput> recalibrationReports = readsWithContext.apply(baseRecalibrator).apply(baseRecalibrator.toBaseRecalOutput());
        final PCollectionView<BaseRecalOutput> mergedRecalibrationReport = recalibrationReports.apply(View.<BaseRecalOutput>asSingleton());

        final ApplyBQSRArgumentCollection applyArgs = new ApplyBQSRArgumentCollection();
        final PCollection<GATKRead> finalReads = markedReads.apply(new ApplyBQSRTransform(headerSingleton, mergedRecalibrationReport, applyArgs));
        SmallBamWriter.writeToFile(pipeline, finalReads, readsHeader, output);
    }

    private void checkSequenceDictionaries(final SAMSequenceDictionary refDictionary, SAMSequenceDictionary readsDictionary) {
        Utils.nonNull(refDictionary);
        Utils.nonNull(readsDictionary);
        SequenceDictionaryUtils.validateDictionaries("reference", refDictionary, "reads", readsDictionary);
    }

}
