package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.transforms.*;
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
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.dev.DoFnWLog;
import org.broadinstitute.hellbender.engine.dataflow.*;
import org.broadinstitute.hellbender.engine.dataflow.datasources.*;
import org.broadinstitute.hellbender.engine.dataflow.transforms.composite.AddContextDataToRead;
import org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates.MarkDuplicates;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.tools.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.SmallBamWriter;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.sql.Ref;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;


@CommandLineProgramProperties(
        summary = "Takes aligned reads (likely from BWA) and runs MarkDuplicates and BQSR. The final result is analysis-ready reads.",
        oneLineSummary = "Takes aligned reads (likely from BWA) and runs MarkDuplicates and BQSR. The final result is analysis-ready reads.",
        usageExample = "Hellbender ReadsPreprocessingPipeline -I single.bam -R referenceName -BQSRKnownVariants variants.vcf -O output.bam",
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

    @Argument(doc = "the reference name", shortName = StandardArgumentDefinitions.REFERENCE_SHORT_NAME,
            fullName = StandardArgumentDefinitions.REFERENCE_LONG_NAME, optional = false)
    protected String referenceName;

    @Argument(doc = "the known variants", shortName = "BQSRKnownVariants", fullName = "baseRecalibrationKnownVariants", optional = true)
    protected List<String> baseRecalibrationKnownVariants;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Override
    protected void setupPipeline( Pipeline pipeline ) {
        // Load the reads.
        final ReadsDataflowSource readsDataflowSource = new ReadsDataflowSource(bam, pipeline);
        final SAMFileHeader readsHeader = readsDataflowSource.getHeader();
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(readsHeader.getSequenceDictionary())
                : IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());

        final PCollectionView<SAMFileHeader> headerSingleton = ReadsDataflowSource.getHeaderView(pipeline, readsHeader);
        final PCollection<GATKRead> initialReads = readsDataflowSource.getReadPCollection(intervals);

        // Apply MarkDuplicates to produce updated GATKReads.
        final PCollection<GATKRead> markedReads = initialReads.apply(new MarkDuplicates(headerSingleton));

        // Load the Variants and the Reference and join them to reads.
        final VariantsDataflowSource variantsDataflowSource = new VariantsDataflowSource(baseRecalibrationKnownVariants, pipeline);

        Map<String, String> referenceNameToIdTable = RefAPISource.buildReferenceNameToIdTable(pipeline.getOptions(), referenceName);
        RefAPIMetadata refAPIMetadata = new RefAPIMetadata(referenceName, referenceNameToIdTable);

        final PCollection<KV<GATKRead, ReadContextData>> readsWithContext = AddContextDataToRead.add(markedReads, refAPIMetadata, variantsDataflowSource);

        // Apply BQSR.
        final PCollection<RecalibrationTables> recalibrationReports = readsWithContext.apply(new BaseRecalibratorStub(headerSingleton));
        final PCollectionView<RecalibrationTables> mergedRecalibrationReport = recalibrationReports.apply(View.<RecalibrationTables>asSingleton());

        final PCollection<GATKRead> finalReads = markedReads.apply(new ApplyBQSRStub(headerSingleton, mergedRecalibrationReport));
        SmallBamWriter.writeToFile(pipeline, finalReads, readsHeader, output);
    }

    private static class BaseRecalibratorStub extends PTransform<PCollection<KV<GATKRead, ReadContextData>>, PCollection<RecalibrationTables>> {
        private static final long serialVersionUID = 1L;
        private final PCollectionView<SAMFileHeader> header;

        public BaseRecalibratorStub( final PCollectionView<SAMFileHeader> header ) {
            this.header = header;
        }

        @Override
        public PCollection<RecalibrationTables> apply( PCollection<KV<GATKRead, ReadContextData>> input ) {
            return input.apply(ParDo.named("BaseRecalibrator").
                    of(new DoFnWLog<KV<GATKRead, ReadContextData>, RecalibrationTables>("BaseRecalibratorStub") {
                        private static final long serialVersionUID = 1L;

                        @Override
                        public void processElement( ProcessContext c ) throws Exception {
                            c.output(new RecalibrationTables(new StandardCovariateList(new RecalibrationArgumentCollection(), Collections.emptyList())));
                        }
                    }).withSideInputs(header));
        }
    }

    private static class ApplyBQSRStub extends PTransform<PCollection<GATKRead>, PCollection<GATKRead>> {
        private static final long serialVersionUID = 1L;
        private final PCollectionView<SAMFileHeader> header;
        private final PCollectionView<RecalibrationTables> recalibrationReport;

        public ApplyBQSRStub( final PCollectionView<SAMFileHeader> header, final PCollectionView<RecalibrationTables> recalibrationReport ) {
            this.header = header;
            this.recalibrationReport = recalibrationReport;
        }

        @Override
        public PCollection<GATKRead> apply( PCollection<GATKRead> input ) {
            return input.apply(ParDo.named("ApplyBQSR").
                    of(new DoFnWLog<GATKRead, GATKRead>("ApplyBQSRStub") {
                        private static final long serialVersionUID = 1L;

                        @Override
                        public void processElement( ProcessContext c ) throws Exception {
                            c.output(c.element());
                        }
                    }).withSideInputs(header, recalibrationReport));
        }
    }

}
