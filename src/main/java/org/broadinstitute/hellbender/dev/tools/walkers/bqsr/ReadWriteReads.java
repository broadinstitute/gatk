package org.broadinstitute.hellbender.dev.tools.walkers.bqsr;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.transforms.Create;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReadInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.ApplyBQSRTransform;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.BaseRecalOutput;
import org.broadinstitute.hellbender.dev.pipelines.bqsr.SmallBamWriter;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.ReadsSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.List;

@CommandLineProgramProperties(
        usage = "Reads a BAM file and saves it right out. This tests our I/O routines.",
        usageShort = "Reads a BAM file and saves it right out.",
        programGroup = ReadProgramGroup.class
)
public final class ReadWriteReads extends DataflowCommandLineProgram {

    @ArgumentCollection
    public final RequiredReadInputArgumentCollection readArguments = new RequiredReadInputArgumentCollection();

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    /**
     * Output path. Can be local or gs://
     */
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc="Write recalibrated reads to this new BAM file")
    public String OUTPUT;

    protected void setupPipeline(Pipeline pipeline) {
        if (readArguments.getReadFilesNames().size()>1) {
            throw new UserException("Sorry, we only support a single input file for now.");
        }
        String filename = readArguments.getReadFilesNames().get(0);
        ReadsSource readsSource = new ReadsSource(filename, pipeline);
        SAMFileHeader header = readsSource.getHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(sequenceDictionary) :
                IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        PCollection<Read> output = readsSource.getReadPCollection(intervals, ValidationStringency.STRICT);
        SmallBamWriter.writeToFile(pipeline, output, header, OUTPUT);
    }

}
