package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.collect.Sets;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedInterval;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalHeader;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalUtils;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalWriter;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.SimpleAnnotatedIntervalWriter;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        oneLineSummary = "Merge annotated genomic regions within specified distance if annotation value(s) are exactly the same.",
        summary = "This simple tool will merge genomic regions (specified in a tsv) when given annotations (columns) contain exact values in neighboring segments and the segments are within a specified maximum genomic distance.",
        programGroup = CopyNumberProgramGroup.class)

@ExperimentalFeature
public class MergeAnnotatedRegionsByAnnotation extends GATKTool {

    final static public String ANNOTATIONS_TO_MATCH = "annotations-to-match";
    final static public String MAX_MERGE_DISTANCE = "max-merge-distance";

    @Argument(
            doc = "Input segment file or annotated segment file.",
            fullName = CopyNumberStandardArgument.SEGMENTS_FILE_LONG_NAME
    )
    private File segmentFile;

    @Argument(
            doc = "Annotation(s) to merge on.",
            fullName = ANNOTATIONS_TO_MATCH
    )
    private List<String> annotationsToMatch;

    @Argument(
            doc = "Maximum distance to merge in bases",
            fullName = MAX_MERGE_DISTANCE,
            optional = true
    )
    private Integer maxMergeDistance = 1000000;

    @Argument(
            doc = "Output TSV file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private File outputFile;

    /** Separator to use when there are conflicts of values not in the input annotation list. */
    final static String DEFAULT_SEPARATOR = "__";

    @Argument(
            doc = "Output contig column name.",
            fullName = "output-contig-column",
            optional = true
    )
    private String outputContigColumn = "CONTIG";

    @Argument(
            doc = "Output start position column name.",
            fullName = "output-start-column",
            optional = true
    )
    private String outputStartColumn = "START";

    @Argument(
            doc = "Output end column name.",
            fullName = "output-end-column",
            optional = true
    )
    private String outputEndColumn = "END";


    @Override
    public void traverse() {

        // Load the seg file.  Not done with the collection class, in order to keep the sequence dictionary optional on the input.
        final AnnotatedIntervalCollection annotatedIntervalCollection = AnnotatedIntervalCollection.create(segmentFile.toPath(), null);
        final List<AnnotatedInterval> initialSegments = annotatedIntervalCollection.getRecords();

        Utils.validateArg(annotatedIntervalCollection.getAnnotations().containsAll(annotationsToMatch), "Input file did not have all of the specified annotations.  Missing annotations were: " +
                Sets.difference(new HashSet<>(annotationsToMatch), new HashSet<>(annotatedIntervalCollection.getAnnotations())).stream().collect(Collectors.joining(", "))
               );

        Utils.validateArg(getBestAvailableSequenceDictionary() != null, "Sequence dictionary not available in the input file nor specified in a reference parameter.  Please specify a reference with the -R parameter for this input file.");

        // Perform the actual merging
        final List<AnnotatedInterval> mergedSegments = AnnotatedIntervalUtils.mergeRegionsByAnnotation(initialSegments, getBestAvailableSequenceDictionary(), annotationsToMatch,
                l -> progressMeter.update(l), DEFAULT_SEPARATOR, maxMergeDistance);

        final AnnotatedIntervalHeader header = new AnnotatedIntervalHeader(outputContigColumn, outputStartColumn, outputEndColumn, new ArrayList<>(mergedSegments.get(0).getAnnotations().keySet()), null);
        final AnnotatedIntervalWriter writer = new SimpleAnnotatedIntervalWriter(outputFile);
        writer.writeHeader(header);
        mergedSegments.forEach(s -> writer.add(s));
        writer.close();
    }
}
