package org.broadinstitute.hellbender.tools.spark.sv.sga;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.TextCigarCodec;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.*;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.*;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

@CommandLineProgramProperties(summary="Filter breakpoint alignments and call variants.",
        oneLineSummary="Filter breakpoint alignments and call variants",
        programGroup = StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class DiscoverVariantsFromContigAlignmentsSGASpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;
    private static final Logger localLogger = LogManager.getLogger(DiscoverVariantsFromContigAlignmentsSGASpark.class);

    @Argument(doc = "output file name for called variants", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String vcfOutput;

    @Argument(doc = "Input file of assembled contigs", shortName = "inputAssemblies",
            fullName = "inputAssemblies", optional = false)
    private String inputAssemblies;

    @Argument(doc = "Input file of contig alignments", shortName = "inputAlignments",
            fullName = "inputAlignments", optional = false)
    private String inputAlignments;

    // This class requires a reference parameter in 2bit format (to broadcast) and a reference in FASTA format
    // (to get a good sequence dictionary).
    // todo: document this better
    @Argument(doc = "FASTA formatted reference", shortName = "fastaReference",
            fullName = "fastaReference", optional = false)
    private String fastaReference;

    @Argument(doc = "To localLogger simple statistics of the contig alignments or not", shortName = "logContigAlignmentSimpleStats",
            fullName = "logContigAlignmentSimpleStats", optional = true)
    private boolean logContigAlignmentSimpleStats = false;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaRDD<AlignedContig> parsedContigAlignments
                = new SGATextFormatAlignmentParser(ctx, inputAssemblies, inputAlignments, logContigAlignmentSimpleStats ? localLogger : null).getAlignedContigs();

        DiscoverVariantsFromContigAlignmentsSAMSpark.discoverVariantsAndWriteVCF(parsedContigAlignments, fastaReference,
                ctx.broadcast(getReference()), getAuthenticatedGCSOptions(), vcfOutput, localLogger);
    }

    public static final class SGATextFormatAlignmentParser extends AlignedContigGenerator {

        private final JavaSparkContext ctx;
        private final String pathToInputAssemblies;
        private final String pathToInputAlignments;
        private final Logger toolLogger;

        SGATextFormatAlignmentParser(final JavaSparkContext ctx,
                                     final String pathToInputAssemblies,
                                     final String pathToInputAlignments,
                                     final Logger toolLogger) {
            this.ctx = ctx;
            this.pathToInputAssemblies = pathToInputAssemblies;
            this.pathToInputAlignments = pathToInputAlignments;
            this.toolLogger = toolLogger;
        }

        @SuppressWarnings("unchecked")
        @Override
        public JavaRDD<AlignedContig> getAlignedContigs() {

            final JavaPairRDD<String, byte[]> contigNameAndSequence = extractContigNameAndSequenceFromTextFile(ctx, pathToInputAssemblies);

            final JavaPairRDD<String, List<AlignedAssembly.AlignmentInterval>> contigNameAndAlignments = parseAndBreakAlignmentTextRecords(ctx, pathToInputAlignments, toolLogger);

            return contigNameAndAlignments.join(contigNameAndSequence).map(pair -> new AlignedContig(pair._1, pair._2._2, pair._2._1));
        }

        /**
         * Expected file format is such that each line represent one local assembly and
         * each line is formatted as "a numerical assembly id" + TAB + {@link ContigsCollection#toPackedFasta()}.
         * @return a pair rdd where each entry is a pair of
         *               1) contig name, formatted as {@link AlignedAssemblyOrExcuse#formatContigName(int, int)}
         *               2) contig sequence
         */
        @VisibleForTesting
        public static JavaPairRDD<String, byte[]> extractContigNameAndSequenceFromTextFile(final JavaSparkContext ctx, final String pathToInputAssemblies) {
            return ContigsCollection
                    .loadContigsCollectionKeyedByAssemblyId(ctx, pathToInputAssemblies)
                    .flatMapToPair(assemblyIdAndContigsCollection -> {
                        final ContigsCollection contigsCollection = assemblyIdAndContigsCollection._2;
                        return contigsCollection.getContents().stream()
                                .map(pair -> new Tuple2<>(AlignedAssemblyOrExcuse.formatContigName(Integer.valueOf(assemblyIdAndContigsCollection._1), Integer.valueOf(pair._1.toString().split(" ")[0].replace("contig", "").replace(">", "").replace("-", ""))), pair._2.toString().getBytes()))
                                .collect(Collectors.toList()).iterator();
                    });
        }

        @VisibleForTesting
        public static JavaPairRDD<String, List<AlignedAssembly.AlignmentInterval>> parseAndBreakAlignmentTextRecords(final JavaSparkContext ctx, final String pathToInputAlignments, final Logger toolLogger) {
            final JavaPairRDD<String, List<AlignedAssembly.AlignmentInterval>> contigNameAndAlignments
                    = ctx.textFile(pathToInputAlignments)
                    .mapToPair(SGATextFormatAlignmentParser::parseTextFileAlignmentIntervalLines)
                    .mapValues(rawAlignments -> rawAlignments.stream().flatMap(oneInterval -> Utils.stream(GappedAlignmentSplitter.split(oneInterval, SVConstants.DiscoveryStepConstants.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, oneInterval.cigarAlong5to3DirectionOfContig.getReadLength() + SVVariantDiscoveryUtils.getTotalHardClipping(oneInterval.cigarAlong5to3DirectionOfContig)))).collect(Collectors.toList()));
            if (toolLogger!=null) {
                debugStats(contigNameAndAlignments, pathToInputAlignments, toolLogger);
            }
            return contigNameAndAlignments;
        }

        /**
         * Parses fields in the same format as they were output in {@link AlignAssembledContigsSpark#formatAlignedAssemblyAsText(AlignedAssembly)}
         */
        @VisibleForTesting
        public static Tuple2<String, List<AlignedAssembly.AlignmentInterval>> parseTextFileAlignmentIntervalLines(final String textLine) {

            try {

                final String[] fields = textLine.split("\t");
                final String contigName = fields[0];
                if (fields[1].equals(AlignAssembledContigsSpark.UNMAPPED_CONTIG_STRING_REP)) {
                    return new Tuple2<>(contigName, Collections.emptyList());
                }

                final List<AlignedAssembly.AlignmentInterval> intervals = new ArrayList<>(fields.length - 1);
                for (int i = 1; i < fields.length; ++i) {
                    final String[] intervalFields = fields[i].split(AlignAssembledContigsSpark.MAPPED_CONTIG_ALIGNMENT_INTERVAL_STRING_REP_FIELD_SEPARATOR);

                    final int contigStart = Integer.valueOf(intervalFields[0].split("-")[0]);
                    final int contigEnd = Integer.valueOf(intervalFields[0].split("-")[1]);

                    intervals.add(new AlignedAssembly.AlignmentInterval(AlignAssembledContigsSpark.decodeStringAsSimpleInterval(intervalFields[1]), contigStart, contigEnd, TextCigarCodec.decode(intervalFields[2]), intervalFields[3].equals("+"), Integer.valueOf(intervalFields[4]), Integer.valueOf(intervalFields[5])));
                }
                return new Tuple2<>(contigName, intervals);
            } catch (final Exception ex) {
                throw new GATKException(textLine, ex);
            }
        }
    }

    private static void debugStats(final JavaPairRDD<String, List<AlignedAssembly.AlignmentInterval>> contigNameAndAlignments,
                                   final String outPrefix,
                                   final Logger toolLogger) {
        toolLogger.info(contigNameAndAlignments.count() + " contigs");
        final long noARs = contigNameAndAlignments.filter(tuple2 -> tuple2._2.size()==0).count();
        toolLogger.info(noARs + " contigs have no alignments");
        final long oneARs = contigNameAndAlignments.filter(tuple2 -> tuple2._2.size()==1).count();
        toolLogger.info(oneARs + " contigs have only one alignments");

        final JavaPairRDD<String, List<Tuple2<Integer, Integer>>> x = contigNameAndAlignments.filter(tuple2 -> tuple2._2.size()==2).mapToPair(tuple2 -> {
            final Iterator<AlignedAssembly.AlignmentInterval> it = tuple2._2().iterator();
            final AlignedAssembly.AlignmentInterval region1 = it.next(), region2 = it.next();
            return new Tuple2<>(tuple2._1(), Arrays.asList(new Tuple2<>(region1.mapQual, region1.referenceInterval.size()), new Tuple2<>(region2.mapQual, region2.referenceInterval.size())));
        });
        x.coalesce(1).saveAsTextFile(outPrefix+"_withTwoAlignments");
        toolLogger.info(x.count() + " contigs have two alignments");

        final JavaPairRDD<String, List<Tuple2<Integer, Integer>>> y = contigNameAndAlignments.filter(tuple2 -> tuple2._2.size()>2).mapToPair(tuple2 -> {
            final AlignedAssembly.AlignmentInterval region1 = tuple2._2().iterator().next();
            return new Tuple2<>(tuple2._1(), StreamSupport.stream(tuple2._2().spliterator(), false).map(ar -> new Tuple2<>(ar.mapQual, ar.referenceInterval.size())).collect(Collectors.toList()));
        });
        y.coalesce(1).saveAsTextFile(outPrefix+"_withMoreThanTwoAlignments");
        toolLogger.info(y.count() + " contigs have more than two alignments");
    }

}
