package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.*;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.BreakpointComplications;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.inference.NovelAdjacencyAndAltHaplotype;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFileUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFReader;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import javax.annotation.Nonnull;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class SvDiscoveryUtils {



    public static void evaluateIntervalsAndNarls( final SvDiscoveryInputMetaData svDiscoveryInputMetaData,
                                                  final List<NovelAdjacencyAndAltHaplotype> narls ) {
        final SAMSequenceDictionary refDict =
                svDiscoveryInputMetaData.getReferenceData().getReferenceSequenceDictionaryBroadcast().getValue();
        final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection parameters =
                svDiscoveryInputMetaData.getDiscoverStageArgs();
        final SVIntervalTree<String> narlyBreakpoints =
                readBreakpointsFromNarls(narls, refDict, parameters.truthIntervalPadding);

        final SvDiscoveryInputMetaData.SampleSpecificData sampleSpecificData =
                svDiscoveryInputMetaData.getSampleSpecificData();
        final List<SVInterval> assembledIntervals = sampleSpecificData.getAssembledIntervals();
        final List<AlignedAssemblyOrExcuse> intervalAssemblies = sampleSpecificData.intervalAssemblies;
        final SAMFileHeader header = sampleSpecificData.getHeaderBroadcast().getValue();
        if ( parameters.truthVCF == null ) {
            if ( parameters.intervalFile != null && intervalAssemblies != null ) {
                AlignedAssemblyOrExcuse.writeIntervalFile(
                        parameters.intervalFile,
                        header,
                        assembledIntervals,
                        intervalAssemblies,
                        null,
                        narlyBreakpoints );
            }
            if ( parameters.narlsFile != null ) {
                writeNarls(parameters.narlsFile, narls, null, null);
            }
        } else {
            final SVIntervalTree<String> trueBreakpoints =
                    SVVCFReader.readBreakpointsFromTruthVCF(parameters.truthVCF, refDict, parameters.truthIntervalPadding);

            final Logger toolLogger = svDiscoveryInputMetaData.getToolLogger();
            if ( assembledIntervals != null ) {
                evaluateIntervalsAgainstTruth(assembledIntervals, trueBreakpoints, toolLogger);

                if ( parameters.intervalFile != null && intervalAssemblies != null ) {
                    AlignedAssemblyOrExcuse.writeIntervalFile(
                            parameters.intervalFile,
                            header,
                            assembledIntervals,
                            intervalAssemblies,
                            trueBreakpoints,
                            narlyBreakpoints );
                }
                if ( parameters.narlsFile != null ) {
                    writeNarls(parameters.narlsFile, narls, trueBreakpoints, refDict);
                }
            }

            evaluateNarlsAgainstTruth(narlyBreakpoints, trueBreakpoints, toolLogger);
        }
    }

    private static SVIntervalTree<String> readBreakpointsFromNarls(final List<NovelAdjacencyAndAltHaplotype> narls,
                                                                   final SAMSequenceDictionary dictionary,
                                                                   final int breakpointPadding) {
        final SVIntervalTree<String> breakpoints = new SVIntervalTree<>();
        for ( final NovelAdjacencyAndAltHaplotype narl : narls ) {
            final int padding = breakpointPadding + getAmbiguity(narl.getComplication());

            final SimpleInterval si1 = narl.getLeftJustifiedLeftRefLoc();
            breakpoints.put(
                    new SVInterval(dictionary.getSequenceIndex(si1.getContig()), si1.getStart()-padding, si1.getStart()+padding), null);

            final SimpleInterval si2 = narl.getLeftJustifiedRightRefLoc();
            breakpoints.put(
                    new SVInterval(dictionary.getSequenceIndex(si2.getContig()), si2.getStart()-padding, si2.getStart()+padding), null);
        }
        return breakpoints;
    }

    private static void evaluateNarlsAgainstTruth(final SVIntervalTree<String> narlyBreakpoints,
                                                  final SVIntervalTree<String> trueBreakpoints,
                                                  final Logger localLogger) {
        final float falsePos = 1.f - narlyBreakpoints.overlapFraction(trueBreakpoints);
        final int nNarly = narlyBreakpoints.size();
        localLogger.info("Breakpoint false positive rate = " + falsePos + " (" + Math.round(falsePos*nNarly) + "/" + nNarly + ")");
        final float falseNeg = 1.f - trueBreakpoints.overlapFraction(narlyBreakpoints);
        final int nTrue = trueBreakpoints.size();
        localLogger.info("Breakpoint false negative rate = " + falseNeg + " (" + Math.round(falseNeg*nTrue) + "/" + nTrue + ")");
    }

    private static void evaluateIntervalsAgainstTruth(final List<SVInterval> assembledIntervals,
                                                      final SVIntervalTree<String> trueBreakpoints,
                                                      final Logger localLogger) {
        final SVIntervalTree<Integer> intervals = new SVIntervalTree<>();
        final int nIntervals = assembledIntervals.size();
        for ( int idx = 0; idx != nIntervals; ++idx ) {
            intervals.put(assembledIntervals.get(idx), idx);
        }
        final float falsePos = 1.f - intervals.overlapFraction(trueBreakpoints);
        localLogger.info("Interval false positive rate = " + falsePos + " (" + Math.round(falsePos*nIntervals) + "/" + nIntervals + ")");
        final float falseNeg = 1.f - trueBreakpoints.overlapFraction(intervals);
        final int nTrue = trueBreakpoints.size();
        localLogger.info("Interval false negative rate = " + falseNeg + " (" + Math.round(falseNeg*nTrue) + "/" + nTrue + ")");
    }

    // TODO: consider dups and inserts as well as micro-homology (note this was moved from BC)
    /** The uncertainty in location due to complications. */
    private static int getAmbiguity(final BreakpointComplications complications) {
        return complications.getHomologyForwardStrandRep().length();
    }

    //==================================================================================================================

    public static Set<String> getCanonicalChromosomes(final String nonCanonicalContigNamesFile,
                                                      @Nonnull final SAMSequenceDictionary dictionary) {
        final LinkedHashSet<String> allContigs = Utils.nonNull(dictionary).getSequences().stream().map(SAMSequenceRecord::getSequenceName)
                .collect(Collectors.toCollection(LinkedHashSet::new));
        if (nonCanonicalContigNamesFile == null)
            return allContigs;

        try (final Stream<String> nonCanonical = Files.lines(IOUtils.getPath(( Utils.nonNull(nonCanonicalContigNamesFile) )))) {
            nonCanonical.forEach(allContigs::remove);
            return allContigs;
        } catch ( final IOException ioe ) {
            throw new UserException("Can't read nonCanonicalContigNamesFile file "+nonCanonicalContigNamesFile, ioe);
        }
    }

    //==================================================================================================================

    /**
     * write SAM file for provided {@code filteredContigs}
     * by extracting original alignments from {@code originalAlignments},
     * to directory specified by {@code outputDir}.
     */
    public static void writeSAMRecords(final JavaRDD<GATKRead> originalAlignments, final Set<String> readNameToInclude,
                                       final String outputPath, final SAMFileHeader header) {
        final List<GATKRead> reads = originalAlignments.collect();
        writeSAMRecords(reads, readNameToInclude, outputPath, header);
    }

    public static void writeSAMRecords(final List<GATKRead> reads, final Set<String> readNameToInclude,
                                       final String outputPath, final SAMFileHeader header) {
        final SAMFileHeader cloneHeader = header.clone();
        final SAMRecordComparator localComparator;
        if (outputPath.toLowerCase().endsWith("bam")) {
            cloneHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            localComparator = new SAMRecordCoordinateComparator();
        } else if (outputPath.toLowerCase().endsWith("sam")) {
            cloneHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
            localComparator = new SAMRecordQueryNameComparator();
        } else {
            throw new IllegalArgumentException("Unsupported output format " + outputPath);
        }

        final List<SAMRecord> samRecords = new ArrayList<>();
        reads.forEach(gatkRead -> {
            if ( readNameToInclude.contains(gatkRead.getName()) ) {
                samRecords.add(gatkRead.convertToSAMRecord(cloneHeader));
            }
        });

        samRecords.sort(localComparator);
        SVFileUtils.writeSAMFile( outputPath, samRecords.iterator(), cloneHeader, true);
    }

    public static void writeNarls( final String narlsFile,
                                   final List<NovelAdjacencyAndAltHaplotype> narls,
                                   final SVIntervalTree<String> trueBreakpoints,
                                   final SAMSequenceDictionary dict ) {
        try ( final BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(narlsFile))) ) {
            for ( final NovelAdjacencyAndAltHaplotype narl : narls ) {
                writer.write(narl.toString());
                if ( trueBreakpoints != null ) {
                    if ( trueBreakpoints.hasOverlapper(toSVInterval(narl.getLeftJustifiedLeftRefLoc(), dict)) ||
                            trueBreakpoints.hasOverlapper((toSVInterval(narl.getLeftJustifiedRightRefLoc(), dict))) ) {
                        writer.write("\tT");
                    } else {
                        writer.write("\tF");
                    }
                }
                writer.newLine();
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Unable to write narls file.", ioe);
        }
    }

    private static SVInterval toSVInterval( final SimpleInterval simpleInterval, final SAMSequenceDictionary dict ) {
        return new SVInterval(dict.getSequenceIndex(simpleInterval.getContig()),
                                simpleInterval.getStart(), simpleInterval.getEnd()+1);
    }
}
