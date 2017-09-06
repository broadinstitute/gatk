package org.broadinstitute.hellbender.tools.spark.sv.evidence.experimental;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVType;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.*;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.*;
import java.util.*;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection;

@CommandLineProgramProperties(summary="Experimental code to dump statistics about regions containing small indels due to fragment length anomalies.",
        oneLineSummary="Experimental do not use.",
        usageExample="gatk-launch FindSmallIndelRegions -O hdfs://cluster-name:8020/path/to/statusDir -I hdfs://cluster-name:8020/path/to/bam --alignerIndexImage notUsed --kmersToIgnore notUsed",
        programGroup=StructuralVariationSparkProgramGroup.class)
@BetaFeature
public final class FindSmallIndelRegions extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private static final FindBreakpointEvidenceSparkArgumentCollection params =
                new FindBreakpointEvidenceSparkArgumentCollection();

    @Argument(doc = "sam file for aligned contigs", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputFile;

    @Argument(doc = "Five-column text file indicating location of true SVs.",
                fullName = "truthFileName", optional = true)
    private static String truthFile;

    @Argument(doc = "Directory into which to write CDFs.", fullName = "cdfDir", optional = true)
    private static String cdfDir;

    private static final int EVIDENCE_SIZE_GUESS = 1000;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {

        final SAMFileHeader header = getHeaderForReads();
        final JavaRDD<GATKRead> allReads = getUnfilteredReads();
        final Set<Integer> crossContigIgnoreSet = Collections.emptySet();
        final int maxTrackedFragmentLength = params.maxTrackedFragmentLength;
        final SVReadFilter filter = new SVReadFilter(params);
        final ReadMetadata readMetadata =
                new ReadMetadata(crossContigIgnoreSet, header, maxTrackedFragmentLength, allReads, filter);
        if ( params.metadataFile != null ) {
            ReadMetadata.writeMetadata(readMetadata, params.metadataFile);
        }
        final String cdfDirname = cdfDir;
        if ( cdfDir != null ) {
            for ( Map.Entry<String, LibraryStatistics> entry : readMetadata.getAllLibraryStatistics().entrySet() ) {
                final String fileName = cdfDir + "/" + entry.getKey() + ".allreads.cdf";
                dumpLibraryCDF(fileName, entry.getValue().getCDF());
            }
        }
        final Broadcast<ReadMetadata> broadcastReadMetadata = ctx.broadcast(readMetadata);
        final SVIntervalTree<SVTypeLen> trueIntervals = readTruthFile(truthFile, readMetadata.getContigNameMap());
        final Broadcast<SVIntervalTree<SVTypeLen>> broadcastTrueIntervals = ctx.broadcast(trueIntervals);
        allReads
            .mapPartitions(readItr -> {
                final List<String> statList = new ArrayList<>(EVIDENCE_SIZE_GUESS);
                final KSWindowStatusFinder KSWindowStatusFinder =
                        new KSWindowStatusFinder(broadcastReadMetadata.value(), filter,
                                                    broadcastTrueIntervals.value(), cdfDirname);
                while ( readItr.hasNext() ) {
                    KSWindowStatusFinder.testReadAndGatherStatus(readItr.next(), statList);
                }
                KSWindowStatusFinder.addKSDataStringsForCurrentWindow(statList);
                return statList.iterator();
            })
            .saveAsTextFile(outputFile);
    }

    @DefaultSerializer(SVTypeLen.Serializer.class)
    public final static class SVTypeLen {
        private SimpleSVType.TYPES type;
        private int len;

        public SVTypeLen( final String svType, final int len ) {
            this.type = SimpleSVType.TYPES.valueOf(svType);
            this.len = len;
        }

        public SVTypeLen( final Kryo kryo, final Input input ) {
            type = SimpleSVType.TYPES.values()[input.readInt()];
            len = input.readInt();
        }

        public void serialize( final Kryo kryo, final Output output ) {
            output.writeInt(type.ordinal());
            output.writeInt(len);
        }

        public SimpleSVType.TYPES getType() { return type; }
        public int getLen() { return len; }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<SVTypeLen> {
            @Override
            public void write( final Kryo kryo, final Output output, final SVTypeLen interval ) {
                interval.serialize(kryo, output);
            }

            @Override
            public SVTypeLen read( final Kryo kryo, final Input input, final Class<SVTypeLen> klass ) {
                return new SVTypeLen(kryo, input);
            }
        }
    }

    private void dumpLibraryCDF( final String fileName, final IntHistogram.CDF libraryCDF ) {
        try ( final BufferedWriter writer =
                      new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(fileName))) ) {
            writer.write(String.format("#POP\t%d\n", libraryCDF.getTotalObservations()));
            final int cdfLen = libraryCDF.size();
            for ( int idx = 0; idx != cdfLen; ++idx ) {
                writer.write(String.format("%.3f\n", libraryCDF.getFraction(idx)));
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Unable to write sample CDF " + fileName, ioe);
        }
    }

    /** Reads a 5-column, tab-delimited file showing where true SV-length indels are located:
     * ContigName  EventStart  EventEnd  EventType  EventLength
     * TODO: The truthfile should probably just be a VCF
     */
    private static SVIntervalTree<SVTypeLen> readTruthFile( final String truthFile,
                                                            final Map<String,Integer> contigNameToIDMap ) {
        final SVIntervalTree<SVTypeLen> eventMap = new SVIntervalTree<>();
        if ( truthFile == null ) return eventMap;
        try ( final BufferedReader reader = new BufferedReader(new InputStreamReader(BucketUtils.openFile(truthFile))) ) {
            String line;
            while ( (line = reader.readLine()) != null ) {
                String[] tokens = line.split("\t");
                if ( tokens.length != 5 ) {
                    throw new IOException("This line has only " + tokens.length + " columns: " + line);
                }
                final Integer contigID = contigNameToIDMap.get(tokens[0]);
                if ( contigID == null ) {
                    throw new IOException("This line has a bogus contig name: " + line);
                }
                final int start;
                try {
                    start = Integer.parseInt(tokens[1]);
                } catch ( final NumberFormatException nfe ) {
                    throw new IOException("This line has a bogus start coordinate: " + line);
                }
                final int end;
                try {
                    int tmpEnd = Integer.parseInt(tokens[2]);
                    if ( tmpEnd != start ) tmpEnd += 1;
                    end = tmpEnd;
                } catch ( final NumberFormatException nfe ) {
                    throw new IOException("This line has a bogus end coordinate: " + line);
                }
                final int len;
                try {
                    len = Integer.parseInt(tokens[4]);
                } catch ( final NumberFormatException nfe ) {
                    throw new IOException("This line has a bogus event length: " + line);
                }
                eventMap.put(new SVInterval(contigID,start,end),new SVTypeLen(tokens[3],len));
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't read truth file " + truthFile, ioe);
        }
        return eventMap;
    }

    static final class KSWindowStatusFinder extends KSWindowFinder {
        private final SVIntervalTree<SVTypeLen> trueIntervals;
        private final String cdfDir;

        public KSWindowStatusFinder( final ReadMetadata readMetadata,
                                     final SVReadFilter filter,
                                     final SVIntervalTree<SVTypeLen> trueIntervals,
                                     final String cdfDir ) {
            super(readMetadata, filter);
            this.trueIntervals = trueIntervals;
            this.cdfDir = cdfDir;
        }

        /** Not used in production.  Experimental/debug method. */
        public void testReadAndGatherStatus( final GATKRead read, final List<String> statusList ) {
            if ( !isTestable(read) ) return;

            final boolean sameContig;
            if ( !isSameBlock(read) ) {
                addKSDataStringsForCurrentWindow(statusList);
                advanceBlock(read);
            }
            addObservation(read);
        }

        /** Not used in production.  Experimental/debug method. */
        public void addKSDataStringsForCurrentWindow( final List<String> statusList ) {
            if ( curContig == null ) return;
            final int oldIdx = fillIdx ^ 1; // bit magic sets oldIdx to 1 if fillIdx is 0, and vice versa
            final int start = Math.max(1, curEnd - 2 * BLOCK_SIZE);
            for ( final Map.Entry<String, IntHistogram[]> entry : libraryToHistoPairMap.entrySet() ) {
                final IntHistogram[] histoPair = entry.getValue();
                final IntHistogram oldHisto = histoPair[oldIdx];
                oldHisto.addObservations(histoPair[fillIdx]);
                final SVInterval curInterval = new SVInterval(readMetadata.getContigID(curContig), start, curEnd);
                final Iterator<SVIntervalTree.Entry<SVTypeLen>> overlapperItr = trueIntervals.overlappers(curInterval);
                final String libName = entry.getKey();
                final long counts = oldHisto.getTotalObservations();
                if ( counts == 0L )
                    statusList.add(String.format("%s\t%s\t%s\t%d\t%d\t0",
                                            libName, overlapperItr.hasNext()?"FN":"TN", curContig, start, curEnd));
                else {
                    final LibraryStatistics stats = readMetadata.getLibraryStatistics(libName);
                    final KSData ksData = evalKS(stats.getCDF(), oldHisto.getCDF());
                    final long maxObservations = (long)(MAX_APPARENT_BASES_PER_WINDOW *stats.getReadStartFrequency());
                    final boolean isCalledInterval =
                            ksData.getSignificance() <= KS_SIGNIFICANCE &&
                            oldHisto.getTotalObservations() <= maxObservations;
                    final boolean isTrue = overlapperItr.hasNext();
                    final String truthStatus = (isCalledInterval==isTrue ? "T" : "F") +
                            (isCalledInterval ? "P" : "N");
                    final int relativeLocation = stats.getMedian() - ksData.getLocationOfMaxDiff();
                    final int firstOverlappingEventSize;
                    if ( !overlapperItr.hasNext() ) firstOverlappingEventSize = 0;
                    else {
                        final SVTypeLen overlappingTypeLen = overlapperItr.next().getValue();
                        switch ( overlappingTypeLen.getType() ) {
                            case DEL: firstOverlappingEventSize = -overlappingTypeLen.getLen(); break;
                            case INS: firstOverlappingEventSize = overlappingTypeLen.getLen(); break;
                            default: firstOverlappingEventSize = 0; break;
                        }
                    }
                    final String status = String.format("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.1e",
                                                        libName, truthStatus, curContig, start, curEnd, counts,
                                                        relativeLocation, firstOverlappingEventSize,
                                                        ksData.getMaxDiff(), ksData.getMaxArea(),
                                                        ksData.getSignificance());
                    statusList.add(status);
                    if ( cdfDir != null && (isCalledInterval || isTrue) ) {
                        final String fileName =
                                cdfDir + "/" + libName + "." + curContig + "." + start + "." + truthStatus + ".cdf";
                        dumpCDF(fileName, oldHisto.getCDF(), status, trueIntervals.overlappers(curInterval));
                    }
                }
                oldHisto.clear();
            }
        }

        private static KSData evalKS( final IntHistogram.CDF librarySample, final IntHistogram.CDF windowSample ) {
            int locationOfMaxDiff = -1;
            float maxDiff = 0.f;
            float maxArea = 0.f;
            float curArea = 0.f;
            final int size = librarySample.size();
            for (int idx = 0; idx < size; ++idx) {
                final float curDiff = windowSample.getFraction(idx) - librarySample.getFraction(idx);

                if ( curDiff*curArea >= 0.f ) {
                    curArea += curDiff;
                } else {
                    maxArea = Math.max(Math.abs(curArea), Math.abs(maxArea));
                    curArea = curDiff;
                }

                if ( Math.abs(curDiff) > Math.abs(maxDiff) ) {
                    maxDiff = curDiff;
                    locationOfMaxDiff = idx;
                }

                maxArea = Math.max(Math.abs(curArea), Math.abs(maxArea));
            }
            final long mCounts = windowSample.getTotalObservations();
            final long nCounts = librarySample.getTotalObservations();
            final double significance = 2.*Math.exp(-2.*nCounts*mCounts*maxDiff*maxDiff/(mCounts+nCounts));
            return new KSData(locationOfMaxDiff, maxDiff, maxArea, significance);
        }

        private void dumpCDF( final String fileName, final IntHistogram.CDF sampleCDF, final String summaryStatus,
                              final Iterator<SVIntervalTree.Entry<SVTypeLen>> overlappers ) {
            try ( final BufferedWriter writer =
                          new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(fileName))) ) {
                writer.write(String.format("#CDF\t%s\n", summaryStatus));
                while ( overlappers.hasNext() ) {
                    final SVIntervalTree.Entry<SVTypeLen> entry = overlappers.next();
                    final SVInterval overlappingInterval = entry.getInterval();
                    final SVTypeLen typeLen = entry.getValue();
                    final String overlapper = String.format("#OVL\t%s\t%d\t%d\t%s\t%d\n",
                            readMetadata.getContigName(overlappingInterval.getContig()),
                            overlappingInterval.getStart(), overlappingInterval.getEnd(),
                            typeLen.getType().toString(), typeLen.getLen());
                    writer.write(overlapper);
                }
                final int histoLen = sampleCDF.size();
                for ( int idx = 0; idx != histoLen; ++idx ) {
                    writer.write(String.format("%.3f\n", sampleCDF.getFraction(idx)));
                }
            } catch ( final IOException ioe ) {
                throw new UserException("Unable to write CDF " + fileName, ioe);
            }
        }
    }

    public final static class KSData {
        private int locationOfMaxDiff;
        private float maxDiff;
        private float maxArea;
        private double significance;

        public KSData( final int locationOfMaxDiff, final float maxDiff, final float maxArea,
                       final double significance ) {
            this.locationOfMaxDiff = locationOfMaxDiff;
            this.maxDiff = maxDiff;
            this.maxArea = maxArea;
            this.significance = significance;
        }

        public int getLocationOfMaxDiff() { return locationOfMaxDiff; }
        public float getMaxDiff() { return maxDiff; }
        public float getMaxArea() { return maxArea; }
        public double getSignificance() { return significance; }
    }
}
