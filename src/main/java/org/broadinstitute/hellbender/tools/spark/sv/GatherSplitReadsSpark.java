package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.*;
import org.apache.hadoop.io.NullWritable;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.seqdoop.hadoop_bam.BAMRecordWriter;
import org.seqdoop.hadoop_bam.SAMRecordWritable;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;

@CommandLineProgramProperties(summary="Gather clustered split reads using spark",
        oneLineSummary="Gather clustered split reads using spark",
        programGroup = SparkProgramGroup.class)
public class GatherSplitReadsSpark extends GATKSparkTool
{
    private static final long serialVersionUID = 1l;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    // A shim to allow writing a stream of SAMRecords as a BAM.
    private static final class BAMWriter extends BAMRecordWriter<NullWritable> {
        public BAMWriter( OutputStream output, SAMFileHeader header ) throws IOException
        { super(output, header, true); }

        public void write( SAMRecord rec ) { writeAlignment(rec); }
        public void close() throws IOException { close(null); }

        @Override public void write( NullWritable ignored, SAMRecordWritable rec )
        { throw new UnsupportedOperationException(); }
    }

    // A shim to transform an iterator over GATKReads into an iterator over SAMRecords
    private static final class GATKRead2SAMRecordIterator implements Iterator<SAMRecord>
    {
        GATKRead2SAMRecordIterator( Iterator<GATKRead> someIterator ) { iterator = someIterator; }

        @Override public boolean hasNext() { return iterator.hasNext(); }
        @Override public SAMRecord next() { return iterator.next().convertToSAMRecord(null); }

        private final Iterator<GATKRead> iterator;
    }

    private static final class EventEvidence implements Comparable<EventEvidence>
    {
        public EventEvidence( final int someContigIdx, final int someLocus )
        {
            contigIdx = someContigIdx;
            locus = someLocus;
            counter = 0;
        }

        public EventEvidence( final int someContigIdx, final int someLocus, final long someCounter )
        {
            contigIdx = someContigIdx;
            locus = someLocus;
            counter = someCounter;
        }

        public int getContigIdx() { return contigIdx; }
        public int getLocus() { return locus; }

        public int compareTo( final EventEvidence ee )
        {
            int result = Integer.compare(contigIdx,ee.contigIdx);
            if ( result == 0 ) result = Integer.compare(locus,ee.locus);
            if ( result == 0 ) result = Long.compare(counter,ee.counter);
            return result;
        }

        private final int contigIdx;
        private final int locus;
        private final long counter; // used as a uniquifier
    }

    private static final class SplitReadClusterer implements Iterator<SAMRecord>
    {
        public SplitReadClusterer( final Iterator<SAMRecord> someInputIterator )
        {
            inputIterator = someInputIterator;
        }

        @Override
        public boolean hasNext()
        {
            while ( !outputIterator.hasNext() )
            {
                if ( !inputIterator.hasNext() )
                    return false;
                outputIterator = processRead(inputIterator.next());
            }
            return true;
        }

        @Override
        public SAMRecord next()
        {
            if ( !hasNext() )
                throw new NoSuchElementException("Iterator<GATKRead> is exhausted.");
            return outputIterator.next();
        }

        private Iterator<SAMRecord> processRead( final SAMRecord read )
        {
            if ( !read.getReadFailsVendorQualityCheckFlag() && !read.getReadUnmappedFlag() &&
                    !read.getDuplicateReadFlag() && read.getMappingQuality() > 0 &&
                    read.getReferenceIndex() != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX )
            {
                final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
                int matchLen = 0;
                for ( CigarElement ele : cigarElements )
                    if (ele.getOperator() == CigarOperator.MATCH_OR_MISMATCH)
                        matchLen += ele.getLength();
                if ( matchLen >= MIN_MATCH_LEN )
                {
                    final ListIterator<CigarElement> itr0 = cigarElements.listIterator();
                    if ( itr0.hasNext() )
                    {
                        CigarElement firstEle = itr0.next();
                        if ( firstEle.getOperator() == CigarOperator.HARD_CLIP && itr0.hasNext() )
                            firstEle = itr0.next();
                        if ( firstEle.getOperator() == CigarOperator.SOFT_CLIP &&
                                firstEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                                highQualityRegion(read.getBaseQualities(), 0) )
                            return cluster(read, read.getStart()); // NON-STRUCTURED return
                    }
                    final ListIterator<CigarElement> itrN = cigarElements.listIterator(cigarElements.size());
                    if ( itrN.hasPrevious() )
                    {
                        CigarElement lastEle = itrN.previous();
                        if ( lastEle.getOperator() == CigarOperator.HARD_CLIP && itrN.hasPrevious() )
                            lastEle = itrN.previous();
                        if ( lastEle.getOperator() == CigarOperator.SOFT_CLIP &&
                                lastEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                                highQualityRegion(read.getBaseQualities(), read.getReadLength()-lastEle.getLength()) )
                            return cluster(read, read.getEnd()); // NON-STRUCTURED return
                    }
                }
            }
            return Collections.emptyIterator();
        }

        private Iterator<SAMRecord> cluster( final SAMRecord read, final int locus )
        {
            final EventEvidence newEvidence = new EventEvidence(read.getReferenceIndex(),locus,++counter);
            locMap.put(newEvidence,read);

            // clean out old stuff that can't possibly be interesting anymore
            final Iterator<Map.Entry<EventEvidence,SAMRecord>> itr = locMap.entrySet().iterator();
            while ( itr.hasNext() )
            {
                final EventEvidence oldEvidence = itr.next().getKey();
                if ( oldEvidence.getContigIdx() == newEvidence.getContigIdx() &&
                        oldEvidence.getLocus() + MAX_LOCUS_DIST > newEvidence.getLocus() )
                    break;
                itr.remove();
            }

            // find all the evidence in a window surrounding the locus of the new evidence
            final SortedMap<EventEvidence,SAMRecord> windowMap = locMap
                    .tailMap(new EventEvidence(newEvidence.getContigIdx(), newEvidence.getLocus()-CLUSTER_WINDOW))
                    .headMap(new EventEvidence(newEvidence.getContigIdx(), newEvidence.getLocus()+CLUSTER_WINDOW + 1));
            if ( windowMap.size() >= MIN_EVIDENCE )
            {
                final List<SAMRecord> list = new LinkedList<>();
                for ( Map.Entry<EventEvidence,SAMRecord> entry : windowMap.entrySet() )
                {
                    if (entry.getValue() != null)
                    {
                        list.add(entry.getValue());
                        entry.setValue(null);
                    }
                }
                return list.iterator(); // NON-STRUCTURED return
            }
            return Collections.emptyIterator();
        }

        private static boolean highQualityRegion( final byte[] quals, int idx )
        {
            for ( final int end = idx+MIN_SOFT_CLIP_LEN; idx != end; ++idx )
                if ( quals[idx] < MIN_QUALITY )
                    return false; // NON-STRUCTURED return
            return true;
        }

        private long counter = 0;
        private final Iterator<SAMRecord> inputIterator;
        private final SortedMap<EventEvidence,SAMRecord> locMap = new TreeMap<>();
        private Iterator<SAMRecord> outputIterator = Collections.emptyIterator();

        private static final int MIN_MATCH_LEN = 45; // minimum length of matched portion of an interesting alignment
        private static final int MIN_SOFT_CLIP_LEN = 30; // minimum length of an interesting soft clip
        private static final byte MIN_QUALITY = 15; // minimum acceptable quality in a soft-clip window
        private static final int MAX_LOCUS_DIST = 500; // stale evidence distance, should be somewhat longer than a read
        private static final int CLUSTER_WINDOW = 2; // size of locus window in which we cluster event evidence
        private static final int MIN_EVIDENCE = 3; // minimum evidence count in a cluster
    }

    private static final class SAMWindowSorter implements Iterator<SAMRecord>, Iterable<SAMRecord>
    {
        public SAMWindowSorter( final Iterator<SAMRecord> someInputIterator, final int someWindowSize )
        {
            inputIterator = someInputIterator;
            windowSize = someWindowSize;
            fillBuffer();
        }

        @Override
        public boolean hasNext()
        {
            return !recordSet.isEmpty();
        }

        @Override
        public SAMRecord next()
        {
            final Iterator<SAMRecord> itr = recordSet.iterator();
            final SAMRecord result = itr.next();
            itr.remove();
            fillBuffer();
            return result;
        }

        @Override
        public Iterator<SAMRecord> iterator() { return this; }

        private void fillBuffer()
        {
            while ( getSetSpan() <= windowSize && inputIterator.hasNext() )
                recordSet.add(inputIterator.next());
        }

        private int getSetSpan()
        {
            if ( recordSet.isEmpty() )
                return 0;
            final SAMRecord first = recordSet.first();
            final SAMRecord last = recordSet.last();
            if ( !first.getReferenceIndex().equals(last.getReferenceIndex()) )
                return Integer.MAX_VALUE;
            return last.getAlignmentStart() - first.getAlignmentStart();
        }

        private static final class SAMComparator implements Comparator<SAMRecord>
        {
            @Override
            public int compare( final SAMRecord rec1, final SAMRecord rec2 )
            {
                int result = rec1.getReferenceIndex().compareTo(rec2.getReferenceIndex());
                if ( result == 0 ) result = Integer.compare(rec1.getAlignmentStart(), rec2.getAlignmentStart());
                if ( result == 0 ) result = rec1.getReadName().compareTo(rec2.getReadName());
                if ( result == 0 ) result = Integer.compare(rec1.getFlags(),rec2.getFlags());
                return result;
            }
        }

        private final SortedSet<SAMRecord> recordSet = new TreeSet<SAMRecord>(new SAMComparator());
        private final Iterator<SAMRecord> inputIterator;
        private final int windowSize;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx )
    {
        final Broadcast<SAMFileHeader> header = ctx.broadcast(getHeaderForReads());
        if ( header.value().getSortOrder() != SAMFileHeader.SortOrder.coordinate )
            throw new IllegalStateException("The BAM must be coordinate sorted.");

        getReads()
            .mapPartitions(readItr ->
            {
                final SplitReadClusterer itr = new SplitReadClusterer(new GATKRead2SAMRecordIterator(readItr));
                return new SAMWindowSorter(itr,SAM_WINDOW_SIZE);
            }, true)
            .coalesce(1)
            .foreachPartition(readItr ->
            {
                final BAMWriter writer = new BAMWriter(new FileOutputStream(output), header.value());
                while (readItr.hasNext()) writer.write(readItr.next());
                writer.close();
            });
    }

    private static final int SAM_WINDOW_SIZE = 500;
}
