package org.broadinstitute.hellbender.tools.walkers.variantrecalling;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;

import java.nio.channels.SeekableByteChannel;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;

/**
 * a service class for HaplotypeBasedVariableRecaller that reads SAM/BAM files.
 *
 * For each given location (query location) it returns reads from all files that fall into the region after
 * trimming them accordign to span of a given variant context
 */
public class TrimmedReadsReader {

    private final List<SamReader>         samReaders = new LinkedList<>();
    private CountingReadFilter            readFilter;
    private final Map<String, Integer>    readGroupMaxClass = new LinkedHashMap<>();
    private final Map<String, String>     readGroupFlowOrder = new LinkedHashMap<>();
    private final FlowBasedArgumentCollection fbArgs = new FlowBasedArgumentCollection();

    public TrimmedReadsReader(final List<Path> readsFiles, final Path referencePath, final int cloudPrefetchBuffer) {

        final Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer);
        final Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer);

        for ( Path readsFile : readsFiles ) {
            samReaders.add(SamReaderFactory.makeDefault().referenceSequence(referencePath).open(readsFile, cloudWrapper, cloudIndexWrapper));
        }
    }

    public SAMSequenceDictionary getSamSequenceDictionary(final SamReader samReader) {
        return ((samReader != null) ? samReader : samReaders.get(0)).getFileHeader().getSequenceDictionary();
    }

    public Map<SamReader, Collection<FlowBasedRead>>  getReads(final Locatable span, final Locatable vcLoc) {

        final Map<SamReader, Collection<FlowBasedRead>>   readsByReader = new LinkedHashMap<>();
        for ( SamReader samReader : samReaders ) {
            final List<FlowBasedRead>     reads = new LinkedList<>();
            final SAMRecordIterator iter = samReader.query(span.getContig(), span.getStart(), span.getEnd(), false);
            while (iter.hasNext()) {

                // establish record. ignore if variant context is not covered by this read?
                SAMRecord record = iter.next();
                if (!record.contains(vcLoc)) {
                    continue;
                }

                // convert to gatk read
                final String readGroup = record.getReadGroup().getId();
                GATKRead gatkRead = new SAMRecordToGATKReadAdapter(record);

                // filter out?
                if (readFilter != null && !readFilter.test(gatkRead)) {
                    continue;
                }

                // soft/hard clipped bases
                gatkRead = ReadClipper.hardClipSoftClippedBases(gatkRead);
                gatkRead = ReadClipper.hardClipToRegion(gatkRead, span.getStart(), span.getEnd());
                if (gatkRead.isUnmapped() || gatkRead.getCigar().isEmpty())
                    continue;

                // convert to a flow based read
                FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(samReader.getFileHeader(), gatkRead);
                final FlowBasedRead fbr = new FlowBasedRead(gatkRead, rgInfo.flowOrder, rgInfo.maxClass, fbArgs);
                fbr.applyAlignment();

                // clip to given span
                final int read_start = fbr.getStart();
                final int read_end = fbr.getEnd();
                final int diff_left = span.getStart() - read_start;
                final int diff_right = read_end - span.getEnd();
                fbr.applyBaseClipping(Math.max(0, diff_left), Math.max(diff_right, 0), true);

                // check if read is valid. it is possible that read was made invalid by applyBaseClipping
                // if so, ignore it (see FlowBasedRead.java:478 valid_key=false;
                if (!fbr.isValid()) {
                    continue;
                }

                // add to output collection
                reads.add(fbr);
            }
            iter.close();
            readsByReader.put(samReader, reads);
        }

        return readsByReader;
    }

    public SAMFileHeader getHeader(final SamReader samReader) {
        return ((samReader != null) ? samReader : samReaders.get(0)).getFileHeader();
    }

    public void setReadFilter(final CountingReadFilter readFilter) {
        this.readFilter = readFilter;
    }

}
