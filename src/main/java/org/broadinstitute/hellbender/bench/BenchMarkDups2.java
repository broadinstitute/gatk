package org.broadinstitute.hellbender.bench;

import com.google.api.services.genomics.model.Read;
import com.google.common.base.Stopwatch;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.tools.dataflow.transforms.markduplicates.MarkDuplicatesFromShardsDataflowTransform;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

/**
 *
 */
public class BenchMarkDups2 {

    public static void main(String[] args) {
        String fname = "../../sample-data/CEUTrio.HiSeq.WGS.b37.ch20.4m-8m.NA12878.bam";
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
            .enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES);
        final SamReader samReader = samReaderFactory.open(new File(fname));

        SAMFileHeader header = samReader.getFileHeader();


        int COUNT = 1_000_000;
        System.out.println("Loading "+COUNT+" reads from "+fname);

        SAMRecord[] records = new SAMRecord[COUNT];
        ArrayList<GATKRead> gatkReads = new ArrayList<GATKRead>(COUNT);
        Read[] reads = new Read[COUNT];
        int i=0;
        Stopwatch readFromDisk = Stopwatch.createStarted();
        for (SAMRecord r : samReader) {
            records[i++] = r;
            if (i>=COUNT) break;
        }
        if (i!=COUNT) {
            System.out.println("Only "+i+" records in that file, pick something larger.");
            return;
        }
        readFromDisk.stop();
        System.out.println("" + readFromDisk.elapsed(TimeUnit.MILLISECONDS) + " ms to read from disk.");
        Stopwatch adapt = Stopwatch.createStarted();
        for (i=0; i<COUNT; i++) {
            gatkReads.add(new SAMRecordToGATKReadAdapter(records[i]));
        }
        System.out.println("" + adapt.elapsed(TimeUnit.MILLISECONDS) + " ms to wrap adapters around the SAMRecords.");

        Stopwatch markDups  = Stopwatch.createStarted();
        MarkDuplicatesFromShardsDataflowTransform.localApply(header, null, gatkReads);
        System.out.println("" + markDups.elapsed(TimeUnit.MILLISECONDS) + " ms to mark duplicates.");


    }
}
