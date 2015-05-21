package org.broadinstitute.hellbender.dev.pipelines.bqsr;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.coders.StringUtf8Coder;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.join.CoGbkResult;
import com.google.cloud.dataflow.sdk.transforms.join.CoGroupByKey;
import com.google.cloud.dataflow.sdk.transforms.join.KeyedPCollectionTuple;
import com.google.cloud.dataflow.sdk.util.GcsUtil;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import com.google.cloud.dataflow.sdk.values.*;
import com.google.cloud.genomics.dataflow.coders.GenericJsonCoder;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.base.Stopwatch;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Iterables;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringLineReader;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibrationArgumentCollection;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.tools.walkers.bqsr.RecalibrationEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.io.*;
import java.nio.channels.Channels;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 * Base Quality Score Recalibration, phase 1
 * This contains the Dataflow-specific bits.
 */
public final class BQSR_Dataflow implements Serializable {

  // how many bases until we move on to the next shard
  static final int BLOCKSIZE = 1000000;
  public static final TupleTag<Read> readTag = new TupleTag<>();
  public static final TupleTag<SimpleInterval> intervalTag = new TupleTag<>();
  public static final TupleTag<RecalibrationTables> tablesTag = new TupleTag<>();

  /**
   * Get a single RecalibrationTables object that represents the output of phase 1 of BQSR.
   * <p>
   * The reference file (*.fasta) must also have a .dict and a .fasta.fai next to it.
   */
  public static PCollection<RecalibrationTables> GetRecalibrationTables(SAMFileHeader readsHeader, PCollection<Read> reads, String referenceFileName, BaseRecalibrationArgumentCollection toolArgs, PCollection<SimpleInterval> placesToIgnore) {
    PCollection<RecalibrationTables> stats = computeBlockStatistics(readsHeader, referenceFileName, toolArgs, groupByBlock(reads, placesToIgnore));
    PCollection<RecalibrationTables> oneStat = aggregateStatistics(stats);
    return oneStat;
  }

  /**
   * Throws an exception if any of the files is missing.
   * (offlineauth can be null if the files are local)
   */
  public static void ensureReferenceIsReadable(final PipelineOptions popts, String filename) throws IOException, GeneralSecurityException {
    for (String fn : getRelatedFiles(filename)) {
      if (BucketUtils.isCloudStorageUrl(fn)) {
        // make sure we can access those files on Google Cloud Storage.
        GcsPath path = GcsPath.fromUri(fn);
        // this will throw if we can't get to the file.
        new GcsUtil.GcsUtilFactory().create(popts).fileSize(path);
      } else {
        // make sure we can access those files on the local filesystem.
        if (!new File(fn).canRead()) {
          if (!new File(fn).exists()) {
            throw new IOException("File not found: " + fn);
          } else {
            throw new IOException("File present, but not readable: " + fn);
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------------------------------
  // non-public methods

  private static String[] getRelatedFiles(String fastaFilename) {
    return new String[] { fastaFilename, ReferenceUtils.getFastaDictionaryFileName(fastaFilename), ReferenceUtils.getFastaIndexFileName(fastaFilename) };
  }

  // a key to group intervals by (aka sharding)
  private static String posKey(final Locatable loc) {
    return loc.getContig() + (loc.getStart() / BLOCKSIZE);
  }

  // a key to group intervals by (aka sharding)
  private static String posKey(final Read r) {
    // Reads are 0-based, SimpleIntervals are 1-based
    int start = r.getAlignment().getPosition().getPosition().intValue() + 1;
    return posKey(new SimpleInterval(r.getAlignment().getPosition().getReferenceName(),
            start,
            start));
  }

  // the shard that the last base of the interval lands in,
  // if different from the first base's. Null otherwise.
  private static String posKeyEnd(final Locatable loc) {
    if (loc.getStart()/BLOCKSIZE == (loc.getEnd())/BLOCKSIZE) return null;
    return loc.getContig() + (loc.getEnd() / BLOCKSIZE);
  }

  /**
   * returns union(reads,placestoIgnore).groupBy(x->posKey(x))
   */
  private static PCollection<KV<String, CoGbkResult>> groupByBlock(final PCollection<Read> reads, final PCollection<SimpleInterval> placesToIgnore) {
    // shard reads
    PCollection<KV<String, Read>> shardedReads = reads.apply(ParDo
        .named("shard reads")
        .of(
            new DoFn<Read, KV<String, Read>>() {
              @Override
              public void processElement(ProcessContext c) {
                Read r = c.element();
                c.output(KV.of(posKey(r), r));
              }

            }))
        // Dataflow boilerplate
        .setCoder(KvCoder.of(StringUtf8Coder.of(), GenericJsonCoder.of(Read.class)));
    // send ignores to every shard that overlaps with them (so, possibly more than one).
    // (for now we assume they can't span more than two)
    PCollection<KV<String, SimpleInterval>> shardedIgnore = placesToIgnore.apply(ParDo
        .named("shard known intervals")
        .of(
            new DoFn<SimpleInterval, KV<String, SimpleInterval>>() {
              @Override
              public void processElement(ProcessContext c) {
                SimpleInterval i = c.element();
                String firstPos = posKey(i);
                c.output(KV.of(firstPos, i));
                String secondPos = posKeyEnd(i);
                if (null!=secondPos) {
                  c.output(KV.of(secondPos, i));
                }
              }
            }))
        // Dataflow boilerplate
        .setCoder(KvCoder.of(StringUtf8Coder.of(), SerializableCoder.of(SimpleInterval.class)));

    return KeyedPCollectionTuple.of(readTag, shardedReads)
        .and(intervalTag, shardedIgnore)
        .apply(CoGroupByKey.<String>create());
  }

  /**
   * For each input block, compute its statistics.
   * <p>
   * At a high level, just delegate to the CalibrationTablesBuilder.
   * More specifically, this downloads the reference files, sorts the known intervals for each shard,
   * and then delegates to CalibrationTablesBuilder. Then it outputs a log message about timings.
   */
  private static PCollection<RecalibrationTables> computeBlockStatistics(final SAMFileHeader readsHeader, String referenceFileName, final BaseRecalibrationArgumentCollection toolArgs, final PCollection<KV<String, CoGbkResult>> readsAndIgnores) {
    PCollection<RecalibrationTables> ret = readsAndIgnores.apply(ParDo
            .named("computeBlockStatistics")
            .of(new DoFn<KV<String, CoGbkResult>, RecalibrationTables>() {
              CalibrationTablesBuilder ct;
              Stopwatch timer;
              int nBlocks = 0;
              int nReads = 0;

              @Override
              public void startBundle(DoFn.Context c) throws Exception {
                timer = Stopwatch.createStarted();
                SAMFileHeader header = readsHeader;

                String localReference = referenceFileName;
                if (BucketUtils.isCloudStorageUrl(referenceFileName)) {
                  // the reference is on GCS, download all 3 files locally first.
                  System.out.println("Downloading reference files");
                  boolean first = true;
                  for (String fname : getRelatedFiles(referenceFileName)) {
                    String localName = "reference";
                    int slash = fname.lastIndexOf('/');
                    if (slash >= 0) {
                      localName = fname.substring(slash + 1);
                    }
                    // download reference if necessary
                    if (new File(localName).exists()) {
                      if (first) localReference = localName;
                    } else {
                      try (
                              InputStream in = Channels.newInputStream(new GcsUtil.GcsUtilFactory().create(c.getPipelineOptions()).open(GcsPath.fromUri(fname)));
                              FileOutputStream fout = new FileOutputStream(localName)) {
                        final byte[] buf = new byte[1024 * 1024];
                        int count;
                        while ((count = in.read(buf)) > 0) {
                          fout.write(buf, 0, count);
                        }
                      }
                      if (first) localReference = localName;
                    }
                    first = false;
                  }
                  System.out.printf("Done downloading reference files (%s ms)\n", timer.elapsed(TimeUnit.MILLISECONDS));
                }

                ct = new CalibrationTablesBuilder(header, localReference, toolArgs);
              }

              @Override
              public void processElement(ProcessContext c) throws Exception {
                nBlocks++;
                // get the reads
                KV<String, CoGbkResult> e = c.element();
                List<Read> reads = new ArrayList<Read>();
                Iterable<Read> readsIter = e.getValue().getAll(BQSR_Dataflow.readTag);
                Iterables.addAll(reads, readsIter);
                nReads += reads.size();
                // get the skip intervals
                List<SimpleInterval> skipIntervals = new ArrayList<SimpleInterval>();
                Iterables.addAll(skipIntervals, e.getValue().getAll(BQSR_Dataflow.intervalTag));
                Collections.sort(skipIntervals, new Comparator<SimpleInterval>() {
                  @Override
                  public int compare(SimpleInterval o1, SimpleInterval o2) {
                    return ComparisonChain.start()
                            .compare(o1.getContig(), o2.getContig())
                            .compare(o1.getStart(), o2.getStart())
                            .result();
                  }
                });
                // update our statistics
                ct.add(reads, skipIntervals);
              }

              @Override
              public void finishBundle(DoFn.Context c) throws Exception {
                ct.done();
                c.output(ct.getRecalibrationTables());
                System.out.println("Finishing a block statistics bundle. It took " + timer.elapsed(TimeUnit.MILLISECONDS) + " ms to process " + nBlocks + " blocks, " + nReads + " reads.");
            timer = null;
          }
        }));

    ret.setCoder(SerializableCoder.of(RecalibrationTables.class));
    return ret;
  }

  /**
   * Merge the statistics from each block. The resulting "collection" contains a single element, with the answer.
   */
  private static PCollection<RecalibrationTables> aggregateStatistics(final PCollection<RecalibrationTables> tables) {
    return tables
        // aggregate
        .apply(Combine.globally(new RecalibrationTablesMerger()))
        // call finalize on the result
        .apply(ParDo
            .named("finalizeRecalTables")
            .of(new DoFn<RecalibrationTables, RecalibrationTables>() {
              @Override
              public void processElement(ProcessContext c) throws Exception {
                RecalibrationTables tables = c.element();
                RecalibrationEngine.finalizeRecalibrationTables(tables);
                c.output(tables);
              }
            }));
  }

}


