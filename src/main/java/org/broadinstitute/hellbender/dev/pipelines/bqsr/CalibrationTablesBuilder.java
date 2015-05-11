
package org.broadinstitute.hellbender.dev.pipelines.bqsr;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.genomics.gatk.common.GenomicsConverter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.tribble.Feature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibratorUprooted;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

/*
 * Feed me reads and skipIntervals, and I'll get you a RecalibrationTables.
 *
 * (OK, yes, I need the reference genome too).
 *
 * This class wraps around BaseRecalibratorUprooted, giving it a cozy environment outside of the normal
 * CommandLine framework. It's also completely independent of Dataflow, so it could be tested outside of that framework.
 */
public class CalibrationTablesBuilder {

  private BaseRecalibratorUprooted br;
  private SAMFileHeader header;
  private ReferenceDataSource reference;

  public CalibrationTablesBuilder(SAMFileHeader header, String referenceFileName, RecalibrationArgumentCollection toolArgs) {
    // 1. we'll paste the header onto the new SAMRecord objects, the code needs it.
    this.header = header;
    // 2. the reference we copied
    File refFile = new File(referenceFileName);
    reference = new ReferenceDataSource(refFile);
    // 3. create the class that'll do the actual work
    br = new BaseRecalibratorUprooted(header, toolArgs);
    br.onTraversalStart(refFile);
  }

  /**
   * For the reference, we need three files: dict, fasta, and fai. Typically the user only specifies one file name, e.g.:
   * input-hs37d5.fa
   * and then we just know to also look for
   * input-hs37d5.dict and input-hs37d5.fa.fai
   *
   * This function takes the user-provided name and returns the three names we have to care about.
   * The first one is guaranteed to be the one given in input.
   * This function also works if the file names start with gs://
   */
  public static String[] expandReferenceFilename(String fastaFilename) {
    String[] ret = new String[3];
    ret[0] = fastaFilename;
    ret[1] = fastaFilename + ".fai";
    int lastDot = fastaFilename.lastIndexOf('.');
    ret[2] = fastaFilename.substring(0, lastDot) + ".dict";
    return ret;
  }

    /**
     * call this as many times as you want.
     * The skipIntervals have to be 1-based, close-ended.
     */
  public void add(Iterable<Read> reads, List<SimpleInterval> skipIntervals) {
    if (null==br) throw new RuntimeException("Can't call add after done");

    for (Read r : reads) {
      SAMRecord sr = GenomicsConverter.makeSAMRecord(r, header);
      final SimpleInterval readInterval = sr.getReadUnmappedFlag() ? null :
          new SimpleInterval(sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd());
      // TODO: this could probably be sped up by taking advantage of a sorted order.
      List<SimpleInterval> knownSitesOverlappingReadInterval =
              skipIntervals.stream().filter(x -> x.overlaps(readInterval))
                      .collect(Collectors.toList());
      br.apply(sr, new ReferenceContext(reference, readInterval), knownSitesOverlappingReadInterval);
    }
  }

    /**
     * call this once before destroying this object.
     */
  public RecalibrationTables done() {
    RecalibrationTables ret = br.getTables();
    br = null;
    return ret;
  }

}
