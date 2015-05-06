package org.broadinstitute.hellbender.dev.utils;

import com.google.api.services.genomics.model.LinearAlignment;
import com.google.api.services.genomics.model.Read;
import com.google.api.services.genomics.model.Position;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.Test;


/**
 * A few tests for the ReadRecord class. We need more, eventually.
 */
public class ReadRecordTest {

  @Test
  public void simpleRead() {
    Read r = new Read()
        .setNumberReads(1)
        .setReadNumber(0)
        .setReadGroupId("RG1")
        .setAlignedSequence("ATCG")
        .setAlignment(new LinearAlignment().setPosition(new Position().setPosition(10L)));
    SAMFileHeader header = new SAMFileHeader();
    // convert from 0-based Read to 1-based SAMRecord
    SAMRecord record = new ReadRecord(r, header);

    Assert.assertEquals(record.getAlignmentStart(), 11);
    Assert.assertEquals(record.getReadLength(), 4);
    Assert.assertEquals(record.getReadBases(), new byte[] { 'A', 'T', 'C', 'G' });
    Assert.assertEquals(record.getSupplementaryAlignmentFlag(), false);
  }

  @Test
  public void paired() {
    Read r = new Read()
        .setReadNumber(0)
        .setNumberReads(2);
    SAMFileHeader header = new SAMFileHeader();

    SAMRecord record = new ReadRecord(r, header);

    Assert.assertEquals(record.getReadPairedFlag(), true);
    Assert.assertEquals(record.getFirstOfPairFlag(), true);
  }

  @Test
  public void secondOfPair() {
    Read r = new Read()
        .setReadNumber(1)
        .setNumberReads(2);
    SAMFileHeader header = new SAMFileHeader();

    SAMRecord record = new ReadRecord(r, header);

    Assert.assertEquals(record.getReadPairedFlag(), true);
    Assert.assertEquals(record.getFirstOfPairFlag(), false);
  }

  @Test
  public void unpaired() {
    Read r = new Read()
        .setReadNumber(0)
        .setNumberReads(1);
    SAMFileHeader header = new SAMFileHeader();

    SAMRecord record = new ReadRecord(r, header);

    Assert.assertEquals(record.getReadPairedFlag(), false);
  }

  @Test
  public void supplementaryAlignment() {
    Read r = new Read()
        .setNumberReads(1)
        .setReadNumber(0)
        .setReadGroupId("RG1")
        .setAlignedSequence("ATCG")
        .setAlignment(new LinearAlignment().setPosition(new Position().setPosition(10L)))
        .setSupplementaryAlignment(true);
    SAMFileHeader header = new SAMFileHeader();

    SAMRecord record = new ReadRecord(r, header);

    Assert.assertEquals(record.getSupplementaryAlignmentFlag(), true);
  }

}