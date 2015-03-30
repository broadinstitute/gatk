package org.broadinstitute.hellbender.utils;


import com.google.api.services.genomics.model.Read;
import com.google.cloud.genomics.dataflow.readers.bam.ReadConverter;
import com.google.cloud.genomics.gatk.common.GenomicsConverter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

public class ReadConverterTest extends BaseTest {



    @DataProvider(name = "files")
    public Object[][] files(){
      return new Object[][] {
              {"org/broadinstitute/hellbender/tools/flag_stat.sam"},
              {"org/broadinstitute/hellbender/tools/flag_stat.bam"},
              {"org/broadinstitute/hellbender/tools/count_reads_sorted.bam"}
      };
    }

    @Test(dataProvider = "files")
    public void SamToReadToSamTest(String filePath) throws IOException {
        File samInput = new File(publicTestDir, filePath);
        ReadsDataSource reads = new ReadsDataSource(samInput);
        SAMFileHeader header = reads.getHeader();
       // File output = createTempFile(samInput.toString(),"" );

        //try(SAMFileWriter writer = new SAMFileWriterFactory().makeSAMWriter(header, true, output)){
            for (SAMRecord sam : reads){
                Read read = ReadConverter.makeRead(sam);
                SAMRecord newSam = GenomicsConverter.makeSAMRecord(read, header );
              //  writer.addAlignment(newSam);
              //  System.out.println(sam.getSAMString() + newSam.getSAMString());
                Assert.assertEquals(newSam.getSAMString(), sam.getSAMString());
                if(!sam.equals(newSam)){
                  System.out.println("sam != newSam");
                }
                if(!newSam.equals(sam)){
                  System.out.println("newSam != sam");
                }
                Assert.assertEquals(newSam, sam);
            }
       // }
        //IntegrationTestSpec.assertMatchingFiles(Collections.singletonList(output), Collections.singletonList(samInput.toString()));

    }

  @Test
  public void samBAM() throws IOException {
    File samInput = new File(publicTestDir, "org/broadinstitute/hellbender/tools/flag_stat.sam");
    File bamInput = new File(publicTestDir, "org/broadinstitute/hellbender/tools/flag_stat.bam");
    SamReader samReads = SamReaderFactory.makeDefault().open(samInput);
    SamReader bamReads = SamReaderFactory.makeDefault().open(bamInput);

    Iterator<SAMRecord> bamIterator = bamReads.iterator();

    for (SAMRecord sam : samReads){
      SAMRecord bam = bamIterator.next();
      Assert.assertEquals(sam.toString(), bam.toString());
      Assert.assertEquals(sam, bam);

    }
  }

}



