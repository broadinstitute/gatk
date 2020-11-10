package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import java.io.Closeable;
import java.io.IOException;

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.conf.Configuration;
import org.apache.orc.*;
import org.apache.orc.storage.ql.exec.vector.BytesColumnVector;
import org.apache.orc.storage.ql.exec.vector.LongColumnVector;
import org.apache.orc.storage.ql.exec.vector.VectorizedRowBatch;

public class PetOrcWriter implements Closeable {
    private Writer writer;
    private VectorizedRowBatch batch;
    private LongColumnVector locationCV;
    private LongColumnVector sampleCV;
    private BytesColumnVector stateCV;

    static private TypeDescription schema = TypeDescription.fromString("struct<location:int,sample:int,state:char(1)>");

    public PetOrcWriter(String outputFile) throws IOException{
        Configuration conf = new Configuration();
        
        this.writer = OrcFile.createWriter(new Path(outputFile),
                      OrcFile.writerOptions(conf)
                             .setSchema(schema)
                             .overwrite(true)
                             .compress(CompressionKind.SNAPPY)
                             .stripeSize( 16 * 1024 * 1024 )
                             .blockSize(   1 * 1024 * 1024 )
                             );
        this.batch = schema.createRowBatch();
        this.locationCV = (LongColumnVector) batch.cols[0];
        this.sampleCV = (LongColumnVector) batch.cols[1];
        this.stateCV = (BytesColumnVector) batch.cols[2];         
    }

    public void addRow(long location, long sampleId, String state) throws IOException {       
        int row = batch.size++;
        locationCV.vector[row] = location;
        sampleCV.vector[row] = sampleId;
        stateCV.setVal(row, state.getBytes());

        // If the batch is full, write it out and reset the counter
        if (batch.size == batch.getMaxSize()) {            
            writer.addRowBatch(batch);
            batch.reset();
        }
    }

    public void close() throws IOException {
        if (batch.size != 0) {
            writer.addRowBatch(batch);
            batch.reset();
        }
        writer.close();
    }
}
