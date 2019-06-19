package org.broadinstitute.hellbender.tools.evoquer;

import htsjdk.samtools.util.CloseableIterator;
import org.apache.avro.generic.GenericRecord;

public interface GATKAvroReader extends Iterable<GenericRecord>, CloseableIterator<GenericRecord> {
    org.apache.avro.Schema getSchema();
}