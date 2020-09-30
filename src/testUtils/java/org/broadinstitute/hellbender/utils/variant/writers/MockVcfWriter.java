package org.broadinstitute.hellbender.utils.variant.writers;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import java.util.ArrayList;
import java.util.List;

final class MockVcfWriter implements VariantContextWriter {
    final List<VariantContext> emitted = new ArrayList<>();
    boolean headerWritten = false;
    boolean closed = false;
    boolean error = false;
    boolean headerSet = false;

    @Override
    public void writeHeader(VCFHeader header) {
        headerSet = true;
        headerWritten = true;
    }

    @Override
    public void close() {
        closed = true;
    }

    @Override
    public boolean checkError() {
        return error;
    }

    @Override
    public void add(VariantContext vc) {
        emitted.add(vc);
    }

    @Override
    public void setHeader(VCFHeader header) {
        headerSet = true;
    }
}
