package org.broadinstitute.hellbender.engine;

import htsjdk.samtools.SAMRecord;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.Spliterators;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

//HACK to get the code to compile and attempt to get the tests to run
public final class HACKRefMetaDataTracker {
    private final GenomeLocParser glp;

    public HACKRefMetaDataTracker(GenomeLocParser glp){
        this.glp = glp;
    }

    public List<Feature> getValues(File knownSitesVCF, SAMRecord read) {
        GenomeLoc loc = glp.createGenomeLoc(read);
        Iterator<VariantContext> variantC = new VCFFileReader(knownSitesVCF).query(loc.getContig(), loc.getStart(), loc.getStop());
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(variantC, 0), false)
                .collect(Collectors.toList());
    }
}
