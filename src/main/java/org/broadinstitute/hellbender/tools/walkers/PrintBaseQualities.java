package org.broadinstitute.hellbender.tools.walkers;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.hadoop.yarn.webapp.hamlet.Hamlet;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

@CommandLineProgramProperties(
        summary = "asdf",
        oneLineSummary = "asdf",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
public class PrintBaseQualities extends ReadWalker {
    PrintStream ps;
    int count = 0;

    @Override
    public void onTraversalStart(){
        final File file = new File("bqs.tsv");
        try {
            ps = new PrintStream(file);
            ps.println("index\tread\tcycle\tquality");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }


    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        if (count > 1000){
            return;
        }
        final String read1 = read.isFirstOfPair() ? "read1" : "read2";
        final boolean reverseRead = read.isReverseStrand();

        final byte[] baseQualities = read.getBaseQualities();
        if (reverseRead){
            ArrayUtils.reverse(baseQualities);
        }

        for (int i = 0; i < baseQualities.length; i++){
            ps.println(count + "\t" + read1 + "\t" + i + "\t" + (int)baseQualities[i]);
        }
        count++;
    }

    @Override
    public void closeTool(){
        ps.close();
    }
}
