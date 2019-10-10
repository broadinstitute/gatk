package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class UMIReadFilter extends ReadFilter {
    static final long serialVersionUID = 1L;
    public static final String UMI_TAG = "RX";

    @Argument(fullName= ReadFilterArgumentDefinitions.UMI_NAME, doc = "umi must have the format 'XXX-XXX'", optional=false)
    public String umi;

    public UMIReadFilter() {};

    @Override
    public boolean test(GATKRead read) {
        final String readUMI = read.getAttributeAsString(UMI_TAG);
        final String umi1 = readUMI.split("-", 2)[0];
        final String umi2 = readUMI.split("-", 2)[1];

        return umi.equals(umi1 + "-" + umi2) || umi.equals(umi2 + "-" + umi1);
    }
}
