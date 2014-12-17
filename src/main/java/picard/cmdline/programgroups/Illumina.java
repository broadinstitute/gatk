package picard.cmdline.programgroups;

import picard.cmdline.CommandLineProgramGroup;

/**
* @author nhomer
*/
public class Illumina implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Illumina Tools"; }
    @Override
    public String getDescription() { return "Tools for manipulating data specific to Illumina sequencers."; }
}
