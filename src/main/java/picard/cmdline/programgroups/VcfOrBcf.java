package picard.cmdline.programgroups;

import picard.cmdline.CommandLineProgramGroup;

/**
* @author nhomer.
*/
public class VcfOrBcf implements CommandLineProgramGroup {
    @Override
    public String getName() { return "VCF/BCF"; }
    @Override
    public String getDescription() { return "Tools for manipulating VCF, BCF, or related data."; }
}
