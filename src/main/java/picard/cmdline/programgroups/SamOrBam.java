package picard.cmdline.programgroups;

import picard.cmdline.CommandLineProgramGroup;

/**
* @author nhomer
*/
public class SamOrBam implements CommandLineProgramGroup {
    @Override
    public String getName() { return "SAM/BAM"; }
    @Override
    public String getDescription() { return "Tools for manipulating SAM, BAM, or related data."; }
}
