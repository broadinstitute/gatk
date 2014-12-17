package picard.cmdline.programgroups;

import picard.cmdline.CommandLineProgramGroup;

/**
* @author nhomer
*/
public class Metrics implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Metrics"; }
    @Override
    public String getDescription() { return "Tools for reporting metrics on various data types."; }
}
