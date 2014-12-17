package picard.cmdline.programgroups;

import picard.cmdline.CommandLineProgramGroup;

/**
 * @author nhomer
 */
public class Testing implements CommandLineProgramGroup {
    @Override
    public String getName() { return "Unit Testing"; }
    @Override
    public String getDescription() { return "Unit testing"; }
}
