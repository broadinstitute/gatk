package picard.cmdline;

import java.util.Comparator;

/**
 * Interface for groups of CommandLinePrograms.
 * @author Nils Homer
 */
public interface CommandLineProgramGroup {

    /** Gets the name of this program. **/
    public String getName();
    /** Gets the description of this program. **/
    public String getDescription();
    /** Compares two program groups by name. **/
    static public Comparator<CommandLineProgramGroup> comparator = new Comparator<CommandLineProgramGroup>() {
        public int compare(final CommandLineProgramGroup a, final CommandLineProgramGroup b) {
            return a.getName().compareTo(b.getName());
        }
    };
}
