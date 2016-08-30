package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Unit tests for {@link ContigGermlinePloidyAnnotationTableReader}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class ContigGermlinePloidyAnnotationTableReaderUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/sexgenotyper/";
    private static final File BAD_CONTIG_PLOIDY_ANNOTS_FILE_1 = new File(TEST_SUB_DIR, "contig_annots_bad_autosomal_annot.tsv");
    private static final File BAD_CONTIG_PLOIDY_ANNOTS_FILE_2 = new File(TEST_SUB_DIR, "contig_annots_bad_class.tsv");
    private static final File BAD_CONTIG_PLOIDY_ANNOTS_FILE_3 = new File(TEST_SUB_DIR, "contig_annots_bad_missing_some_annots.tsv");

    /* autosomal contigs can not have different plodies on different classes */
    @Test(expectedExceptions = UserException.BadInput.class)
    public void testBadAutosomalContigPloidyValues() {
        try {
            final ContigGermlinePloidyAnnotationTableReader reader = new ContigGermlinePloidyAnnotationTableReader(
                    new FileReader(BAD_CONTIG_PLOIDY_ANNOTS_FILE_1));
            reader.stream().count();
        } catch (IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read test resource file");
        }
    }

    /* only AUTOSOMAL or ALLOSOMAL is acceptable */
    @Test(expectedExceptions = UserException.BadInput.class)
    public void testBadContigClass() {
        try {
            final ContigGermlinePloidyAnnotationTableReader reader = new ContigGermlinePloidyAnnotationTableReader(
                    new FileReader(BAD_CONTIG_PLOIDY_ANNOTS_FILE_2));
            reader.stream().count();
        } catch (IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read test resource file");
        }
    }

    /* all lines in the table must have all values for the ploidy classes defined in the header */
    @Test(expectedExceptions = UserException.BadInput.class)
    public void testBadMissingSomeAnnotations() {
        try {
            final ContigGermlinePloidyAnnotationTableReader reader = new ContigGermlinePloidyAnnotationTableReader(
                    new FileReader(BAD_CONTIG_PLOIDY_ANNOTS_FILE_3));
            reader.stream().count();
        } catch (IOException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read test resource file");
        }
    }

}
