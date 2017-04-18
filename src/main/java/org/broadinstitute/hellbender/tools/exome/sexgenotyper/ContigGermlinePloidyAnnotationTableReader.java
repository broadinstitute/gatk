package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import com.google.common.collect.Sets;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import javax.annotation.Nonnull;
import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Reads contig germline ploidy annotations from tab-separated files and readers.
 *
 * For example, a basic annotation file for homo sapiens is as follows:
 *
 *     <pre>
 *         CONTIG    CLASS           SEX_XX    SEX_XY
 *         1         AUTOSOMAL       2          2
 *         2         AUTOSOMAL       2          2
 *         ...       ...             ...        ...
 *         X         ALLOSOMAL       2          0
 *         Y         ALLOSOMAL       1          1
 *     </pre>
 *
 * CLASS is either AUTOSOMAL or ALLOSOMAL. "SEX_XX" and "SEX_YY" are arbitrary sex genotype string identifiers
 * along with their ploidies (= number of homologs) for each CONTIG. AUTOSOMAL contigs must have the same ploidy
 * for all sexes. One may include additional sex genotypes (along with contig ploidies) for species having more
 * than two sexes or for detecting aneuploidy by adding additional columns.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class ContigGermlinePloidyAnnotationTableReader extends TableReader<ContigGermlinePloidyAnnotation> {

    private static final Logger logger = LogManager.getLogger(ContigGermlinePloidyAnnotationTableReader.class);

    /**
     * The set of ploidy tags (= sex genotype identifier strings)
     */
    private final Set<String> ploidyTagsSet;

    /**
     * Public constructor from a reader.
     *
     * @param sourceName name of the source
     * @param sourceReader an instance of {@link Reader}
     * @throws IOException if a reading error occurs
     */
    public ContigGermlinePloidyAnnotationTableReader(final String sourceName, @Nonnull final Reader sourceReader)
            throws IOException {
        super(sourceName, sourceReader);
        TableUtils.checkMandatoryColumns(columns(), ContigGermlinePloidyAnnotationTableColumn.MANDATORY_CONTIG_ANNOTATION_COLUMNS,
                UserException.BadInput::new);
        ploidyTagsSet = Sets.difference(new HashSet<>(columns().names()),
                ContigGermlinePloidyAnnotationTableColumn.MANDATORY_CONTIG_ANNOTATION_COLUMNS_SET);
        if (ploidyTagsSet.isEmpty()) {
            throw new UserException.BadInput("At least one ploidy column is required!");
        }
        logger.info("Ploidy tags: " + ploidyTagsSet.stream().collect(Collectors.joining(", ")));
    }

    /**
     * Public constructor from a reader.
     *
     * @param sourceReader an instance of {@link Reader}
     * @throws IOException if a reading error occurs
     */
    public ContigGermlinePloidyAnnotationTableReader(@Nonnull final Reader sourceReader) throws IOException {
        this(null, sourceReader);
    }

    /**
     * Creates a {@link ContigGermlinePloidyAnnotation} instance from a {@link DataLine}
     *
     * @param dataLine a data line
     * @return an instance of {@link ContigGermlinePloidyAnnotation}
     */
    @Override
    protected ContigGermlinePloidyAnnotation createRecord(@Nonnull final DataLine dataLine) {
        final String contigName = dataLine.get(ContigGermlinePloidyAnnotationTableColumn.CONTIG);

        final ContigClass contigClass;
        final String contigClassString = dataLine.get(ContigGermlinePloidyAnnotationTableColumn.CLASS);
        if (!ContigClass.CONTIG_CLASS_NAMES_SET.contains(contigClassString)) {
            throw new UserException.BadInput("Bad contig class: provided value: " + contigClassString + ", acceptable values: " +
                    ContigClass.CONTIG_CLASS_NAMES_SET.stream().collect(Collectors.joining(", ", "[", "]")));
        } else {
            contigClass = ContigClass.valueOf(contigClassString);
        }

        final Map<String, Integer> ploidyMap = new HashMap<>();
        try {
            /* all lines must have all ploidy annotations defined in the header */
            ploidyTagsSet.forEach(tag -> ploidyMap.put(tag, dataLine.getInt(tag)));
        } catch (final IllegalArgumentException ex) {
            throw new UserException.BadInput("All lines in the contig germline ploidy annotation table must have values for all" +
                    " ploidy classes; " + ex.getMessage());
        }

        return new ContigGermlinePloidyAnnotation(contigName, contigClass, ploidyMap);
    }

    /**
     * Reads contig ploidy annotations from a file.
     *
     * @param contigGermlinePloidyAnnotationsFile a file containing contig germline ploidy annotations
     * @return a list of contig ploidy annotations
     * @throws IOException if a read error occurs
     */
    public static List<ContigGermlinePloidyAnnotation> readContigGermlinePloidyAnnotationsFromFile(@Nonnull final File contigGermlinePloidyAnnotationsFile) {
        IOUtils.canReadFile(contigGermlinePloidyAnnotationsFile);
        try {
            return readContigGermlinePloidyAnnotationsFromReader(contigGermlinePloidyAnnotationsFile.getAbsolutePath(),
                    new FileReader(contigGermlinePloidyAnnotationsFile));
        } catch (final FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile("Could not read contig germline ploidy annotations file " +
                    contigGermlinePloidyAnnotationsFile.getAbsolutePath());
        }
    }

    /**
     * Reads contig germline ploidy annotations from a {@link Reader}.
     *
     * @param contigAnnotationSourceName a string identifier for the reader.
     * @param contigAnnotationReader an instance of {@link Reader}.
     * @return list of contig ploidy annotations
     * @throws IOException if a read error occurs
     */
    public static List<ContigGermlinePloidyAnnotation> readContigGermlinePloidyAnnotationsFromReader(@Nonnull final String contigAnnotationSourceName,
                                                                                                     @Nonnull final Reader contigAnnotationReader) {
        /* read contig annotations */
        try (final ContigGermlinePloidyAnnotationTableReader reader =
                     new ContigGermlinePloidyAnnotationTableReader(contigAnnotationSourceName, contigAnnotationReader)) {
            return reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(contigAnnotationSourceName, e);
        }
    }
}
