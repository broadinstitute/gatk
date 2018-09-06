package org.broadinstitute.hellbender.utils.recalibration.covariates;

import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.ClassFinder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents a list of BQSR covariates.
 *
 * Note: the first two covariates ({@link ReadGroupCovariate} and {@link QualityScoreCovariate})
 * are special in the way that they are represented in the BQSR recalibration table, and are
 * always required.
 *
 * The remaining covariates are called "additional covariates". The default additional covariates
 * are the context and cycle covariates, but the client can request others and/or disable the default
 * additional covariates.
 */
public final class BQSRCovariateList implements Iterable<Covariate>, Serializable {
    private static final long serialVersionUID = 1L;

    private static final Logger logger = LogManager.getLogger(BQSRCovariateList.class);

    private final ReadGroupCovariate readGroupCovariate;
    private final QualityScoreCovariate qualityScoreCovariate;
    private final List<Covariate> additionalCovariates;
    private final List<Covariate> allCovariates;

    private final Map<Class<? extends Covariate>, Integer> indexByClass;

    private static final List<String> REQUIRED_COVARIATE_NAMES =
            Collections.unmodifiableList(Arrays.asList("ReadGroupCovariate", "QualityScoreCovariate"));

    private static final List<String> STANDARD_COVARIATE_NAMES =
            Collections.unmodifiableList(Arrays.asList("ContextCovariate", "CycleCovariate"));

    public static final List<String> COVARIATE_PACKAGES =
            Collections.unmodifiableList(Arrays.asList("org.broadinstitute.hellbender.utils.recalibration.covariates"));
    public static final Set<Class<?>> DISCOVERED_COVARIATES;

    static {
        final ClassFinder classFinder = new ClassFinder();

        for ( final String covariatePackage : COVARIATE_PACKAGES ) {
            classFinder.find(covariatePackage, Covariate.class);
        }

        DISCOVERED_COVARIATES = Collections.unmodifiableSet(classFinder.getConcreteClasses());
    }

    public static List<String> getAllDiscoveredCovariateNames() {
        return DISCOVERED_COVARIATES.stream().map(Class::getSimpleName).collect(Collectors.toList());
    }

    /**
     * Creates a new list of BQSR covariates and initializes each covariate.
     */
    public BQSRCovariateList(final RecalibrationArgumentCollection rac, final SAMFileHeader header) {
        this(rac, ReadGroupCovariate.getReadGroupIDs(header));
    }

    /**
     * Creates a new list of BQSR covariates and initializes each covariate.
     */
    public BQSRCovariateList(final RecalibrationArgumentCollection rac, final List<String> allReadGroups) {
        readGroupCovariate = new ReadGroupCovariate();
        readGroupCovariate.initialize(rac, allReadGroups);
        qualityScoreCovariate = new QualityScoreCovariate();
        qualityScoreCovariate.initialize(rac, allReadGroups);

        this.additionalCovariates = Collections.unmodifiableList(createAdditionalCovariates(rac, allReadGroups));

        List<Covariate> allCovariatesList = new ArrayList<>();
        allCovariatesList.add(readGroupCovariate);
        allCovariatesList.add(qualityScoreCovariate);
        additionalCovariates.forEach(allCovariatesList::add);
        this.allCovariates = Collections.unmodifiableList(allCovariatesList);

        //precompute for faster lookup (shows up on profile)
        indexByClass = new LinkedHashMap<>();
        for(int i = 0; i < allCovariates.size(); i++){
            indexByClass.put(allCovariates.get(i).getClass(), i);
        }
    }

    public static boolean isRequiredCovariate(final String covariateName) {
        return REQUIRED_COVARIATE_NAMES.contains(covariateName);
    }

    public static boolean isStandardCovariate(final String covariateName) {
        return STANDARD_COVARIATE_NAMES.contains(covariateName);
    }

    private List<Covariate> createAdditionalCovariates(final RecalibrationArgumentCollection rac, final List<String> allReadGroups) {
        final List<Covariate> additionalCovariatesToAdd = new ArrayList<>();

        if ( ! rac.DO_NOT_USE_STANDARD_COVARIATES ) {
            additionalCovariatesToAdd.addAll(createStandardCovariates(rac, allReadGroups));
        }

        for ( final String requestedAdditionalCovariate : rac.COVARIATES ) {
            if ( isRequiredCovariate(requestedAdditionalCovariate) ) {
                logger.warn("Covariate " + requestedAdditionalCovariate + " is a required covariate that is always on. Ignoring explicit request for it.");
            }
            else if ( ! rac.DO_NOT_USE_STANDARD_COVARIATES && isStandardCovariate(requestedAdditionalCovariate) ) {
                logger.warn("Covariate " + requestedAdditionalCovariate + " is a standard covariate that is always on when not running with --no-standard-covariates. Ignoring explicit request for it.");
            }
            else {
                if ( additionalCovariatesToAdd.stream().anyMatch(cov -> cov.getClass().getSimpleName().equals(requestedAdditionalCovariate)) ) {
                    throw new UserException("Covariate " + requestedAdditionalCovariate + " was requested multiple times");
                }
                
                additionalCovariatesToAdd.add(createAdditionalCovariate(requestedAdditionalCovariate, rac, allReadGroups));
            }
        }

        return additionalCovariatesToAdd;
    }

    private List<Covariate> createStandardCovariates(final RecalibrationArgumentCollection rac, final List<String> allReadGroups) {
        List<Covariate> standardCovariates = new ArrayList<>();

        for ( final String standardCovariateName : STANDARD_COVARIATE_NAMES ) {
            standardCovariates.add(createAdditionalCovariate(standardCovariateName, rac, allReadGroups));
        }

        return standardCovariates;
    }

    private Covariate createAdditionalCovariate(final String covariateName, final RecalibrationArgumentCollection rac, final List<String> allReadGroups) {
        for ( final Class<?> covariateClass : DISCOVERED_COVARIATES ) {
            if ( covariateName.equals(covariateClass.getSimpleName()) ) {
                try {
                    @SuppressWarnings("unchecked")
                    final Covariate covariate = ((Class<? extends Covariate>)covariateClass).newInstance();
                    covariate.initialize(rac, allReadGroups);
                    return covariate;
                }
                catch ( Exception e ) {
                    throw new GATKException("Error instantiating covariate class " + covariateClass.getSimpleName());
                }
            }
        }

        throw new UserException("No covariate with the name " + covariateName + " was found. " +
                "Available covariates are: " + DISCOVERED_COVARIATES);
    }

    /**
     * Returns 2. ReadGroupCovariate and QualityScoreCovariate are always required
     */
    public static int numberOfRequiredCovariates() {
        return REQUIRED_COVARIATE_NAMES.size();
    }

    /**
     * Returns the list of simple class names of our covariates. The returned list is unmodifiable.
     * For example "CycleCovariate".
     */
    public List<String> getCovariateClassNames() {
        return Collections.unmodifiableList(allCovariates.stream().map(cov -> cov.getClass().getSimpleName()).collect(Collectors.toList()));
    }

    /**
     * Returns the size of the list of standard covariates.
     */
    public int size(){
        return allCovariates.size();
    }

    /**
     * Returns a new iterator over all covariates in this list.
     * Note: the list is unmodifiable and the iterator does not support modifying the list.
     */
    @Override
    public Iterator<Covariate> iterator() {
        return allCovariates.iterator();
    }

    public ReadGroupCovariate getReadGroupCovariate() {
        return readGroupCovariate;
    }

    public QualityScoreCovariate getQualityScoreCovariate() {
        return qualityScoreCovariate;
    }

    /**
     * returns an unmodifiable view of the additional covariates stored in this list.
     */
    public Iterable<Covariate> getAdditionalCovariates() {
        return additionalCovariates;
    }

    /**
     * Return a human-readable string representing the used covariates
     *
     * @return a non-null comma-separated string
     */
    public String covariateNames() {
        return String.join(",", getCovariateClassNames());
    }

    /**
     * Get the covariate by the index.
     * @throws IndexOutOfBoundsException if the index is out of range
     *         (<tt>index &lt; 0 || index &gt;= size()</tt>)
     */
    public Covariate get(final int covIndex) {
        return allCovariates.get(covIndex);
    }

    /**
     * Returns the index of the covariate by class name or -1 if not found.
     */
    public int indexByClass(final Class<? extends Covariate> clazz){
        return indexByClass.getOrDefault(clazz, -1);
    }

    /**
     * For each covariate compute the values for all positions in this read and
     * record the values in the provided storage object.
      */
    public void recordAllValuesInStorage(final GATKRead read, final SAMFileHeader header, final ReadCovariates resultsStorage, final boolean recordIndelValues) {
        for (int i = 0, n = allCovariates.size(); i < n; i++) {
            final Covariate cov = allCovariates.get(i);
            resultsStorage.setCovariateIndex(i);
            cov.recordValues(read, header, resultsStorage, recordIndelValues);
        }
    }

    /**
     * Retrieves a covariate by the parsed name {@link Covariate#parseNameForReport()} or null
     * if no covariate with that name exists in the list.
     */
    public Covariate getCovariateByParsedName(final String covName) {
        return allCovariates.stream().filter(cov -> cov.parseNameForReport().equals(covName)).findFirst().orElse(null);
    }

}
