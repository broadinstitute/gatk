package org.broadinstitute.hellbender.tools.exome.sexgenotyper;

import com.google.common.collect.Sets;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.Target;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import javax.annotation.Nonnull;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * A collection of {@link Target} instances along with helper methods for generating their genotype ploidy
 * annotations on-the-fly using provided contig ploidy annotations.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlinePloidyAnnotatedTargetCollection implements TargetCollection<Target> {

    /**
     * Map from targets to their germline ploidy annotations (based on target contigs)
     */
    private Map<Target, ContigGermlinePloidyAnnotation> targetToContigPloidyAnnotationMap;

    /**
     * List of autosomal targets
     */
    private List<Target> autosomalTargetList;

    /**
     * List of allosomal targets
     */
    private List<Target> allosomalTargetList;

    /**
     * List of all targets
     */
    private List<Target> fullTargetList;

    /**
     * Set of all targets
     */
    private Set<Target> fullTargetSet;


    /**
     * Public constructor.
     *
     * @param contigAnnotsList list of contig germline ploidy annotations.
     * @param targetList list of targets
     */
    public GermlinePloidyAnnotatedTargetCollection(@Nonnull final List<ContigGermlinePloidyAnnotation> contigAnnotsList,
                                                   @Nonnull final List<Target> targetList) {
        performValidityChecks(targetList, contigAnnotsList);

        fullTargetList = Collections.unmodifiableList(targetList);
        fullTargetSet = Collections.unmodifiableSet(new HashSet<>(fullTargetList));

        /* map targets to ploidy annotations */
        final Map<String, ContigGermlinePloidyAnnotation> contigNameToContigPloidyAnnotationMap = contigAnnotsList.stream()
                .collect(Collectors.toMap(ContigGermlinePloidyAnnotation::getContigName, Function.identity()));
        targetToContigPloidyAnnotationMap = Collections.unmodifiableMap(
                fullTargetList.stream().collect(Collectors.toMap(Function.identity(),
                        target -> contigNameToContigPloidyAnnotationMap.get(target.getContig()))));

        /* autosomal and allosomal target lists */
        autosomalTargetList = Collections.unmodifiableList(fullTargetList.stream()
                .filter(target -> targetToContigPloidyAnnotationMap.get(target).getContigClass() == ContigClass.AUTOSOMAL)
                .collect(Collectors.toList()));

        allosomalTargetList = Collections.unmodifiableList(fullTargetList.stream()
                .filter(target -> targetToContigPloidyAnnotationMap.get(target).getContigClass() == ContigClass.ALLOSOMAL)
                .collect(Collectors.toList()));
    }

    /**
     * Returns an unmodifiable list of autosomal targets contained in the collection
     * @return unmodifiable list of contained autosomal targets
     */
    public List<Target> getAutosomalTargetList() {
        return Collections.unmodifiableList(autosomalTargetList);
    }

    /**
     * Returns an unmodifiable list of allosomal targets contained in the collection
     * @return unmodifiable list of contained allosomal targets
     */
    public List<Target> getAllosomalTargetList() {
        return Collections.unmodifiableList(allosomalTargetList);
    }

    /**
     * Returns the ploidy of a target for a given ploidy tag (= genotype identifier string)
     *
     * @param target target in question
     * @param genotypeName genotype identifier string
     * @return integer ploidy
     */
    public int getTargetGermlinePloidyByGenotype(@Nonnull final Target target, @Nonnull final String genotypeName) {
        if (!fullTargetSet.contains(target)) {
            throw new IllegalArgumentException("Target \"" + target.getName() + "\" can not be found");
        }
        return targetToContigPloidyAnnotationMap.get(target).getGermlinePloidy(genotypeName);
    }

    /**
     * Perform a number of checks on the arguments passed to the constructor:
     * <dl>
     *     <dt> Assert both lists are non-empty </dt>
     *     <dt> Assert targets have unique names </dt>
     *     <dt> Assert each contig is annotated only once </dt>
     *     <dt> Assert all targets have annotated contigs </dt>
     * </dl>
     *
     * @param targetList list of targets
     * @param contigAnnotsList list of contig ploidy annotations
     */
    private void performValidityChecks(final List<Target> targetList, final List<ContigGermlinePloidyAnnotation> contigAnnotsList) {
        /* assert the lists are non-empty */
        if (targetList.isEmpty()) {
            throw new UserException.BadInput("Target list can not be empty");
        }
        if (contigAnnotsList.isEmpty()) {
            throw new UserException.BadInput("Contig germline ploidy annotation list can not be empty");
        }

        /* assert targets have unique names */
        final Map<String, Long> targetNameCounts = targetList.stream()
                .map(Target::getName)
                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));
        if (targetNameCounts.keySet().size() < targetList.size()) {
            throw new UserException.BadInput("Targets must have unique names. Non-unique target names: " +
                    targetNameCounts.keySet().stream()
                            .filter(name -> targetNameCounts.get(name) > 1)
                            .collect(Collectors.joining(", ")));
        }

        /* assert contigs are not annotated multiple times */
        final Map<String, Long> contigAnnotsCounts = contigAnnotsList.stream()
                .map(ContigGermlinePloidyAnnotation::getContigName)
                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));
        if (contigAnnotsCounts.keySet().size() < contigAnnotsList.size()) {
            throw new UserException.BadInput("Some contigs are multiply annotated: " +
                    contigAnnotsCounts.keySet().stream()
                            .filter(contig -> contigAnnotsCounts.get(contig) > 1) /* multiply annotated contigs */
                            .collect(Collectors.joining(", ")));
        }

        /* assert all contigs present in the target list are annotated */
        final Set<String> contigNamesFromTargets = targetList.stream()
                .map(Target::getContig).collect(Collectors.toSet());
        final Set<String> contigNamesFromAnnots = contigAnnotsList.stream()
                .map(ContigGermlinePloidyAnnotation::getContigName).collect(Collectors.toSet());
        final Set<String> missingContigs = Sets.difference(contigNamesFromTargets, contigNamesFromAnnots);
        if (missingContigs.size() > 0) {
            throw new UserException.BadInput("All contigs must be annotated. Annotations are missing for: " +
                    missingContigs.stream().collect(Collectors.joining(", ")));
        }

        /* assert all contigs have annotations for all ploidy classes */
        final Set<String> firstAnnotPloidyTagSet = contigAnnotsList.get(0).getGenotypesSet();
        if (contigAnnotsList.stream().filter(annot -> !annot.getGenotypesSet().equals(firstAnnotPloidyTagSet)).count() > 0) {
            throw new UserException.BadInput("Not all entries in the contig germline ploidy annotation list have the same " +
                    "set of genotypes");
        }
    }

    @Override
    public int targetCount() {
        return fullTargetList.size();
    }

    @Override
    public Target target(int index) {
        return fullTargetList.get(index);
    }

    @Override
    public String name(Target target) {
        return target.getName();
    }

    @Override
    public SimpleInterval location(int index) {
        return location(target(index));
    }

    @Override
    public SimpleInterval location(Target target) {
        return new SimpleInterval(target.getContig(), target.getStart(), target.getEnd());
    }

    @Override
    public List<Target> targets() {
        return Collections.unmodifiableList(fullTargetList);
    }


    /**
     * Not implemented yet (currently is not required in the use case of this class)
     */
    @Override
    public int index(String name) {
        throw new UnsupportedOperationException();
    }

    /**
     * Not implemented yet (currently is not required in the use case of this class)
     */
    @Override
    public IndexRange indexRange(Locatable location) {
        throw new UnsupportedOperationException();
    }

}
