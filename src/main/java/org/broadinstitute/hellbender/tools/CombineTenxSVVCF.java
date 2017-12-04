package org.broadinstitute.hellbender.tools;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.*;

@CommandLineProgramProperties(summary = "Combine multiple 10x VCF files (normalized by NormalizeTenxVCFs",
        oneLineSummary = "Combine multiple 10x VCF files (normalized by NormalizeTenxVCFs",
        programGroup = VariantProgramGroup.class)
public class CombineTenxSVVCF extends MultiVariantWalker {

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The combined VCF output file", optional=false)
    private File outputFile;

    @Argument(fullName="breakpointMergeThreshold", shortName="breakpointMergeThreshold", doc = "Distance at which to merge coherent breakpoints")
    protected int breakpointMergeThreshold = 500;

    LinkedList<VariantContext> currentBreakpoints = new LinkedList<>();

    LinkedList<SVClique> currentCliques = new LinkedList<>();
    Map<VariantContext, SVClique> cliqueAssignments = new HashMap<>();

    Map<String, String> idMappings = new HashMap<>();

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        final String contig = variant.getContig();
        final int loc = variant.getStart();

        // if any cliques are out of scope, emit them by making an exemplar, reassigning the other variants to be genotypes of the exemplar,
        // removing the variants from current breakpoints, clique assignments, and currentCliques
        for (final Iterator<SVClique> iterator = currentCliques.iterator(); iterator.hasNext(); ) {
            final SVClique currentClique = iterator.next();
            if (currentClique.outOfScope(breakpointMergeThreshold, variant)) {
                System.err.print("emitting clique: " + currentClique);
            }
        }
        
        // go back through the list. if variant matches previous
        Iterator<SVClique> svCliqueIterator = currentCliques.descendingIterator();
        while (svCliqueIterator.hasNext()) {
            SVClique prevClique =  svCliqueIterator.next();
            
        }

        // if this vc does match a current


    }

    static boolean concordant(final int breakpointMergeThreshold, final VariantContext vc1, final VariantContext vc2) {
        if (vc1.getAttribute("SVTYPE").equals("BND") && vc2.getAttribute("SVTYPE").equals("BND")) {
            // todo
            Utils.validate(vc1.getAlternateAlleles().size() == 1 && vc2.getAlternateAlleles().size() == 1, "can't handle BND variants with more than one alt");

            final String vc1Contig = vc1.getContig();
            final int vc1Pos = vc1.getStart();
            final BreakendAdjacency bnd1 = parseBreakendAllele(vc1.getAlternateAllele(0).getDisplayString());

            final String vc2Contig = vc2.getContig();
            final int vc2Pos = vc2.getStart();
            final BreakendAdjacency bnd2 = parseBreakendAllele(vc2.getAlternateAllele(0).getDisplayString());

            // todo: local bases?
            return (vc1Contig.equals(vc2Contig) &&
                    Math.abs(vc1Pos - vc2Pos) <= breakpointMergeThreshold &&
                    bnd1.before == bnd2.before &&
                    bnd1.revComp == bnd2.revComp &&
                    bnd1.contig.equals(bnd2.contig) &&
                    Math.abs(bnd1.position - bnd2.position) <= breakpointMergeThreshold);

        } else {
            return false;
        }
    }


    static BreakendAdjacency parseBreakendAllele(final String breakendAllele) {
        final String localBases;
        final boolean before;
        final boolean revComp;
        final SimpleInterval adjacentPos;

        final String[] fields = breakendAllele.split("(?<=(\\[|\\])|(?=(\\[|\\])))");
        if (fields[0].equals("[") || fields[0].equals("]")) {
            before = true;
            if (fields[0].equals("]")) {
                revComp = false;
            } else {
                revComp = true;
            }
            adjacentPos = new SimpleInterval(fields[1]);
            localBases = fields[fields.length - 1];
        } else if (fields[fields.length -1].equals("[") || fields[fields.length -1].equals("]")) {
            before = false;
            if (fields[fields.length - 1].equals("[")) {
                revComp = false;
            } else {
                revComp = true;
            }
            adjacentPos = new SimpleInterval(fields[fields.length - 2]);
            localBases = fields[0];
        } else {
            // todo: unpaired breakends
            throw new GATKException("Don't know how to parse BND allele " + breakendAllele);
        }
        return new BreakendAdjacency(localBases, adjacentPos.getContig(), adjacentPos.getStart(), before, revComp);
    }

    static class BreakendAdjacency {
        final String localBases;
        final String contig;
        final int position;
        final boolean before;
        final boolean revComp;

        public BreakendAdjacency(final String localBases, final String contig, final int position, final boolean before, final boolean revComp) {
            this.localBases = localBases;
            this.contig = contig;
            this.position = position;
            this.before = before;
            this.revComp = revComp;
        }
    }

    static class SVClique {

        final Set<VariantContext> members = new HashSet<>();
        VariantContext max= null;

        public void add(VariantContext vc) {
            if (max == null || vc.getStart() > max.getStart()) {
                max = vc;
            }
            members.add(vc);
        }

        public boolean outOfScope(final int breakpointMergeThreshold, final VariantContext other) {
            return ! concordant(breakpointMergeThreshold, other, max);
        }
        
        public SVClique shiftToInclude(final VariantContext other, final int breakpointMergeThreshold) {
            final SVClique newClique = new SVClique();
            for (final VariantContext member : members) {
                if (concordant(breakpointMergeThreshold, member, other)) {
                    newClique.add(member);
                }
            }
            newClique.add(other);
            return newClique;
        }
    }
}
