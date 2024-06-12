package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.netflix.servo.util.Objects;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.tuple.Triple;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.BiPredicate;
import java.util.stream.Collectors;

public class PartiallyDeterminedHaplotypeComputationEngineUnitTest extends GATKBaseTest {

    Event SNP_C_90 = new Event("20",90, Allele.REF_A,Allele.ALT_C);
    Event DEL_AAAAAAA_98 = new Event("20",98, Allele.create("AAAAAAA", true),Allele.ALT_A);
    Event SNP_C_100 = new Event("20",100, Allele.REF_A,Allele.ALT_C);
    Event SNP_G_101 = new Event("20",101, Allele.REF_A,Allele.ALT_G);
    Event SNP_G_102 = new Event("20",102, Allele.REF_A,Allele.ALT_G);
    Event SNP_C_104 = new Event("20",104, Allele.REF_A,Allele.ALT_C);
    Event SNP_C_105 = new Event("20",105, Allele.REF_A,Allele.ALT_C);
    Event SNP_G_105 = new Event("20",105, Allele.REF_A,Allele.ALT_G);
    Event SNP_C_106 = new Event("20",106, Allele.REF_A,Allele.ALT_C);
    Event SNP_T_106 = new Event("20",106, Allele.REF_A,Allele.ALT_T);
    Event SNP_C_109 = new Event("20",109, Allele.REF_A,Allele.ALT_C);
    Event SNP_C_107 = new Event("20",107, Allele.REF_A,Allele.ALT_C);

    Event DEL_AA_105 = new Event("20",105, Allele.create("AA", true),Allele.ALT_A);
    Event DEL_AA_100 = new Event("20",100, Allele.create("AA", true),Allele.ALT_A);
    Event DEL_AAA_102 = new Event("20",102, Allele.create("AAA", true),Allele.ALT_A);
    Event DEL_AAAAAAA_102 = new Event("20",102, Allele.create("AAAAAAA", true),Allele.ALT_A);


    Event INS_TT_105 = new Event("20",105, Allele.REF_A, Allele.create("AT"));
    Event INS_TT_103 = new Event("20",103, Allele.REF_A, Allele.create("AT"));
    Event INS_TT_100 = new Event("20",100, Allele.REF_A, Allele.create("AT"));
    Event INS_GGG_106 = new Event("20",106, Allele.REF_A, Allele.create("AGG"));

    // TODO THESE ARE FOR INVALID TEST CASES
    Event SNP_C_99 = new Event("20",99, Allele.REF_A,Allele.ALT_C);
    Event SNP_C_120 = new Event("20",120, Allele.REF_A,Allele.ALT_C);
    Event SNP_G_120 = new Event("20",120, Allele.REF_A,Allele.ALT_G);

    @DataProvider
    public Object[][] makeEventGroupClustersDataProvider() {
        // format:
        // 1) list of all events
        // 2) list of 2- and 3-element groups of events that are mutually excluded due to the Smith-Waterman heuristic
        //      (note that there is no reference sequence in this test, hence we can set whatever  exclusions we want here)
        //      (note also that this is in addition to any mutexes due to overlapping loci)
        //  3) list of sets of events that make up the desired partition into event groups
        // for convenience, 2) and 3) are representing by indices within the list 1).
        return new Object[][] {
                // no mutexes; singleton event groups result
                { List.of(SNP_C_90), List.of(), List.of(List.of(0))},
                { List.of(SNP_C_90, SNP_C_100), List.of(), List.of(List.of(0), List.of(1))},
                { List.of(SNP_C_90, SNP_C_100, SNP_C_105), List.of(), List.of(List.of(0), List.of(1), List.of(2))},
                { List.of(SNP_C_90, SNP_C_100, INS_TT_105, SNP_C_109), List.of(), List.of(List.of(0), List.of(1), List.of(2), List.of(3))},

                // all events are connected by a path of overlaps; everything belongs to a single event group
                { List.of(SNP_C_105, SNP_G_105), List.of(), List.of(List.of(0,1))},
                { List.of(DEL_AAAAAAA_102, SNP_C_105, SNP_G_105), List.of(), List.of(List.of(0,1,2))},
                { List.of(DEL_AAAAAAA_102, SNP_C_105, SNP_G_105, SNP_C_106), List.of(), List.of(List.of(0,1,2,3))},

                // multiple event groups due to independent overlaps -- note that insertions have 0.5 added to their start for DRAGEN
                { List.of(DEL_AAAAAAA_102, SNP_C_105, SNP_G_105, SNP_C_120), List.of(), List.of(List.of(0,1,2), List.of(3))},
                { List.of(SNP_C_105, SNP_G_105, INS_TT_105), List.of(), List.of(List.of(0,1), List.of(2))},
                { List.of(SNP_C_105, SNP_G_105, INS_GGG_106, SNP_C_107), List.of(), List.of(List.of(0,1), List.of(2), List.of(3))},
                { List.of(DEL_AA_100, SNP_G_101, DEL_AA_105, SNP_C_106), List.of(), List.of(List.of(0,1), List.of(2,3))},
                { List.of(DEL_AA_100, SNP_G_101, DEL_AA_105, SNP_C_106, SNP_C_120), List.of(), List.of(List.of(0,1), List.of(2,3), List.of(4))},

                // Smith-Waterman pair mutex joining event groups that would otherwise be independent
                { List.of(SNP_C_90, SNP_C_100), List.of(List.of(0,1)), List.of(List.of(0,1))},
                { List.of(SNP_C_90, SNP_C_100, SNP_C_105), List.of(List.of(0,1)), List.of(List.of(0,1), List.of(2))},
                { List.of(DEL_AA_100, SNP_G_101, DEL_AA_105, SNP_C_106), List.of(List.of(1,2)), List.of(List.of(0,1,2,3))},
                { List.of(DEL_AA_100, SNP_G_101, DEL_AA_105, SNP_C_106, SNP_C_120), List.of(List.of(1,2)), List.of(List.of(0,1, 2,3), List.of(4))},

                // this example is unrealistic. but illustrates that the mutex between the outer SNPs forces the inner SNPinto the same EventGroup
                { List.of(SNP_C_90, SNP_C_100, SNP_C_105), List.of(List.of(0,2)), List.of(List.of(0,1,2))},


                // two Smith-Waterman pair mutexes transitively combining three event groups
                { List.of(DEL_AA_100, SNP_G_101, DEL_AA_105, SNP_C_106, SNP_C_120), List.of(List.of(1,2), List.of(3,4)), List.of(List.of(0,1, 2, 3, 4))},

                // Smith-Waterman pair mutex doing nothing because it is redundant with an overlap mutex
                { List.of(DEL_AA_100, SNP_G_101, DEL_AA_105, SNP_C_106, SNP_C_120), List.of(List.of(2,3)), List.of(List.of(0,1), List.of(2,3), List.of(4))},

                // Smith-Waterman trio mutex transitively combining three event groups
                { List.of(DEL_AA_100, SNP_G_101, DEL_AA_105, SNP_C_106, SNP_C_120), List.of(List.of(1,2,4)), List.of(List.of(0,1, 2, 3, 4))},
        };
    }

    @Test(dataProvider = "makeEventGroupClustersDataProvider")
    public void testMakeEventGroupClusters(List<Event> eventsInOrder, List<List<Integer>> swMutexes, List<List<Integer>> expectedEventGroups) {
        // convert indices to events
        final List<List<Event>> mutexes = swMutexes.stream()
                .map(mutexIndices -> mutexIndices.stream().map(eventsInOrder::get).toList())
                .toList();

        // set of actual partition -- each list in the set is in HAPLOTYPE_SNP_FIRST_COMPARATOR order because that's how the code works
        final Set<List<Event>> actualPartition = PartiallyDeterminedHaplotypeComputationEngine.getEventGroupClusters(eventsInOrder, mutexes).stream()
                .map(eventGroup -> eventGroup.eventsInOrderForTesting())
                .collect(Collectors.toSet());

        // set of expected partition -- each list in the set is in HAPLOTYPE_SNP_FIRST_COMPARATOR order because we sort it explicitly
        final Set<List<Event>> expectedPartition = expectedEventGroups.stream()
                .map(eventsList -> eventsList.stream().map(eventsInOrder::get).sorted(PartiallyDeterminedHaplotypeComputationEngine.HAPLOTYPE_SNP_FIRST_COMPARATOR).toList())
                .collect(Collectors.toSet());

        Assert.assertEquals(expectedPartition, actualPartition);
    }

    @DataProvider
    public Object[][] makeBranchesDataProvider() {
        // format:
        // 1) list of all events
        // 2) list of 2- and 3-element groups of events that are mutually excluded due to the Smith-Waterman heuristic
        //      (note that there is no reference sequence in this test, hence we can set whatever  exclusions we want here)
        //      (note also that this is in addition to any mutexes due to overlapping loci)
        //  3) expected determined branches as list of list of sets -- the outer index is for the event groups.  Each event group
        //  (in order of position) has a list of undetermined branches, each undetermined branch being a set of events
        //  external to the event group (when that event group is the determined one)
        //  4) expected determined event sets as list of list of sets -- the outer index is again for the event groups.  Each event group
        //  has a list of *all* allowed subsets of events from that event group, which are used in the determined part of PD haplotypes.
        // for convenience, 2) and 3) are representing by indices within the list 1).
        return new Object[][] {
                // no mutexes, hence singleton event groups and for each event group (which comprises only one event)
                // the undetermined branch includes every other event
                // the determined branch includes every subset since there are no mutexes
                { List.of(SNP_C_90), List.of(),
                        List.of(List.of(Set.of())),
                        List.of(List.of(Set.of(), Set.of(0)))},
                { List.of(SNP_C_90, SNP_C_100), List.of(),
                        List.of(List.of(Set.of(1)), List.of(Set.of(0))),
                        List.of(List.of(Set.of(), Set.of(0)), List.of(Set.of(), Set.of(1)))},
                { List.of(SNP_C_90, SNP_C_100, SNP_C_105), List.of(),
                        List.of(List.of(Set.of(1,2)), List.of(Set.of(0,2)), List.of(Set.of(0,1))),
                        List.of(List.of(Set.of(), Set.of(0)), List.of(Set.of(), Set.of(1)), List.of(Set.of(), Set.of(2)))},
                { List.of(SNP_C_90, SNP_C_100, INS_TT_105, SNP_C_109), List.of(),
                        List.of(List.of(Set.of(1,2,3)), List.of(Set.of(0,2,3)), List.of(Set.of(0,1,3)), List.of(Set.of(0,1,2))),
                        List.of(List.of(Set.of(), Set.of(0)), List.of(Set.of(), Set.of(1)), List.of(Set.of(), Set.of(2)), List.of(Set.of(), Set.of(3)))},

                // all events are connected by a path of overlaps; everything belongs to a single event group,
                // hence the undetermined event set is empty
                // not all subsets are allowed as determined subsets  due to overlaps
                { List.of(SNP_C_105, SNP_G_105), List.of(),
                        List.of(List.of(Set.of())),
                        List.of(List.of(Set.of(), Set.of(0), Set.of(1)))},
                { List.of(DEL_AAAAAAA_102, SNP_C_105, SNP_G_105), List.of(),
                        List.of(List.of(Set.of())),
                        List.of(List.of(Set.of(), Set.of(0), Set.of(1), Set.of(2)))},
                { List.of(DEL_AAAAAAA_102, SNP_C_105, SNP_G_105, SNP_C_106), List.of(),
                        List.of(List.of(Set.of())),
                        List.of(List.of(Set.of(), Set.of(0), Set.of(1), Set.of(2), Set.of(3), Set.of(1,3), Set.of(2,3)))},

                // first event group is spanning deletion and two spanned SNPs; the only undetermined branch is the other SNP
                // second event group is the other SNP and it admits two undetermined branches from the first event group:
                // 1) the spanning deletion alone, and 2) the SNPs, which though overlapping can coexist as undetermiend events
                { List.of(DEL_AAAAAAA_102, SNP_C_105, SNP_G_105, SNP_C_120), List.of(),
                    List.of(List.of(Set.of(3)), List.of(Set.of(0), Set.of(1,2))),
                    List.of(List.of(Set.of(), Set.of(0), Set.of(1), Set.of(2)), List.of(Set.of(), Set.of(3)))},

                // similar to above but now there are two overlapping SNPs at 120, which coexist for undetermined sets
                { List.of(DEL_AAAAAAA_102, SNP_C_105, SNP_G_105, SNP_C_120, SNP_G_120), List.of(),
                        List.of(List.of(Set.of(3,4)), List.of(Set.of(0), Set.of(1,2))),
                        List.of(List.of(Set.of(), Set.of(0), Set.of(1), Set.of(2)), List.of(Set.of(), Set.of(3), Set.of(4)))},

                //  insertion at 106 and SNP at 107 don't overlap so there are 3 event groups: the SNPs at 105, the insertion, and the SNP
                // at 107.  In each case all events outside the determined event group coexist as undetermined events
                { List.of(SNP_C_105, SNP_G_105, INS_GGG_106, SNP_C_107), List.of(),
                    List.of(List.of(Set.of(2,3)), List.of(Set.of(0,1,3)), List.of(Set.of(0,1,2))),
                    List.of(List.of(Set.of(), Set.of(0), Set.of(1)), List.of(Set.of(), Set.of(2)), List.of(Set.of(), Set.of(3)))},

                // two deletion/spanned SNP pairs and a lone SNP.  Hence three event groups.  Each deletion/spanned SNP pair
                // induces two branches when undetermined.  Thus we have 2x2=4 branches when the SNP at 120 is determined.
                { List.of(DEL_AA_100, SNP_G_101, DEL_AA_105, SNP_C_106, SNP_C_120), List.of(),
                    List.of(List.of(Set.of(2,4), Set.of(3,4)), List.of(Set.of(0,4), Set.of(1,4)), List.of(Set.of(0,2), Set.of(0,3), Set.of(1,2), Set.of(1,3))),
                    List.of(List.of(Set.of(), Set.of(0), Set.of(1)), List.of(Set.of(), Set.of(2), Set.of(3)), List.of(Set.of(), Set.of(4)))},

                // two event groups, the horrid mess and the lone SNP;  When the horrid mess is determined there is a lone undetermined
                // branch from the lone SNP.  When the lone SNP is determined, well, it's just total chaos.
                { List.of(DEL_AAAAAAA_98, DEL_AA_100, SNP_G_101, DEL_AAAAAAA_102, DEL_AA_105, SNP_C_106, SNP_C_120), List.of(),
                        List.of(List.of(Set.of(6)),
                            List.of(Set.of(0,4), Set.of(0,5), Set.of(1,3), Set.of(1,4), Set.of(1,5), Set.of(2,3), Set.of(2,4), Set.of(2,5))),
                        List.of(List.of(Set.of(), Set.of(0), Set.of(1), Set.of(2), Set.of(3), Set.of(4), Set.of(5),
                                Set.of(0,4), Set.of(0,5), Set.of(1,3), Set.of(1,4), Set.of(1,5), Set.of(2,3), Set.of(2,4), Set.of(2,5)),
                                    List.of(Set.of(), Set.of(6)))},
        };
    }

    // TODO: update this test: previously the method computed the entire branch, both determined and undetermined
    // TODO: now the method computes only the event sets for the undetermined part of PD haplotypes
    @Test(dataProvider = "makeBranchesDataProvider")
    public void testMakeBranches(List<Event> eventsInOrder, List<List<Integer>> swMutexes,
                                 final List<List<Set<Integer>>> expectedUndeterminedBranchIndices,
                                 final List<List<Set<Integer>>> expectedDeterminedBranchIndices) {
        // convert indices to events
        final List<List<Event>> mutexes = swMutexes.stream()
                .map(mutexIndices -> mutexIndices.stream().map(eventsInOrder::get).toList())
                .toList();

        final List<PartiallyDeterminedHaplotypeComputationEngine.EventGroup> eventGroups =
                PartiallyDeterminedHaplotypeComputationEngine.getEventGroupClusters(eventsInOrder, mutexes);

        Assert.assertEquals(eventGroups.size(), expectedUndeterminedBranchIndices.size(), "wrong number of event groups");

        for (int determinedEventGroupIndex = 0; determinedEventGroupIndex < eventGroups.size(); determinedEventGroupIndex++) {
            final Set<Set<Event>> actualUndeterminedBranches = PartiallyDeterminedHaplotypeComputationEngine
                    .computeUndeterminedBranches(eventGroups, determinedEventGroupIndex).stream().collect(Collectors.toSet());

            final Set<Set<Event>> expectedUndeterminedBranches = expectedUndeterminedBranchIndices.get(determinedEventGroupIndex).stream()
                    .map(indexSet -> indexSet.stream().map(eventsInOrder::get).collect(Collectors.toSet()))
                    .collect(Collectors.toSet());

            Assert.assertEquals(actualUndeterminedBranches, expectedUndeterminedBranches);

            final Set<Set<Event>> actualDeterminedBranches = eventGroups.get(determinedEventGroupIndex).determinedEventSets()
                    .stream().collect(Collectors.toSet());

            final Set<Set<Event>> expectedDeterminedBranches = expectedDeterminedBranchIndices.get(determinedEventGroupIndex).stream()
                    .map(indexSet -> indexSet.stream().map(eventsInOrder::get).collect(Collectors.toSet()))
                    .collect(Collectors.toSet());

            Assert.assertEquals(actualDeterminedBranches, expectedDeterminedBranches);
        }
    }

    @DataProvider
    public Object[][] testConstructHaplotypeFromVariantsDataProvider() {
        return new Object[][] {
                { Collections.emptyList(), "AAAAAAAAAA", "10M", 0 },
                { List.of(SNP_C_105), "AAAAACAAAA", "5M1X4M", 0 },
                { List.of(SNP_C_100), "CAAAAAAAAA", "1X9M", 0 },
                { List.of(SNP_C_109), "AAAAAAAAAC", "9M1X", 0 },
                { List.of(SNP_C_105, SNP_C_106), "AAAAACCAAA", "5M2X3M", 0 },

                { List.of(DEL_AA_105), "AAAAAAAAA", "6M1D3M", 0 },
                { List.of(DEL_AA_100), "AAAAAAAAA", "1M1D8M", 0 },
                { List.of(DEL_AA_105, SNP_C_109), "AAAAAAAAC", "6M1D2M1X", 0 },
                { List.of(DEL_AA_105, SNP_C_107, SNP_C_109), "AAAAAACAC", "6M1D1X1M1X", 0 },

                { List.of(INS_TT_105),  "AAAAAATAAAA", "6M1I4M", 0 },
                { List.of(INS_GGG_106), "AAAAAAAGGAAA", "7M2I3M", 0 },
                { List.of(DEL_AA_100, INS_GGG_106, SNP_C_109), "AAAAAAGGAAC", "1M1D5M2I2M1X", 0 },

                //this tests that SNPS can be inserted immediately prior to (and following) indels
                { List.of( SNP_C_105, DEL_AA_105 ), "AAAAACAAA", "5M1X1D3M", 1 },
                { List.of( SNP_C_100, DEL_AA_100 ), "CAAAAAAAA", "1X1D8M", 1 },
                { List.of( SNP_C_100, DEL_AA_100, SNP_G_102 ), "CGAAAAAAA", "1X1D1X7M", 1 },
                { List.of( SNP_C_105, INS_TT_105 ), "AAAAACTAAAA", "5M1X1I4M", 1 },
                { List.of( SNP_C_100, INS_TT_100, SNP_G_101 ), "CTGAAAAAAAA", "1X1I1X8M", 1 },
                { List.of( SNP_C_100, INS_TT_100, SNP_G_101, SNP_C_105, DEL_AA_105 ), "CTGAAACAAA", "1X1I1X3M1X1D3M", 2 },

                //testing that the logic around anchor bases isn't resulting in variants being dropped accidentally
                { List.of( SNP_C_104, DEL_AA_105 ), "AAAACAAAA", "4M1X1M1D3M", 0 },
                { List.of( SNP_C_104, INS_TT_105 ), "AAAACATAAAA", "4M1X1M1I4M", 0 },

        };
    }

    @Test(dataProvider = "testConstructHaplotypeFromVariantsDataProvider")
    public void basicConstructHaplotypeFromVariants(List<Event> events, String expectedBases, String expectedCigar, int numberOfCompounds) {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, events, true);
        Assert.assertEquals(result.getBases(), expectedBases.getBytes());
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode(expectedCigar));

        // Assert that the resulting event map matches the input variants:
        EventMap resultEMap = result.getEventMap();
        // NOTE, because of representation in VCF lines, the compound alleles get compressed into a single in the event map, here we assert that this is correct.
        Assert.assertEquals(resultEMap.getNumberOfEvents(), events.size() - numberOfCompounds);
    }

    // encapsulate the determined events, constituent event, and determined span
    private static class PDHaplotypeSummary {
        final Set<Event> determinedEvents;
        final Set<Event> undeterminedEvents;
        final SimpleInterval determinedSpan;

        public PDHaplotypeSummary(Set<Event> determinedEvents, Set<Event> undeterminedEvents, SimpleInterval determinedSpan) {
            this.determinedEvents = determinedEvents;
            this.undeterminedEvents = undeterminedEvents;
            this.determinedSpan = determinedSpan;
        }

        public static PDHaplotypeSummary fromPDHaplotye(final Haplotype haplotype) {
            final PartiallyDeterminedHaplotype pdHaplotype = (PartiallyDeterminedHaplotype) haplotype;
            final Set<Event> determinedEvents = pdHaplotype.getDeterminedEvents();
            final Set<Event> undeterminedEvents = pdHaplotype.getConstituentEvents().stream()
                    .filter(event -> !determinedEvents.contains(event))
                    .collect(Collectors.toSet());
            return new PDHaplotypeSummary(determinedEvents, undeterminedEvents, pdHaplotype.getDeterminedSpan());
        }

        @Override
        public boolean equals( final Object other ) {
            return other instanceof PDHaplotypeSummary
                    && determinedSpan.equals(((PDHaplotypeSummary) other).determinedSpan)
                    && determinedEvents.equals(((PDHaplotypeSummary) other).determinedEvents)
                    && undeterminedEvents.equals(((PDHaplotypeSummary) other).undeterminedEvents);
        }

        @Override
        public int hashCode() {
            return Objects.hash(determinedSpan.hashCode(), determinedEvents.hashCode(), undeterminedEvents.hashCode());
        }
    }

    @DataProvider
    public Object[][] testGeneratePDHaplotyesDataProvider() {
        // FORMAT: byte[] refBases, List<Event> eventsInOrder, List<List<Integer>> swMutexIndices,
        //                                        Set<Triple<List<Integer>, List<Integer>, SimpleInterval>> expectedIndices of PDHaplotypeSummaries
        return new Object[][]{
                {"CTTGAAGCTGAG".getBytes(), List.of(SNP_C_105), List.of(),      // lone SNP
                    Set.of(Triple.of(List.of(), List.of(), new SimpleInterval("20", 105, 105)), // ref haplotype
                            Triple.of(List.of(0), List.of(), new SimpleInterval("20", 105, 105)))},  // SNP haplotype
                {"CTTGAAGCTGAG".getBytes(), List.of(SNP_C_105, SNP_G_105), List.of(),      // two overlapping SNPs
                        Set.of(Triple.of(List.of(), List.of(), new SimpleInterval("20", 105, 105)), // ref haplotype
                                Triple.of(List.of(0), List.of(), new SimpleInterval("20", 105, 105)),   // first SNP haplotype
                                Triple.of(List.of(1), List.of(), new SimpleInterval("20", 105, 105)))},  // second SNP haplotype
                {"CTTGAAAGCTGAG".getBytes(), List.of(SNP_C_105, DEL_AA_105), List.of(),      // overlapping deletion and SNP
                        Set.of(Triple.of(List.of(), List.of(), new SimpleInterval("20", 105, 106)), // ref haplotype
                                Triple.of(List.of(0), List.of(), new SimpleInterval("20", 105, 106)),   // first SNP haplotype
                                Triple.of(List.of(1), List.of(), new SimpleInterval("20", 105, 106)))},  // second SNP haplotype
                {"CTTGAAGACTGAG".getBytes(), List.of(SNP_C_105, SNP_C_107), List.of(),      // two non-overlapping SNPs
                        Set.of(Triple.of(List.of(), List.of(1), new SimpleInterval("20", 105, 105)), // ref haplotype (1st SNP's determined span), 2nd SNP undetermined
                                Triple.of(List.of(0), List.of(1), new SimpleInterval("20", 105, 105)),   // first SNP haplotype, 2nd SNP undetermined
                                Triple.of(List.of(), List.of(0), new SimpleInterval("20", 107, 107)),   // ref haplotype (2nd SNP's determined span), 1st SNP undetermined
                                Triple.of(List.of(1), List.of(0), new SimpleInterval("20", 107, 107)))},  // second SNP haplotype, 1st SNP undetermined
                {"CTTGAAGACTGAG".getBytes(), List.of(SNP_C_105, SNP_G_105, SNP_C_107), List.of(),      // two overlapping SNPs and one other SNP
                        Set.of(Triple.of(List.of(), List.of(2), new SimpleInterval("20", 105, 105)),
                                Triple.of(List.of(0), List.of(2), new SimpleInterval("20", 105, 105)),
                                Triple.of(List.of(1), List.of(2), new SimpleInterval("20", 105, 105)),
                                Triple.of(List.of(), List.of(0,1), new SimpleInterval("20", 107, 107)),
                                Triple.of(List.of(2), List.of(0,1), new SimpleInterval("20", 107, 107)))},
                {"CTAAAAGACTGAG".getBytes(), List.of(DEL_AAA_102, SNP_G_102, SNP_C_104), List.of(),      // deletion spanning two non-overlapping SNPs
                        Set.of(Triple.of(List.of(), List.of(), new SimpleInterval("20", 102, 104)),
                                Triple.of(List.of(0), List.of(), new SimpleInterval("20", 102, 104)),
                                Triple.of(List.of(1), List.of(), new SimpleInterval("20", 102, 104)),
                                Triple.of(List.of(2), List.of(), new SimpleInterval("20", 102, 104)),
                                Triple.of(List.of(1,2), List.of(), new SimpleInterval("20", 102, 104)))},
                {"CTAAAAAAAAAAA".getBytes(), List.of(DEL_AAA_102, SNP_G_102, SNP_C_104, SNP_C_107), List.of(),      // deletion spanning two non-overlapping SNPs; one other SNP; everything shares an event group due to the STR
                        Set.of(Triple.of(List.of(), List.of(), new SimpleInterval("20", 102, 107)),
                                Triple.of(List.of(0), List.of(), new SimpleInterval("20", 102, 107)),
                                Triple.of(List.of(1), List.of(), new SimpleInterval("20", 102, 107)),
                                Triple.of(List.of(2), List.of(), new SimpleInterval("20", 102, 107)),
                                Triple.of(List.of(3), List.of(), new SimpleInterval("20", 102, 107)),
                                Triple.of(List.of(0,3), List.of(), new SimpleInterval("20", 102, 107)),
                                Triple.of(List.of(1,2), List.of(), new SimpleInterval("20", 102, 107)),
                                Triple.of(List.of(1,3), List.of(), new SimpleInterval("20", 102, 107)),
                                Triple.of(List.of(2,3), List.of(), new SimpleInterval("20", 102, 107)),
                                Triple.of(List.of(1,2,3), List.of(), new SimpleInterval("20", 102, 107)))}
        };
    }

    @Test(dataProvider = "testGeneratePDHaplotyesDataProvider")
    public void testGeneratePDHaplotyes(final byte[] refBases, final List<Event> eventsInOrder, List<List<Integer>> swMutexIndices,
                                        final Set<Triple<List<Integer>, List<Integer>, SimpleInterval>> expectedIndices) {
        final Set<PDHaplotypeSummary> expected = expectedIndices.stream()
                .map(triple -> new PDHaplotypeSummary(
                        triple.getLeft().stream().map(eventsInOrder::get).collect(Collectors.toSet()),
                        triple.getMiddle().stream().map(eventsInOrder::get).collect(Collectors.toSet()),
                        triple.getRight()))
                .collect(Collectors.toSet());

        final Haplotype refHaplotype = new Haplotype(refBases, true, 500, TextCigarCodec.decode(refBases.length + "M"));
        refHaplotype.setGenomeLocation(new SimpleInterval("20", 100, 100 + refBases.length));

        final List<List<Event>> swMutexes = swMutexIndices.stream()
                .map(indices -> indices.stream().map(eventsInOrder::get).toList()).toList();


        final Set<Haplotype> pdHaplotypes = PartiallyDeterminedHaplotypeComputationEngine.generatePDHaplotypes(
                refHaplotype, eventsInOrder, swMutexes, false, false, false);

        final Set<PDHaplotypeSummary> actual = pdHaplotypes.stream()
                .map(PDHaplotypeSummary::fromPDHaplotye)
                .collect(Collectors.toSet());

        Assert.assertEquals(actual, expected);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void TestOutOfOrderInputs() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<Event> variants = List.of(SNP_C_105, SNP_G_105);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, variants, true);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void TestSNPsOverlapping() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<Event> events = List.of(SNP_C_109, DEL_AA_100);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, events, true);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void TestVariantNotOverlappingHap() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<Event> events = List.of(SNP_C_90);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, events, true);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void TestVariantIndelPartiallyOverlapping() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        List<Event> events = List.of(DEL_AAAAAAA_98);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, events, true);
    }

    //This is a test asserting that a real edge case that was prone to cause failures in the PDHMM is handled properly when compound variants are taken into account.
    //(62,Rlen=1,[C])->(82,Rlen=1,[C])->(84,Rlen=13,[C])
    @Test
    public void testMessyAlignmentSite() {
        Haplotype ref = new Haplotype("AAGAAAGATGGAGGCCCAGCCAGATCTGGACCCCACAGGCCGTCTCCCCACACAGCCATTCATGTGGTCTACTTCCAGCCATTCATGTGGTCTATTTCCAAGAAAATAGCCCATCCCCCCAAGATAACACCTTCTCAAAAACTTTACAGCTTTGTGTCTACACTGATATTTAGGTATTTTCTTTCTTTTTTTTTTATGATTAACACATCTAATTCAAGAATATCTTGGCAGGATATTCCCCGCTTAGGAAATG".getBytes(), true, 575, TextCigarCodec.decode("253M"));
        ref.setGenomeLocation(new SimpleInterval("20", 24152646, 24152898));

        final Event e1 = new Event("20", 24152708, Allele.REF_T, Allele.ALT_C);
        final Event e2 = new Event("20", 24152728, Allele.REF_T, Allele.ALT_C);
        final Event e3 = new Event("20", 24152729, Allele.create("CATGTGGTCTATT", true), Allele.ALT_C);

        final List<Event> events = List.of(e1, e2, e3);

        Haplotype result = PartiallyDeterminedHaplotypeComputationEngine.constructHaplotypeFromEvents(ref, events, true);
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode("62M1X19M1X1M12D157M"));

        // Assert that the resulting event map matches the input variants:
        EventMap resultEMap = result.getEventMap();
        Assert.assertEquals(resultEMap.getNumberOfEvents(), events.size());
        for (Event e : events) {
            Event actualEvent = resultEMap.get(e.getStart());
            Assert.assertNotNull(actualEvent);
            Assert.assertEquals(actualEvent, e);
        }
    }

    @DataProvider
    public Object[][] testGeneratePDHaplotypeDataProvider() {
        // format: List<Event> events, Set<Event> determinedEvents, int start of determined span, int end of determined span
        // String expectedBases, byte[] expectedAltArray, String expectedCigar
        return new Object[][] {
                {List.of(SNP_C_105, SNP_C_106), Set.of(SNP_C_106), 106,106, "AAAAAACAAA", new byte[]{0,0,0,0,0,17,0,0,0,0}, "6M1X3M"},
                {List.of(SNP_C_105, SNP_C_106), Set.of(), 106,106, "AAAAAAAAAA", new byte[]{0,0,0,0,0,17,0,0,0,0}, "10M"},

                {List.of(INS_TT_103, SNP_C_105, SNP_C_106), Set.of(INS_TT_103), 103,103, "AAAATAAAAAA", new byte[]{0,0,0,0,0,0,17,17,0,0,0}, "4M1I6M"},
                {List.of(INS_TT_103, SNP_C_105, SNP_C_106), Set.of(), 103,103, "AAAAAAAAAA",  new byte[]{0,0,0,0,0,17,17,0,0,0}, "10M"},
                {List.of(INS_TT_103, SNP_C_105, SNP_C_106), Set.of(SNP_C_105), 105,105, "AAAATACAAAA", new byte[]{0,0,0,0,6,0,0,17,0,0,0}, "4M1I1M1X4M"},
                {List.of(INS_TT_103, SNP_C_105, SNP_C_106), Set.of(),  105,105, "AAAATAAAAAA", new byte[]{0,0,0,0,6,0,0,17,0,0,0}, "4M1I6M"},

                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), Set.of(DEL_AAA_102), 102,104, "AAAAAAAA"  , new byte[]{0,0,0,17,17,0,0,0}, "3M2D5M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), Set.of(), 102, 104, "AAAAAAAAAA", new byte[]{0,0,0,0,0,17,17,0,0,0}, "10M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), Set.of(SNP_C_105), 105,105,  "AAAAACAAAA", new byte[]{0,0,0,2,4,0,17,0,0,0}, "5M1X4M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), Set.of(), 105,105,  "AAAAAAAAAA", new byte[]{0,0,0,2,4,0,17,0,0,0}, "10M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), Set.of(SNP_C_106), 106,106,  "AAAAAACAAA", new byte[]{0,0,0,2,4,17,0,0,0,0}, "6M1X3M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106), Set.of(), 106,106,  "AAAAAAAAAA", new byte[]{0,0,0,2,4,17,0,0,0,0}, "10M"},

                // making sure we support "complex alleles" from DRAGEN
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106, INS_GGG_106), Set.of(SNP_C_105), 105,105,  "AAAAACAGGAAA", new byte[]{0,0,0,2,4,0,17,2,4,0,0,0}, "5M1X1M2I3M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106, SNP_T_106, INS_GGG_106), Set.of(), 105, 105, "AAAAAAAGGAAA", new byte[]{0,0,0,2,4,0,81,2,4,0,0,0}, "7M2I3M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106, INS_GGG_106), Set.of(DEL_AAA_102),  102, 104, "AAAAAGGAAA", new byte[]{0,0,0,17,17,2,4,0,0,0}, "3M2D2M2I3M"},
                {List.of(DEL_AAA_102, SNP_C_105, SNP_C_106, SNP_T_106, INS_GGG_106), Set.of(), 102,104,  "AAAAAAAGGAAA", new byte[]{0,0,0,0,0,17,81,2,4,0,0,0}, "7M2I3M"},
                {List.of(SNP_G_101, SNP_C_105, DEL_AA_105), Set.of(SNP_G_101), 101,101,  "AGAAAAAAAA", new byte[]{0,0,0,0,0,17,6,0,0,0}, "1M1X8M"},
                {List.of(SNP_G_101, SNP_C_105, DEL_AA_105), Set.of(), 101,101,   "AAAAAAAAAA", new byte[]{0,0,0,0,0,17,6,0,0,0}, "10M"},
        };
    }
    @Test(dataProvider = "testGeneratePDHaplotypeDataProvider")
    public void testGeneratePDHaplotypeFromVariants(List<Event> events, Set<Event> determinedEvents, int startOfDeterminedSpan, int endOfDeterminedSpan, String expectedBases, byte[] expectedAltArray, String expectedCigar) {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));
        final SimpleInterval determinedSpan = new SimpleInterval("20", startOfDeterminedSpan, endOfDeterminedSpan);
        final List<Event> constituentEvents = events.stream()
                .filter(event -> determinedEvents.contains(event) || !determinedSpan.overlaps(event))
                .toList();  // haplotypes exclude events in the determined span that are not determined
        final PartiallyDeterminedHaplotype result = PartiallyDeterminedHaplotypeComputationEngine.createNewPDHaplotypeFromEvents(ref, determinedEvents, determinedSpan, constituentEvents);
        Assert.assertEquals(new String(result.getBases()), expectedBases);
        Assert.assertEquals(result.getAlternateBases(), expectedAltArray);
        Assert.assertEquals(result.getDeterminedEvents(), determinedEvents);
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode(expectedCigar));
        Assert.assertEquals(result.getDeterminedSpan(), determinedSpan);
    }

    // NOTE: This is an enforcement of a behavior that I consider to be a bug in DRAGEN. Specifically my assumption that we needn't ever concern
    // ourselves with overlapping variants turns out to be false... As it turns out in DRAGEN, they are entirely accepting of constructing a
    // PD haplotype that is REF at bases that underlie a spanning deletion... This means (for example) that if we have a 10 base undetermined
    // deletion from 100-109 and we have a determined ref deletion at position 105-106, that we should STILL construct the haplotype with
    // PD bases from 100-109 even though it means we are assigning that deletion at position 100 to be ref (essentially enforcing that we
    // don't handle spanning deletions). Joint Detection will likely override this behavior in the future.
    @Test
    public void testDeletionUnderlappingDeterminedBases() {
        Haplotype ref = new Haplotype("AAAAAAAAAA".getBytes(), true, 500, TextCigarCodec.decode("10M"));
        ref.setGenomeLocation(new SimpleInterval("20", 100, 110));

        // TODO: fix all this
        /*PartiallyDeterminedHaplotype result = PartiallyDeterminedHaplotypeComputationEngine.createNewPDHaplotypeFromEvents(ref, Set.of(), 105, List.of(DEL_AAAAAAA_102, DEL_AA_105));
        Assert.assertEquals(new String(result.getBases()), "AAAAAAAAAA");
        Assert.assertEquals(result.getAlternateBases(), new byte[]{0,0,0,2,0,0,0,0,4,0});
        Assert.assertEquals(result.getCigar(), TextCigarCodec.decode("10M"));

        // TODO: replace with determined span
        Assert.assertEquals(result.getDeterminedPosition(), DEL_AA_105.getStart());*/
    }

    @Test
    public void testEventsOverlapForPDHapsCode() {
        final BiPredicate<Event, Event> overlaps = PartiallyDeterminedHaplotypeComputationEngine::eventsOverlapForPDHapsCode;

        // easy SNP cases
        Assert.assertFalse(overlaps.test(SNP_C_100, SNP_G_101));
        Assert.assertFalse(overlaps.test(SNP_C_100, SNP_G_102));
        Assert.assertFalse(overlaps.test(SNP_C_107, SNP_G_105));
        Assert.assertTrue(overlaps.test(SNP_C_105, SNP_G_105));
        Assert.assertTrue(overlaps.test(SNP_T_106, SNP_T_106));

        // overlap of SNP and deletion -- note that we add 1 to deletion start but not to deletion end
        Assert.assertFalse(overlaps.test(DEL_AAA_102, SNP_G_101));
        Assert.assertFalse(overlaps.test(DEL_AAA_102, SNP_G_102));
        Assert.assertTrue(overlaps.test(DEL_AAA_102, SNP_C_104));
        Assert.assertFalse(overlaps.test(DEL_AAA_102, SNP_C_105));

        // overlap of SNP and insertion -- note that we add 0.5 to insertion start and end
        Assert.assertFalse(overlaps.test(SNP_G_102, INS_TT_103));
        Assert.assertFalse(overlaps.test(SNP_C_104, INS_TT_103));
        Assert.assertFalse(overlaps.test(SNP_C_105, INS_TT_105));

        // two insertions should overlap only if they occur at the same position
        Assert.assertTrue(overlaps.test(INS_TT_105, INS_TT_105));
        Assert.assertFalse(overlaps.test(INS_TT_105, INS_GGG_106));

        // two deletions
        Assert.assertTrue(overlaps.test(DEL_AAAAAAA_102, DEL_AAA_102));
        Assert.assertTrue(overlaps.test(DEL_AA_105, DEL_AAAAAAA_102));
        Assert.assertFalse(overlaps.test(DEL_AA_100, DEL_AAA_102));

        // deletion and insertion
        Assert.assertFalse(overlaps.test(INS_TT_105, DEL_AA_105));  // add 1 to deletion start but only 0.5 to insertion end
        Assert.assertFalse(overlaps.test(INS_TT_103, DEL_AA_105));
        Assert.assertTrue(overlaps.test(DEL_AAAAAAA_102, INS_GGG_106));
        Assert.assertTrue(overlaps.test(INS_TT_103, DEL_AAA_102));
    }

    @Test
    public void testFindSTRs() {
        final String contig = "CONTIG";
        final int refStart = 100;
        // 4 units of AC (8 bases), 9 bases of nothing, 10-base T homopolymer, 6 bases of nothing
        final String refBases = "ACACACAC" + "GTAGGTACG" + "TTTTTTTTTT" + "ACGCTG";
        final Haplotype refHaplotype = new Haplotype(refBases.getBytes(), new SimpleInterval(contig, refStart, refStart + refBases.length()));

        final List<SimpleInterval> strs = PartiallyDeterminedHaplotypeComputationEngine.findSTRs(refHaplotype);

        final List<SimpleInterval> expected = List.of(new SimpleInterval(contig, refStart + 0, refStart + 8),
                new SimpleInterval(contig, refStart + 8 + 9, refStart + 8 + 9 + 10));
        Assert.assertEquals(strs, expected);
    }
}