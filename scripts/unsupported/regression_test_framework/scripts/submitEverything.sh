#!/usr/bin/env bash

cromshell submit ./wdl/compareHaplotypeCallerRunsNotCached.wdl ./json/HaplotypeCallerImporvemnt1.json /Users/emeryj/hellbender/Scripts/markDuplicatesTesting/picardComparison/stressTest/withSplittingIndexTrial/PairedSingleSampleWf.options ./wdl/wdlinputs.zip

sleep 1

cromshell submit ./wdl/compareHaplotypeCallerRunsNotCached.wdl ./json/HaplotypeCallerImporvemnt1_chr1.json /Users/emeryj/hellbender/Scripts/markDuplicatesTesting/picardComparison/stressTest/withSplittingIndexTrial/PairedSingleSampleWf.options ./wdl/wdlinputs.zip

sleep 1

cromshell submit ./wdl/compareHaplotypeCallerRunsNotCached.wdl ./json/HaplotypeCallerImporvemnt2.json /Users/emeryj/hellbender/Scripts/markDuplicatesTesting/picardComparison/stressTest/withSplittingIndexTrial/PairedSingleSampleWf.options ./wdl/wdlinputs.zip

sleep 1

cromshell submit ./wdl/compareHaplotypeCallerRunsNotCached.wdl ./json/HaplotypeCallerImporvemnt2_chr1.json /Users/emeryj/hellbender/Scripts/markDuplicatesTesting/picardComparison/stressTest/withSplittingIndexTrial/PairedSingleSampleWf.options ./wdl/wdlinputs.zip

sleep 1

cromshell submit ./wdl/compareHaplotypeCallerRunsNotCached.wdl ./json/HaplotypeCallerImporvemnt3.json /Users/emeryj/hellbender/Scripts/markDuplicatesTesting/picardComparison/stressTest/withSplittingIndexTrial/PairedSingleSampleWf.options ./wdl/wdlinputs.zip

sleep 1

cromshell submit ./wdl/compareHaplotypeCallerRunsNotCached.wdl ./json/HaplotypeCallerImporvemnt3_chr1.json /Users/emeryj/hellbender/Scripts/markDuplicatesTesting/picardComparison/stressTest/withSplittingIndexTrial/PairedSingleSampleWf.options ./wdl/wdlinputs.zip

sleep 1

cromshell submit ./wdl/compareHaplotypeCallerRunsNotCached.wdl ./json/HaplotypeCallerImporvemntAll.json /Users/emeryj/hellbender/Scripts/markDuplicatesTesting/picardComparison/stressTest/withSplittingIndexTrial/PairedSingleSampleWf.options ./wdl/wdlinputs.zip

sleep 1

cromshell submit ./wdl/compareHaplotypeCallerRunsNotCached.wdl ./json/HaplotypeCallerImporvemntAll_chr1.json /Users/emeryj/hellbender/Scripts/markDuplicatesTesting/picardComparison/stressTest/withSplittingIndexTrial/PairedSingleSampleWf.options ./wdl/wdlinputs.zip