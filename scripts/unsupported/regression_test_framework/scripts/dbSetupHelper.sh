#!/usr/bin/env bash

# Path to local files associated with the runs:
LOCAL_ROOT=~/.cromshell/cromwell-v36.dsde-methods.broadinstitute.org

# List of runs:
HC_INITIAL_RUN_LIST="84edd0ed-f084-47bb-af11-d8b47b9f1865 683dad15-3dea-4f35-8826-6d31f0e0c7bc 334ef0b9-ce97-49b3-8728-b6400396cde7 be5c2094-9ed2-4898-903d-cf519128ca48 5e7bc348-4745-4e6b-8231-bae0f57fc0b0 73a9ee75-6006-40f4-8f1b-681c38a501a8 662f5bfb-038c-491d-925b-896cc1038ff2 301fbc8e-2be1-41bd-847f-0dc4aff9f9af 001d8293-6c8c-41b9-8f07-274e49039d9b 1b557c70-a9e8-40f7-a4b2-d33c792c0255 88aba234-c27b-4bad-842e-9662704d64ca"

RUN_LIST=${HC_INITIAL_RUN_LIST}

################################################################################
################################################################################

# Turn on error exiting:
set -e

################################################################################
################################################################################

function printVar() 
{
  echo "$1 = ${!1}"
}

function getFileType()
{
  FT="UNSTRUCTURED_TEXT"
  bn=$( basename $1 )
  if [[ $bn =~ ^.*.bam$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="GENOME_READS"
  elif [[ $bn =~ ^.*.bam$ ]] ; then
    FT="EXOME_READS"
  elif [[ $bn =~ ^.*.bai$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="GENOME_READS_INDEX"
  elif [[ $bn =~ ^.*.bai$ ]] ; then
    FT="EXOME_READS_INDEX"
  elif [[ $bn =~ ^.*.vcf.gz$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="SOMATIC_VCF"
  elif [[ $bn =~ ^.*.vcf.gz$ ]] ; then
    FT="SOMATIC_VCF"
  elif [[ $bn =~ ^.*.vcf.gz.tbi$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="SOMATIC_VCF_INDEX"
  elif [[ $bn =~ ^.*.vcf.gz.tbi$ ]] ; then
    FT="SOMATIC_VCF_INDEX"
  elif [[ $bn =~ ^.*.vcf$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="SOMATIC_VCF"
  elif [[ $bn =~ ^.*.vcf$ ]] ; then
    FT="SOMATIC_VCF"
  elif [[ $bn =~ ^.*.vcf.idx$ ]] && [[ $bn =~ ^G.* ]] ; then
    FT="SOMATIC_VCF_INDEX"
  elif [[ $bn =~ ^.*.vcf.idx$ ]] ; then
    FT="SOMATIC_VCF_INDEX"
  fi

  echo "${FT}"
}

function getOutputFileID() 
{
  local fileBaseName=$1
  local runUUID=$2

  outputFileIndex=0
  if [[ "${runUUID}" == "84edd0ed-f084-47bb-af11-d8b47b9f1865" ]] ; then
    outputFileIndex=1
  elif [[ "${runUUID}" == "683dad15-3dea-4f35-8826-6d31f0e0c7bc" ]] ; then
    outputFileIndex=19
  elif [[ "${runUUID}" == "334ef0b9-ce97-49b3-8728-b6400396cde7" ]] ; then
    outputFileIndex=37
  elif [[ "${runUUID}" == "be5c2094-9ed2-4898-903d-cf519128ca48" ]] ; then
    outputFileIndex=55
  elif [[ "${runUUID}" == "5e7bc348-4745-4e6b-8231-bae0f57fc0b0" ]] ; then
    outputFileIndex=73
  elif [[ "${runUUID}" == "73a9ee75-6006-40f4-8f1b-681c38a501a8" ]] ; then
    outputFileIndex=91
  elif [[ "${runUUID}" == "662f5bfb-038c-491d-925b-896cc1038ff2" ]] ; then
    outputFileIndex=109
  elif [[ "${runUUID}" == "301fbc8e-2be1-41bd-847f-0dc4aff9f9af" ]] ; then
    outputFileIndex=127
  elif [[ "${runUUID}" == "001d8293-6c8c-41b9-8f07-274e49039d9b" ]] ; then
    outputFileIndex=145
  elif [[ "${runUUID}" == "1b557c70-a9e8-40f7-a4b2-d33c792c0255" ]] ; then
    outputFileIndex=163
  elif [[ "${runUUID}" == "88aba234-c27b-4bad-842e-9662704d64ca" ]] ; then
    outputFileIndex=181
  fi

  if [[ "${fileBaseName}" == "NexPond-359781" ]] ; then
    let outputFileIndex=$outputFileIndex+0
  elif [[ "${fileBaseName}" == "NexPond-359877" ]] ; then
    let outputFileIndex=$outputFileIndex+2
  elif [[ "${fileBaseName}" == "NexPond-360361" ]] ; then
    let outputFileIndex=$outputFileIndex+4
  elif [[ "${fileBaseName}" == "NexPond-360457" ]] ; then
    let outputFileIndex=$outputFileIndex+6
  elif [[ "${fileBaseName}" == "NexPond-361337" ]] ; then
    let outputFileIndex=$outputFileIndex+8
  elif [[ "${fileBaseName}" == "NexPond-361433" ]] ; then
    let outputFileIndex=$outputFileIndex+10
  elif [[ "${fileBaseName}" == "NexPond-362428" ]] ; then
    let outputFileIndex=$outputFileIndex+12
  elif [[ "${fileBaseName}" == "NexPond-363907" ]] ; then
    let outputFileIndex=$outputFileIndex+14
  elif [[ "${fileBaseName}" == "NexPond-445394" ]] ; then
    let outputFileIndex=$outputFileIndex+16
  fi

  echo $outputFileIndex
}

function createConcordanceJson() 
{
  echo "{"
  echo "  \"Concordance.gatk_docker\": \"broadinstitute/gatk-nightly:2019-02-26-4.1.0.0-31-g23bd0a2f8-SNAPSHOT\","
  echo ""
  echo "  \"Concordance.eval_vcf\": \"${1}\","
  echo "  \"Concordance.eval_vcf_idx\": \"${1}.idx\","
  echo "  \"Concordance.truth_vcf\": \"${2}\","
  echo "  \"Concordance.truth_vcf_idx\": \"${2}.tbi\","
  echo "  \"Concordance.intervals\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\","
  echo ""
  echo "  \"Concordance.default_disk_space_gb\": 512,"
  echo "  \"Concordance.mem_gb\": 32,"
  echo "  \"Concordance.boot_disk_size_gb\": 64"
  echo "}"
}

function createGenotypeConcordanceJson() 
{
  echo "{"
  echo "  \"GenotypeConcordanceTask.gatk_docker\": \"broadinstitute/gatk-nightly:2019-02-26-4.1.0.0-31-g23bd0a2f8-SNAPSHOT\","
  echo ""
  echo "  \"GenotypeConcordanceTask.call_vcf\": \"${1}\","
  echo "  \"GenotypeConcordanceTask.call_index\": \"${1}.idx\","
  echo "  \"GenotypeConcordanceTask.call_sample\": \"${1}.idx\","
  echo ""
  echo "  \"GenotypeConcordanceTask.truth_vcf\": \"${2}\","
  echo "  \"GenotypeConcordanceTask.truth_index\": \"${2}.tbi\","
  echo "  \"GenotypeConcordanceTask.truth_sample\": \"${2}.tbi\","
  echo ""
  echo "  \"GenotypeConcordanceTask.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\","
  echo ""
  echo "  \"GenotypeConcordanceTask.output_base_name\": \"${3}\","
  echo ""
  echo "  \"GenotypeConcordanceTask.default_disk_space_gb\": 512,"
  echo "  \"GenotypeConcordanceTask.mem_gb\": 32,"
  echo "  \"GenotypeConcordanceTask.boot_disk_size_gb\": 64"
  echo "}"
}

# ============================================================================== 
# DB Functions:

MYSQL_CMD_TEMPLATE='mysql -N -s -h mysql-prd2.broadinstitute.org -u DSPRegressionTestViewer --password=DSPRegressionTestViewer -e'

function getAnalysisRunIDFromDB() 
{
  ${MYSQL_CMD_TEMPLATE} "select max(idAnalysisRun) from DSPRegressionTesting.AnalysisRuns;" 2>/dev/null
}

function getNextOutputFileIDFromDB() 
{
  ${MYSQL_CMD_TEMPLATE} "select max(idOutputFiles)+1 from DSPRegressionTesting.OutputFiles;" 2>/dev/null
}

function getNextMetricIDFromDB() 
{
  ${MYSQL_CMD_TEMPLATE} "select max(idMetrics)+1 from DSPRegressionTesting.Metrics;" 2>/dev/null
}

function getNextMetricTimingIDFromDB() 
{
  ${MYSQL_CMD_TEMPLATE} "select max(idMetric_Timing)+1 from DSPRegressionTesting.Metric_Timing;" 2>/dev/null
}

function getNextMetricCFAIDFromDB() 
{
  ${MYSQL_CMD_TEMPLATE} "select max(idConcordance_FilterAnalysis)+1 from DSPRegressionTesting.Metric_Concordance_FilterAnalysis;" 2>/dev/null
}

function getNextMetricCSIDFromDB() 
{
  ${MYSQL_CMD_TEMPLATE} "select max(idConcordance_Summary)+1 from DSPRegressionTesting.Metric_Concordance_Summary;" 2>/dev/null
}

function getNextMetricGCCIDFromDB() 
{
  ${MYSQL_CMD_TEMPLATE} "select max(id)+1 from DSPRegressionTesting.Metric_GenotypeConcordance_ContingencyMetrics;" 2>/dev/null
}

function getNextMetricGCSIDFromDB() 
{
  ${MYSQL_CMD_TEMPLATE} "select max(idGenotypeConcordance_SummaryMetrics)+1 from DSPRegressionTesting.Metric_GenotypeConcordance_SummaryMetrics;" 2>/dev/null
}

function getNextMetricGCDIDFromDB() 
{
  ${MYSQL_CMD_TEMPLATE} "select max(idGenotypeConcordance_DetailMetrics)+1 from DSPRegressionTesting.Metric_GenotypeConcordance_DetailMetrics;" 2>/dev/null
}

################################################################################
################################################################################

doToolDataCreation=false
doAnalysisDataCreation=true

################################################################################
################################################################################

if $doToolDataCreation ; then 

  echo "Creating DATA file outputs."

  # Setup output files in the order they should be inserted:
  outputFilesQueryFile=D01_outputFilesQueryFile.sql
  testScenarioOutputFilesQueryFile=D02_testScenarioOutputFilesQueryFile.sql

  echo "INSERT INTO DSPRegressionTesting.OutputFiles (fileType, path, timeCreated, md5sum, sourceType) VALUES " > ${outputFilesQueryFile}
  echo "INSERT INTO DSPRegressionTesting.TestScenarioOutputFiles (scenarioRun, outputFile) VALUES " > ${testScenarioOutputFilesQueryFile}

  scenarioRunIdBase=1
  scenarioId=1
  outputFileId=1
  isFirstFile=true
  for r in ${RUN_LIST} ; do 

    echo "Processing run: ${r}"

    gatkDocker=$( cat ${r}/*SNAPSHOT*.json | jq '."ToolComparisonWdl.gatk_docker"'| tr -d '"' )

    tmpFile=$(mktemp)
    echo -n "  Getting outputs...  "
    ~/Development/cromshell/cromshell list-outputs $r 2>/dev/null | grep 'HaplotypeCaller' | grep -v 'timingInformation' | xargs sutil ls -l | grep -v '^[ \t]*$'  | grep -v '^TOTAL' > $tmpFile
    echo 'DONE!'

    while read outFileInfo ; do

      echo "  Processing output file ${outputFileId}"

      outFile=$( echo $outFileInfo | awk '{print $3}' )

      shardIdx=$( echo $outFile | grep -o shard-[0-9]* | sed 's#shard-##g' )
      let scenarioRunId=$scenarioRunIdBase+$shardIdx

      timeStamp=$( echo $outFileInfo | awk '{print $2}' | tr -d 'Z' | tr 'T' ' ')

      FT=$( getFileType $outFile )

      md5sum=$( gsutil hash -hm $outFile 2>/dev/null | grep md5 | awk '{print $NF}'| \grep -o '^[0-9A-Za-z]*$' )

      if ! $isFirstFile ; then
        echo -n ',' >> ${outputFilesQueryFile}
        echo -n ',' >> ${testScenarioOutputFilesQueryFile}
      fi

      echo "('$FT', '$outFile', '$timeStamp', '$md5sum', 'TOOL')" >> ${outputFilesQueryFile}
      echo "(${scenarioRunId}, ${outputFileId})" >> ${testScenarioOutputFilesQueryFile}
      let outputFileId=${outputFileId}+1
      isFirstFile=false
    done < $tmpFile

    let scenarioId=$scenarioId+1
    let scenarioRunIdBase=$scenarioRunIdBase+9

  done

  echo ";" >> ${outputFilesQueryFile}
  echo ';' >> ${testScenarioOutputFilesQueryFile}

  echo "Done processing TOOL output files."

fi

################################################################################

if $doAnalysisDataCreation ; then

  echo "Creating ANALYSIS file outputs."

  # Setup output files in the order they should be inserted:
  analysisRunsQueryFile=A01_analysisRunsQueryFile.sql
  outputFilesAnalysisQueryFile=A02_outputFilesAnalysisQueryFile.sql
  analysisOutputFilesQueryFile=A03_analysisOutputFilesQueryFile.sql
   
  metricsConcordanceFAQueryFile=A04_metricsConcordanceFAQueryFile.sql
  metricsConcordanceSummaryQueryFile=A05_metricsConcordanceSummaryQueryFile.sql

  metricsGCContingencyQueryFile=A06_metricsGCContingencyQueryFile.sql
  metricsGCDetailQueryFile=A07_metricsGCDetailQueryFile.sql
  metricsGCSummaryQueryFile=A08_metricsGCSummaryQueryFile.sql

  metricsTimingQueryFile=A09_metricsTiming.sql
  
  metricsQueryFile=A10_metrics.sql

  echo "INSERT INTO DSPRegressionTesting.AnalysisRuns (analysisInfoID, scenarioOutputForComparison, configuration) VALUES " > ${analysisRunsQueryFile}
  echo "INSERT INTO DSPRegressionTesting.OutputFiles (fileType, path, timeCreated, md5sum, sourceType) VALUES " > ${outputFilesAnalysisQueryFile}
  echo "INSERT INTO DSPRegressionTesting.AnalysisOutputFiles (analysis, outputFile) VALUES " > ${analysisOutputFilesQueryFile}
 
  echo "INSERT INTO DSPRegressionTesting.Metric_Concordance_FilterAnalysis(filter, tn, fn, uniqueTn, uniqueFn) VALUES " > ${metricsConcordanceFAQueryFile}
  echo "INSERT INTO DSPRegressionTesting.Metric_Concordance_Summary(variantType, truePositive, falsePositive, falseNegative, concordanceSensitivity, concordancePrecision) VALUES " > ${metricsConcordanceSummaryQueryFile}

  echo "INSERT INTO DSPRegressionTesting.Metric_GenotypeConcordance_ContingencyMetrics(variantType, truthSample, callSample, tpCount, tnCount, fpCount, fnCount, emptyCount) VALUES " > ${metricsGCContingencyQueryFile}
  echo "INSERT INTO DSPRegressionTesting.Metric_GenotypeConcordance_DetailMetrics(variantType, truthSample, callSample, truthState, callState, count, contingencyValues) VALUES " > ${metricsGCDetailQueryFile}
  echo "INSERT INTO DSPRegressionTesting.Metric_GenotypeConcordance_SummaryMetrics(variantType, truthSample, callSample, hetSensitivity, hetPPV, hetSpecificity, homvarSensitivity, homvarPPV, homvarSpecificity, varSensitivity, varPPV, varSpecificity, genotypeConcordance, nonRefGenotypeConcordance) VALUES " > ${metricsGCSummaryQueryFile}

  echo "INSERT INTO DSPRegressionTesting.Metric_Timing(startTime, endTime, elapsedTime) VALUES " > ${metricsTimingQueryFile}

  echo "INSERT INTO DSPRegressionTesting.Metrics(metricTableName, concreteMetricID, sourceAnalysis) VALUES " > ${metricsQueryFile}

  outputFileId=$( getNextOutputFileIDFromDB )
  isFirstFile=true
  isFirstAnalysisFile=true
  analysisRunID=$( getAnalysisRunIDFromDB )
 
  leadSepMCFA=''
  leadSepMCS=''
  leadSepMGCS=''
  leadSepMGCD=''
  leadSepMGCC=''
  leadSepMT=''
  leadSepM=''

  metricCFAID=$( getNextMetricCFAIDFromDB )
  metricCSID=$( getNextMetricCSIDFromDB )

  metricGCCID=$( getNextMetricGCCIDFromDB )
  metricGCDID=$( getNextMetricGCDIDFromDB )
  metricGCSID=$( getNextMetricGCSIDFromDB )

  metricTID=$( getNextMetricTimingIDFromDB )

  metricID=$( getNextMetricIDFromDB )

  for r in ${RUN_LIST} ; do

    echo "$(date +%Y%m%d-%H%M%S) - Processing run: ${r}"

    #################### 
    
    tmpFile=mainFileList.txt
    echo -n "  Getting outputs..."
    ~/Development/cromshell/cromshell list-outputs $r 2>/dev/null | grep -v 'call-HaplotypeCallerTask/.*vcf[\.idx]*$' | xargs gsutil ls -l | grep -v '^[ \t]*$'  | grep -v '^TOTAL' > $tmpFile
    echo 'DONE!'

    echo -n "  Creating shard / input map..."
    shardInputMapFile=$(mktemp)
    grep 'HaplotypeCallerTask' $tmpFile | grep HaplotypeCallerTask | sed 's#.*shard-#shard-#g' | sed 's#\(.*NexPond-[0-9]*\).*#\1#g' | tr '/' '\t' > ${shardInputMapFile}
    echo 'DONE!'

    #################### 

    echo -n "  Sorting outputs into Concordance, GenotypeConcordance, and Timing sets... "
    concordanceFile=concordanceFileList.txt
    genotypeConcordanceFile=genotypeConcordanceFileList.txt
    timingFile=timingFileList.txt
    grep 'call-Concordance' $tmpFile > $concordanceFile
    grep 'call-GenotypeConcordance' $tmpFile > $genotypeConcordanceFile
    grep 'call-HaplotypeCallerTask' $tmpFile > $timingFile 
    echo 'DONE!'

    #################### 

    echo "  $(date +%Y%m%d-%H%M%S) - Handling Concordance runs...  "
    lastShardIdx=999
    while read outFileInfo ; do
      
      outFile=$( echo $outFileInfo | awk '{print $3}' ) 
      localFile="${LOCAL_ROOT}/$(echo "${outFile}" | sed 's#.*ToolComparisonWdl/##g')"
      shardIdx=$( echo $outFile | grep -o shard-[0-9]* | sed 's#shard-##g' )
      
      echo "    Processing output file ${outputFileId} - $outFile"

      if [[ "${lastShardIdx}" != "${shardIdx}" ]] ; then
        let analysisRunID=$analysisRunID+1
        lastShardIdx=${shardIdx}
      fi

      timeStamp=$( echo $outFileInfo | awk '{print $2}' | tr -d 'Z' | tr 'T' ' ')
      truthFileBaseName=$(grep "shard-${shardIdx}" $shardInputMapFile | awk '{print $NF}' )
      inputFile="gs://broad-dsde-methods/cromwell-execution-36/ToolComparisonWdl/${r}/call-HaplotypeCallerTask/shard-${shardIdx}/${truthFileBaseName}.HC.vcf"
      truthFile="gs://broad-dsp-methods-regression-testing/inputData/${truthFileBaseName}.vcf.gz"

      originalOutFileID=$( getOutputFileID $truthFileBaseName $r )

      FT=$( getFileType $outFile )
      #md5sum=$( gsutil hash -hm $outFile 2>/dev/null | grep md5 | awk '{print $NF}'| \grep -o '^[0-9A-Za-z]*$' )
      md5sum=$( md5sum-lite ${localFile} | awk '{print $1}' )
      
      if ! $isFirstFile ; then
        echo -n ',' >> ${analysisOutputFilesQueryFile}
        echo -n ',' >> ${outputFilesAnalysisQueryFile}
      fi

      set +e
      grep -q "(1, NULL, '$(createConcordanceJson $inputFile $truthFile | jq -c .)')" $analysisRunsQueryFile
      rv=$?
      set -e
      if [ $rv -ne 0 ] ; then
        if ! $isFirstFile ; then 
          echo ",(1, NULL, '$(createConcordanceJson $inputFile $truthFile | jq -c .)')" >> $analysisRunsQueryFile 
        else
          echo "(1, NULL, '$(createConcordanceJson $inputFile $truthFile | jq -c .)')" >> $analysisRunsQueryFile 
        fi
      fi 

      echo "('$FT', '$outFile', '$timeStamp', '$md5sum', 'ANALYSIS')" >> ${outputFilesAnalysisQueryFile}
      echo "($analysisRunID, $outputFileId)" >> $analysisOutputFilesQueryFile
      
      set +e
      grep -q "($analysisRunID, $originalOutFileID)" $analysisOutputFilesQueryFile
      rv=$?
      set -e 
      [ $rv -ne 0 ] && echo ",($analysisRunID, $originalOutFileID)" >> $analysisOutputFilesQueryFile

      let outputFileId=${outputFileId}+1

      # Handle metrics:
      if [[ $outFile =~ .*filter-analysis.txt ]] ; then
        # Get filter analysis metric:
        #q=$( gsutil cat $outFile | tail -n+2 | sed -e 's#^\([0-9a-zA-Z]*\)\([ \t]*.*\)#"\1"\2#g' | tr '\t' ',' )
        q=$( cat $localFile | tail -n+2 | sed -e 's#^\([0-9a-zA-Z]*\)\([ \t]*.*\)#"\1"\2#g' | tr '\t' ',' )
        echo "${leadSepMCFA}($q)" >> $metricsConcordanceFAQueryFile 

        # Add to primary metrics table:
        echo "${leadSepM}('Metric_Concordance_FilterAnalysis', $metricCFAID, $analysisRunID)" >> $metricsQueryFile

        # Update metric counters:
        let metricCFAID=$metricCFAID+1
        let metricID=$metricID+1

        # Update metric lead separators:
        leadSepMCFA=','
        leadSepM=','

      elif [[ $outFile =~ .*/summary.txt ]] ; then
        while read line ; do 

          # get summary metric:
          q=$( echo $line | sed 's#^\([a-zA-Z0-9]*\)\([ \t]*.*\)#"\1"\2#g' | tr ' ' ',' )
          echo "${leadSepMCS}($q)" >> $metricsConcordanceSummaryQueryFile

          # Add to primary metrics table:
          echo "${leadSepM}('Metric_Concordance_Summary', $metricCSID, $analysisRunID)" >> $metricsQueryFile

          # Update metrics counters:
          let metricCSID=$metricCSID+1
          let metricID=$metricID+1

          # Update metric lead separators:
          leadSepMCS=','
          leadSepM=','

        #done < <(gsutil cat $outFile | tail -n+2)
        done < <(cat $localFile | tail -n+2)
      fi

      isFirstFile=false
    done < $concordanceFile

    #################### 

    echo "  $(date +%Y%m%d-%H%M%S) - Handling GenotypeConcordance runs...  "
    lastShardIdx=999
    while read outFileInfo ; do
      
      outFile=$( echo $outFileInfo | awk '{print $3}' ) 
      localFile="${LOCAL_ROOT}/$(echo "${outFile}" | sed 's#.*ToolComparisonWdl/##g')"
      shardIdx=$( echo $outFile | grep -o shard-[0-9]* | sed 's#shard-##g' )

      echo "    Processing output file ${outputFileId} - $outFile"

      if [[ "${lastShardIdx}" != "${shardIdx}" ]] ; then
        let analysisRunID=$analysisRunID+1
        lastShardIdx=${shardIdx}
      fi

      timeStamp=$( echo $outFileInfo | awk '{print $2}' | tr -d 'Z' | tr 'T' ' ')
      truthFileBaseName=$(grep "shard-${shardIdx}" $shardInputMapFile | awk '{print $NF}' )
      inputFile="gs://broad-dsde-methods/cromwell-execution-36/ToolComparisonWdl/${r}/call-HaplotypeCallerTask/shard-${shardIdx}/${truthFileBaseName}.HC.vcf"
      truthFile="gs://broad-dsp-methods-regression-testing/inputData/${truthFileBaseName}.vcf.gz"
      originalOutFileID=$( getOutputFileID $truthFileBaseName $r )
      FT=$( getFileType $outFile )
      #md5sum=$( gsutil hash -hm $outFile 2>/dev/null | grep md5 | awk '{print $NF}'| \grep -o '^[0-9A-Za-z]*$' )
      md5sum=$( md5sum-lite ${localFile} | awk '{print $1}' )
      
      set +e
      grep -q "(2, NULL, '$(createGenotypeConcordanceJson $inputFile $truthFile $truthFileBaseName | jq -c .)')"  $analysisRunsQueryFile
      rv=$?
      set -e
      [ $rv -ne 0 ] && echo ",(2, NULL, '$(createGenotypeConcordanceJson $inputFile $truthFile $truthFileBaseName | jq -c .)')" >> $analysisRunsQueryFile 

      set +e
      grep -q "($analysisRunID, $originalOutFileID)" $analysisOutputFilesQueryFile
      rv=$?
      set -e
      [ $rv -ne 0 ] && echo ",($analysisRunID, $originalOutFileID)" >> $analysisOutputFilesQueryFile
      
      echo ",('$FT', '$outFile', '$timeStamp', '$md5sum', 'ANALYSIS')" >> ${outputFilesAnalysisQueryFile}
      echo ",($analysisRunID, $outputFileId)" >> $analysisOutputFilesQueryFile

      let outputFileId=${outputFileId}+1

      # Handle metrics:
      if [[ $outFile =~ .*genotype_concordance_contingency_metrics ]] ; then
        # Get GC contingency metrics: 
        while read line ; do 
          echo "${leadSepMGCC}($line)" >> $metricsGCContingencyQueryFile

          # Add to primary metrics table:
          echo "${leadSepM}('Metric_GenotypeConcordance_ContingencyMetrics', $metricGCCID, $analysisRunID)" >> $metricsQueryFile

          # Update metric counters:
          let metricGCCID=$metricGCCID+1
          let metricID=$metricID+1

          # Update metric lead separators:
          leadSepMGCC=','
          leadSepM=','
        #done < <(gsutil cat $outFile | tail -n+8 | grep -v '^[ \t]*$' | tr '\t' ',' | sed 's#^\([a-zA-Z0-9]*\),\([a-zA-Z0-9]*\),\([a-zA-Z0-9]*\),#"\1","\2","\3",#g' )
        done < <(cat $localFile | tail -n+8 | grep -v '^[ \t]*$' | tr '\t' ',' | sed 's#^\([a-zA-Z0-9]*\),\([a-zA-Z0-9]*\),\([a-zA-Z0-9]*\),#"\1","\2","\3",#g' )

      elif [[ $outFile =~ .*genotype_concordance_detail_metrics ]] ; then
        # Get GC contingency metrics: 
        while read line ; do 
          echo "${leadSepMGCD}($line)" >> $metricsGCDetailQueryFile

          # Add to primary metrics table:
          echo "${leadSepM}('Metric_GenotypeConcordance_DetailMetrics', $metricGCDID, $analysisRunID)" >> $metricsQueryFile

          # Update metric counters:
          let metricGCDID=$metricGCDID+1
          let metricID=$metricID+1

          # Update metric lead separators:
          leadSepMGCD=','
          leadSepM=','
        #done < <(gsutil cat $outFile | tail -n+8 | grep -v '^[ \t]*$' | tr ',' ';' | tr '\t' ',' | sed -e "s#,#','#g" -e "s#^#'#g" -e "s#\$#'#g" | tr ';' ',' ) 
        done < <(cat $localFile | tail -n+8 | grep -v '^[ \t]*$' | tr ',' ';' | tr '\t' ',' | sed -e "s#,#','#g" -e "s#^#'#g" -e "s#\$#'#g" | tr ';' ',' ) 

      elif [[ $outFile =~ .*genotype_concordance_summary_metrics ]] ; then
        # Get GC Summary metrics: 
        while read line ; do 
          echo "${leadSepMGCS}($line)" >> $metricsGCSummaryQueryFile

          # Add to primary metrics table:
          echo "${leadSepM}('Metric_GenotypeConcordance_SummaryMetrics', $metricGCSID, $analysisRunID)" >> $metricsQueryFile

          # Update metric counters:
          let metricGCSID=$metricGCSID+1
          let metricID=$metricID+1

          # Update metric lead separators:
          leadSepMGCS=','
          leadSepM=','
        #done < <(gsutil cat $outFile | tail -n+8 | grep -v '^[ \t]*$' | tr '\t' ',' | sed -e "s#,#','#g" -e "s#^#'#g" -e "s#\$#'#g" -e "s#'\?'#NULL#g" )
        done < <(cat $localFile | tail -n+8 | grep -v '^[ \t]*$' | tr '\t' ',' | sed -e "s#,#','#g" -e "s#^#'#g" -e "s#\$#'#g" -e "s#'\?'#NULL#g" )
       fi

    done < $genotypeConcordanceFile

    #################### 

    echo "  $(date +%Y%m%d-%H%M%S) - Handling Timing runs...  "
    lastShardIdx=999
    while read outFileInfo ; do
      
      outFile=$( echo $outFileInfo | awk '{print $3}' ) 
      localFile="${LOCAL_ROOT}/$(echo "${outFile}" | sed 's#.*ToolComparisonWdl/##g')"
      shardIdx=$( echo $outFile | grep -o shard-[0-9]* | sed 's#shard-##g' )
      
      echo "    Processing output file ${outputFileId} - $outFile"

      timeStamp=$( echo $outFileInfo | awk '{print $2}' | tr -d 'Z' | tr 'T' ' ')
      truthFileBaseName=$(grep "shard-${shardIdx}" $shardInputMapFile | awk '{print $NF}' )
      inputFile="gs://broad-dsde-methods/cromwell-execution-36/ToolComparisonWdl/${r}/call-HaplotypeCallerTask/shard-${shardIdx}/${truthFileBaseName}.HC.vcf"
      truthFile="gs://broad-dsp-methods-regression-testing/inputData/${truthFileBaseName}.vcf.gz"
      originalOutFileID=$( getOutputFileID $truthFileBaseName $r )
      FT=$( getFileType $outFile )
      #md5sum=$( gsutil hash -hm $outFile 2>/dev/null | grep md5 | awk '{print $NF}'| \grep -o '^[0-9A-Za-z]*$' )
      md5sum=$( md5sum-lite ${localFile} | awk '{print $1}' )

      if [[ "${lastShardIdx}" != "${shardIdx}" ]] ; then
        let analysisRunID=$analysisRunID+1
        echo ",($analysisRunID, $originalOutFileID)" >> $analysisOutputFilesQueryFile
        echo ",(3, NULL, NULL)" >> $analysisRunsQueryFile 
        lastShardIdx=${shardIdx}
      fi

      echo ",('$FT', '$outFile', '$timeStamp', '$md5sum', 'ANALYSIS')" >> ${outputFilesAnalysisQueryFile}
      echo ",($analysisRunID, $outputFileId)" >> $analysisOutputFilesQueryFile

      let outputFileId=${outputFileId}+1

      # Handle metrics:
      if [[ $outFile =~ .*.timingInformation.txt ]] ; then
        # Get timing metric:
        #q=$( gsutil cat $outFile | sed -e 's#.*:[ \t]*##g' | tr '\n' ','|sed 's#.$##g' )
        q=$( cat $localFile | sed -e 's#.*:[ \t]*##g' | tr '\n' ','|sed 's#.$##g' )
        echo "${leadSepMT}($q)" >> $metricsTimingQueryFile

        # Add to primary metrics table:
        echo "${leadSepM}('Metric_Timing', $metricTID, $analysisRunID)" >> $metricsQueryFile

        # Update metric counters:
        let metricTID=$metricTID+1
        let metricID=$metricID+1

        # Update metric lead separators:
        leadSepMT=','
        leadSepM=','
      fi

      isFirstFile=false
    done < $timingFile

    ##################### 

  done

  # End our SQL queries:
  echo ";" >> ${analysisRunsQueryFile}
  echo ";" >> ${outputFilesAnalysisQueryFile}
  echo ";" >> ${analysisOutputFilesQueryFile}
  echo ";" >> ${metricsConcordanceFAQueryFile}
  echo ";" >> ${metricsConcordanceSummaryQueryFile}
  echo ";" >> ${metricsGCContingencyQueryFile}
  echo ";" >> ${metricsGCDetailQueryFile}
  echo ";" >> ${metricsGCSummaryQueryFile}
  echo ";" >> ${metricsTimingQueryFile}
  echo ";" >> ${metricsQueryFile}

  # Consolidate all files together into a single sql file:
  cat A* > Z_MASTER_Analysis_File.SQL

  echo "Done processing ANALYSIS output files."
fi

set +e


