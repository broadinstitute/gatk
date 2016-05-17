#include "headers.h"
#include "jni_common.h"
#include "org_broadinstitute_hellbender_utils_pairhmm_VectorLoglessPairHMM.h"
#include "utils.h"
#include "LoadTimeInitializer.h"

using namespace std;

bool use_double = false;

//Should be called only once for the whole Java process - initializes field ids for the classes JNIReadDataHolderClass
//and JNIHaplotypeDataHolderClass
JNIEXPORT void JNICALL Java_org_broadinstitute_hellbender_utils_pairhmm_VectorLoglessPairHMM_jniInitializeClassFields
  (JNIEnv* env, jobject thisObject, jclass readDataHolderClass, jclass haplotypeDataHolderClass)
{
  assert(readDataHolderClass);
  assert(haplotypeDataHolderClass);
  jfieldID fid;
  fid = env->GetFieldID(readDataHolderClass, "readBases", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for readBases");
  g_load_time_initializer.m_readBasesFID = fid;
  fid = env->GetFieldID(readDataHolderClass, "readQuals", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for readQuals");
  g_load_time_initializer.m_readQualsFID = fid;
  fid = env->GetFieldID(readDataHolderClass, "insertionGOP", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for insertionGOP");
  g_load_time_initializer.m_insertionGOPFID = fid;
  fid = env->GetFieldID(readDataHolderClass, "deletionGOP", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for deletionGOP");
  g_load_time_initializer.m_deletionGOPFID = fid;
  fid = env->GetFieldID(readDataHolderClass, "overallGCP", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for overallGCP");
  g_load_time_initializer.m_overallGCPFID = fid;

  fid = env->GetFieldID(haplotypeDataHolderClass, "haplotypeBases", "[B");
  assert(fid && "JNI pairHMM: Could not get FID for haplotypeBases");
  g_load_time_initializer.m_haplotypeBasesFID = fid;
}

JNIEXPORT void JNICALL initializeHaplotypes
  (JNIEnv * env, jobject& thisObject, jint numHaplotypes, jobjectArray& haplotypeDataArray,
   vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector, vector<unsigned>& haplotypeBasesLengths)
{
  jboolean is_copy = JNI_FALSE;
  haplotypeBasesArrayVector.clear();
  haplotypeBasesLengths.clear();
  haplotypeBasesArrayVector.resize(numHaplotypes);
  haplotypeBasesLengths.resize(numHaplotypes);
  jsize haplotypeBasesLength = 0;
  for(unsigned j=0;j<numHaplotypes;++j)
  {
    jobject haplotypeObject = env->GetObjectArrayElement(haplotypeDataArray, j);
    jbyteArray haplotypeBases = (jbyteArray)env->GetObjectField(haplotypeObject, g_load_time_initializer.m_haplotypeBasesFID);
#ifdef ENABLE_ASSERTIONS
    assert(haplotypeBases && ("haplotypeBases is NULL at index : "+to_string(j)+"\n").c_str());
#endif
    //Need a global reference as this will be accessed across multiple JNI calls to JNIComputeLikelihoods()
    jbyteArray haplotypeBasesGlobalRef = (jbyteArray)env->NewGlobalRef(haplotypeBases);
#ifdef ENABLE_ASSERTIONS
    assert(haplotypeBasesGlobalRef && ("Could not get global ref to haplotypeBases at index : "+to_string(j)+"\n").c_str());
#endif
    env->DeleteLocalRef(haplotypeBases);	//free the local reference
    jbyte* haplotypeBasesArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(haplotypeBasesGlobalRef, &is_copy);
    haplotypeBasesLength = env->GetArrayLength(haplotypeBasesGlobalRef);
#ifdef ENABLE_ASSERTIONS
    assert(haplotypeBasesArray && "haplotypeBasesArray not initialized in JNI"); 
    //assert(haplotypeBasesLength < MCOLS);
#endif
#ifdef DEBUG
    cout << "JNI haplotype length "<<haplotypeBasesLength<<"\n";
#endif
    haplotypeBasesArrayVector[j] = make_pair(haplotypeBasesGlobalRef, haplotypeBasesArray);
    haplotypeBasesLengths[j] = haplotypeBasesLength;
#ifdef DEBUG
    for(unsigned k=0;k<haplotypeBasesLength;++k)
      g_load_time_initializer.debug_dump("haplotype_bases_jni.txt",to_string((int)haplotypeBasesArray[k]),true);
#endif
  }
}

JNIEXPORT void JNICALL releaseHaplotypes(JNIEnv * env, jobject thisObject,
    vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector, vector<unsigned>& haplotypeBasesLengths
    )
{
  //Now release haplotype arrays
  for(int j=haplotypeBasesArrayVector.size()-1;j>=0;--j)	//note the order - reverse of GET
  {
    RELEASE_BYTE_ARRAY_ELEMENTS(haplotypeBasesArrayVector[j].first, haplotypeBasesArrayVector[j].second, JNI_RO_RELEASE_MODE);
    env->DeleteGlobalRef(haplotypeBasesArrayVector[j].first);	//free the global reference
  }
  haplotypeBasesArrayVector.clear();
  haplotypeBasesLengths.clear(); 
}

//Create a vector of testcases for computation - copy the references to bytearrays read/readQuals etc into the appropriate
//testcase struct
inline JNIEXPORT void JNICALL Java_org_broadinstitute_hellbender_utils_pairhmm_VectorLoglessPairHMM_jniInitializeTestcasesVector
  (JNIEnv* env, jint numReads, jint numHaplotypes, jobjectArray& readDataArray,
   vector<vector<pair<jbyteArray,jbyte*> > >& readBasesArrayVector,
   vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector, vector<unsigned>& haplotypeBasesLengths,
   vector<testcase>& tc_array)
{
  jboolean is_copy = JNI_FALSE;
  unsigned tc_idx = 0;
  for(unsigned i=0;i<numReads;++i)
  {
    //Get bytearray fields from read
    jobject readObject = env->GetObjectArrayElement(readDataArray, i);
    jbyteArray readBases = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_readBasesFID);
    jbyteArray insertionGOP = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_insertionGOPFID);
    jbyteArray deletionGOP = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_deletionGOPFID);
    jbyteArray overallGCP = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_overallGCPFID);
    jbyteArray readQuals = (jbyteArray)env->GetObjectField(readObject, g_load_time_initializer.m_readQualsFID);

#ifdef ENABLE_ASSERTIONS
    assert(readBases && ("readBases is NULL at index : "+to_string(i)+"\n").c_str());
    assert(insertionGOP && ("insertionGOP is NULL at index : "+to_string(i)+"\n").c_str());
    assert(deletionGOP && ("deletionGOP is NULL at index : "+to_string(i)+"\n").c_str());
    assert(overallGCP && ("overallGCP is NULL at index : "+to_string(i)+"\n").c_str());
    assert(readQuals && ("readQuals is NULL at index : "+to_string(i)+"\n").c_str());
#endif
    jsize readLength = env->GetArrayLength(readBases);

    jbyte* readBasesArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(readBases, &is_copy);	//order of GET-RELEASE is important
    jbyte* readQualsArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(readQuals, &is_copy);
    jbyte* insertionGOPArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(insertionGOP, &is_copy);
    jbyte* deletionGOPArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(deletionGOP, &is_copy);
    jbyte* overallGCPArray = (jbyte*)GET_BYTE_ARRAY_ELEMENTS(overallGCP, &is_copy);
#ifdef ENABLE_ASSERTIONS
    assert(readBasesArray && "readBasesArray not initialized in JNI"); 
    assert(readQualsArray && "readQualsArray not initialized in JNI"); 
    assert(insertionGOPArray && "insertionGOP array not initialized in JNI");
    assert(deletionGOPArray && "deletionGOP array not initialized in JNI");
    assert(overallGCPArray && "overallGCP array not initialized in JNI");
    assert(readLength == env->GetArrayLength(readQuals));
    assert(readLength == env->GetArrayLength(insertionGOP));
    assert(readLength == env->GetArrayLength(deletionGOP));
    assert(readLength == env->GetArrayLength(overallGCP));
#endif
#ifdef DEBUG
    cout << "JNI read length "<<readLength<<"\n";
#endif
#ifdef DEBUG
    for(unsigned j=0;j<readLength;++j)
    {
      g_load_time_initializer.debug_dump("reads_jni.txt",to_string((int)readBasesArray[j]),true);
      g_load_time_initializer.debug_dump("reads_jni.txt",to_string((int)readQualsArray[j]),true);
      g_load_time_initializer.debug_dump("reads_jni.txt",to_string((int)insertionGOPArray[j]),true);
      g_load_time_initializer.debug_dump("reads_jni.txt",to_string((int)deletionGOPArray[j]),true);
      g_load_time_initializer.debug_dump("reads_jni.txt",to_string((int)overallGCPArray[j]),true);
    }
#endif
    for(unsigned j=0;j<numHaplotypes;++j)
    {
      jsize haplotypeLength = (jsize)haplotypeBasesLengths[j];
      jbyte* haplotypeBasesArray = haplotypeBasesArrayVector[j].second;
      tc_array[tc_idx].rslen = (int)readLength;
      tc_array[tc_idx].haplen = (int)haplotypeLength;
      tc_array[tc_idx].hap = (char*)haplotypeBasesArray;
      tc_array[tc_idx].rs = (char*)readBasesArray;
      tc_array[tc_idx].q = (char*)readQualsArray;
      tc_array[tc_idx].i = (char*)insertionGOPArray;
      tc_array[tc_idx].d = (char*)deletionGOPArray;
      tc_array[tc_idx].c = (char*)overallGCPArray;
#ifdef DUMP_TO_SANDBOX
      g_load_time_initializer.dump_sandbox(tc_array[tc_idx], tc_idx);
#endif
      ++tc_idx;  
    }
    //Store the read array references and release them at the end because they are used by compute_full_prob
    //Maintain order in which GET_BYTE_ARRAY_ELEMENTS called
    readBasesArrayVector[i].clear();
    readBasesArrayVector[i].resize(5);
    readBasesArrayVector[i][0] = make_pair(readBases, readBasesArray);
    readBasesArrayVector[i][1] = make_pair(readQuals, readQualsArray);
    readBasesArrayVector[i][2] = make_pair(insertionGOP, insertionGOPArray);
    readBasesArrayVector[i][3] = make_pair(deletionGOP, deletionGOPArray);
    readBasesArrayVector[i][4] = make_pair(overallGCP, overallGCPArray);
  }
}

//Do compute over vector of testcase structs
inline void compute_testcases(vector<testcase>& tc_array, unsigned numTestCases, double* likelihoodDoubleArray,
    unsigned maxNumThreadsToUse)
{
  #pragma omp parallel for schedule (dynamic,10000) num_threads(maxNumThreadsToUse)
  for(unsigned tc_idx=0;tc_idx<numTestCases;++tc_idx)
  {
    float result_avxf = use_double ? 0 : g_compute_full_prob_float(&(tc_array[tc_idx]), 0);
    double result = 0;
    if (result_avxf < MIN_ACCEPTED) {
      double result_avxd = g_compute_full_prob_double(&(tc_array[tc_idx]), 0);
      result = log10(result_avxd) - log10(ldexp(1.0, 1020.0));
    }
    else
      result = (double)(log10f(result_avxf) - log10f(ldexpf(1.f, 120.f)));
#ifdef DUMP_TO_SANDBOX
    g_load_time_initializer.dump_sandbox_result(result);    
#endif
    likelihoodDoubleArray[tc_idx] = result;
  }
}

//Inform the Java VM that we no longer need access to the read arrays (and free memory)
inline JNIEXPORT void JNICALL Java_org_broadinstitute_hellbender_utils_pairhmm_VectorLoglessPairHMM_jniReleaseReadArrays
  (JNIEnv* env, vector<vector<pair<jbyteArray,jbyte*> > >& readBasesArrayVector)
{
  //Release read arrays first
  for(int i=readBasesArrayVector.size()-1;i>=0;--i)//note the order - reverse of GET
  {
    for(int j=readBasesArrayVector[i].size()-1;j>=0;--j)
      RELEASE_BYTE_ARRAY_ELEMENTS(readBasesArrayVector[i][j].first, readBasesArrayVector[i][j].second, JNI_RO_RELEASE_MODE);
    readBasesArrayVector[i].clear();
  }
  readBasesArrayVector.clear();
}


//JNI function to invoke compute_full_prob_avx
//readDataArray - array of JNIReadDataHolderClass objects which contain the readBases, readQuals etc
//haplotypeDataArray - array of JNIHaplotypeDataHolderClass objects which contain the haplotypeBases
//likelihoodArray - array of doubles to return results back to Java. Memory allocated by Java prior to JNI call
//maxNumThreadsToUse - Max number of threads that OpenMP can use for the HMM computation
JNIEXPORT void JNICALL Java_org_broadinstitute_hellbender_utils_pairhmm_VectorLoglessPairHMM_jniComputeLikelihoods
  (JNIEnv* env, jobject thisObject, jint numReads, jint numHaplotypes, 
   jobjectArray readDataArray, jobjectArray haplotypeDataArray, jdoubleArray likelihoodArray, jint maxNumThreadsToUse)
{
  //Very important to get good performance on Intel processors
  //Function: enabling FTZ converts denormals to 0 in hardware
  //Denormals cause microcode to insert uops into the core causing big slowdown
  //NOTE: To protect against this flag being reset on us, we're going to set it on every call to jniComputeLikelihoods
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

#ifdef DEBUG
  cout << "JNI numReads "<<numReads<<" numHaplotypes "<<numHaplotypes<<"\n";
#endif
  jboolean is_copy = JNI_FALSE;
  struct timespec start_time;
  unsigned numTestCases = numReads*numHaplotypes;
  //vector to store testcases
  vector<testcase> tc_array;
  tc_array.clear();
  tc_array.resize(numTestCases);
  //Store read arrays for release later
  vector<vector<pair<jbyteArray,jbyte*> > > readBasesArrayVector;
  readBasesArrayVector.clear();
  readBasesArrayVector.resize(numReads);
#ifdef DUMP_TO_SANDBOX
  g_load_time_initializer.open_sandbox();
#endif

  vector<pair<jbyteArray, jbyte*> > l_haplotypeBasesArrayVector;
  vector<pair<jbyteArray, jbyte*> >& haplotypeBasesArrayVector = l_haplotypeBasesArrayVector;
  vector<unsigned> l_haplotypeBasesLengths;
  vector<unsigned>& haplotypeBasesLengths = l_haplotypeBasesLengths;
  initializeHaplotypes(env, thisObject, numHaplotypes, haplotypeDataArray, haplotypeBasesArrayVector, haplotypeBasesLengths);

  //Copy byte array references from Java memory into vector of testcase structs
  Java_org_broadinstitute_hellbender_utils_pairhmm_VectorLoglessPairHMM_jniInitializeTestcasesVector(env,
      numReads, numHaplotypes, readDataArray, readBasesArrayVector, haplotypeBasesArrayVector, haplotypeBasesLengths, tc_array);

  //Get double array where results are stored (to pass back to java)
  jdouble* likelihoodDoubleArray = (jdouble*)GET_DOUBLE_ARRAY_ELEMENTS(likelihoodArray, &is_copy);
#ifdef ENABLE_ASSERTIONS
  assert(likelihoodDoubleArray && "likelihoodArray is NULL");
  assert(env->GetArrayLength(likelihoodArray) == numTestCases);
#endif
  compute_testcases(tc_array, numTestCases, likelihoodDoubleArray, maxNumThreadsToUse); //actual computation
#ifdef DUMP_COMPUTE_VALUES
  for(unsigned tc_idx=0;tc_idx<numTestCases;++tc_idx)
    g_load_time_initializer.debug_dump("return_values_jni.txt",to_string(likelihoodDoubleArray[tc_idx]),true);
#endif
  RELEASE_DOUBLE_ARRAY_ELEMENTS(likelihoodArray, likelihoodDoubleArray, 0); //release mode 0, copy back results to Java memory (if copy made)
  Java_org_broadinstitute_hellbender_utils_pairhmm_VectorLoglessPairHMM_jniReleaseReadArrays(env, readBasesArrayVector);
  releaseHaplotypes(env, thisObject, haplotypeBasesArrayVector, haplotypeBasesLengths);

  tc_array.clear();
#ifdef DUMP_TO_SANDBOX
  g_load_time_initializer.close_sandbox();
#endif
}

JNIEXPORT void JNICALL Java_org_broadinstitute_hellbender_utils_pairhmm_VectorLoglessPairHMM_jniClose
  (JNIEnv* env, jobject thisObject)
{
#ifdef DUMP_COMPUTE_VALUES
  g_load_time_initializer.debug_close();
#endif
}
