#include <iomanip>
#include "utils.h"
#include "LoadTimeInitializer.h"
using namespace std;
const char* LoadTimeInitializerStatsNames[] = 
{
  "num_regions",
  "num_reads",
  "num_haplotypes",
  "num_testcases",
  "num_double_invocations",
  "haplotype_length",
  "readlength",
  "product_read_length_haplotype_length",
  "dummy"
};

LoadTimeInitializer g_load_time_initializer;

LoadTimeInitializer::LoadTimeInitializer()		//will be called when library is loaded
{
  ConvertChar::init();
  //Very important to get good performance on Intel processors
  //Function: enabling FTZ converts denormals to 0 in hardware
  //Denormals cause microcode to insert uops into the core causing big slowdown
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

  //Profiling: times for compute and transfer (either bytes copied or pointers copied)
  m_compute_time = 0;
  m_data_transfer_time = 0;
  m_bytes_copied = 0;

  //Initialize profiling counters
  for(unsigned i=0;i<TOTAL_NUMBER_STATS;++i)
  {
    m_sum_stats[i] = 0;
    m_sum_square_stats[i] = 0;
    m_max_stats[i] = 0;
    m_min_stats[i] = 0xFFFFFFFFFFFFFFFFull;
  }

  //for debug dump
  m_filename_to_fptr.clear();
  m_written_files_set.clear();

  initialize_function_pointers();

  //Initialize static members of class
  Context<float>::initializeStaticMembers();
  Context<double>::initializeStaticMembers();

  cout.flush();
}

void LoadTimeInitializer::print_profiling()
{
  double mean = 0;
  double variance = 0;
  uint64_t denominator = 1;
  cout << "Time spent in compute_testcases "<<m_compute_time*1e-9<<"\n";
  cout << "Time spent in data transfer (Java <--> C++) "<<m_data_transfer_time*1e-9<<"\n";

  cout << "\nHC input stats\nstat_name,sum,sum_square,mean,variance,min,max\n";
  for(unsigned i=0;i<TOTAL_NUMBER_STATS;++i)
  {
    cout << LoadTimeInitializerStatsNames[i];
    cout << "," << m_sum_stats[i];
    cout << "," << std::scientific << m_sum_square_stats[i];
    denominator = 1;
    switch(i)
    {
      case NUM_READS_IDX:
      case NUM_HAPLOTYPES_IDX:
      case NUM_TESTCASES_IDX:
        denominator = m_sum_stats[NUM_REGIONS_IDX];
        break;
      case HAPLOTYPE_LENGTH_IDX:
        denominator = m_sum_stats[NUM_HAPLOTYPES_IDX];
        break;
      case READ_LENGTH_IDX:
        denominator = m_sum_stats[NUM_READS_IDX];
        break;
      case PRODUCT_READ_LENGTH_HAPLOTYPE_LENGTH_IDX:
        denominator = m_sum_stats[NUM_TESTCASES_IDX];
        break;
      default:
        denominator = 1;
        break;
    }
    mean = ((double)m_sum_stats[i])/denominator;
    cout << "," << std::scientific << mean;
    variance = (m_sum_square_stats[i]/denominator) - (mean*mean);       //E(X^2)-(E(X))^2
    cout << "," << std::scientific << variance;
    cout << "," << m_min_stats[i];
    cout << "," << m_max_stats[i];
    cout << "\n";
  }
  cout << "\n";
  cout.flush();
}

void LoadTimeInitializer::debug_dump(string filename, string s, bool to_append, bool add_newline)
{
  map<string, ofstream*>::iterator mI = m_filename_to_fptr.find(filename);
  ofstream* fptr = 0;
  if(mI == m_filename_to_fptr.end())
  {
    m_filename_to_fptr[filename] = new ofstream();
    fptr = m_filename_to_fptr[filename];
    //File never seen before
    if(m_written_files_set.find(filename) == m_written_files_set.end())
    {
      to_append = false;
      m_written_files_set.insert(filename);
    }
    fptr->open(filename.c_str(), to_append ? ios::app : ios::out);
    assert(fptr->is_open());
  }
  else
    fptr = (*mI).second;
  //ofstream fptr;
  //fptr.open(filename.c_str(), to_append ? ofstream::app : ofstream::out);
  (*fptr) << s;
  if(add_newline)
    (*fptr) << "\n";
  //fptr.close();
}
void LoadTimeInitializer::debug_close()
{
  for(map<string,ofstream*>::iterator mB = m_filename_to_fptr.begin(), mE = m_filename_to_fptr.end();
      mB != mE;mB++)
  {
    (*mB).second->close();
    delete (*mB).second;
  }
  m_filename_to_fptr.clear();
}

void LoadTimeInitializer::dump_sandbox(testcase& tc, unsigned tc_idx)
{
  unsigned haplotypeLength = tc.haplen;
  unsigned readLength = tc.rslen;
  ofstream& dumpFptr = m_sandbox_fptr;
  for(unsigned k=0;k<haplotypeLength;++k)
    dumpFptr<<(char)(tc.hap[k]);
  dumpFptr<<" ";
  for(unsigned k=0;k<readLength;++k)
    dumpFptr<<(char)(tc.rs[k]);
  dumpFptr<<" ";
  for(unsigned k=0;k<readLength;++k)
    dumpFptr<<(char)(tc.q[k]+33);
  dumpFptr<<" ";
  for(unsigned k=0;k<readLength;++k)
    dumpFptr<<(char)(tc.i[k]+33);
  dumpFptr<<" ";
  for(unsigned k=0;k<readLength;++k)
    dumpFptr<<(char)(tc.d[k]+33);
  dumpFptr<<" ";
  for(unsigned k=0;k<readLength;++k)
    dumpFptr<<(char)(tc.c[k]+33);
  dumpFptr<<" ";
}

void LoadTimeInitializer::dump_sandbox_result(double result)
{
  ofstream& dumpFptr = m_sandbox_fptr;
  dumpFptr << std::scientific << std::setprecision(6) << result << "\n";
}

void LoadTimeInitializer::update_stat(LoadTimeInitializerStatsEnum stat_idx, uint64_t value)
{
  m_sum_stats[stat_idx] += value;
  double v = value;
  m_sum_square_stats[stat_idx] += (v*v);
  m_max_stats[stat_idx] = std::max(m_max_stats[stat_idx], value);
  m_min_stats[stat_idx] = std::min(m_min_stats[stat_idx], value);
}
