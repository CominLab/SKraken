#include "kraken_headers.hpp"
#include "quickfile.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "seqreader.hpp"
#include <math.h> 

#define SKIP_LEN 50000

using namespace std;
using namespace kraken;

void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
void process_files();
void process_single_file();
void process_file(string filename, uint32_t taxid);
void set_lcas(uint32_t taxid, string &seq, size_t start, size_t finish);
uint32_t count_children (uint32_t taxid);
void calculate_percentage();
uint32_t at_rank(uint32_t taxid, uint8_t rank);

string DB_filename, Index_filename, Nodes_filename, Score_filename,
  File_to_taxon_map_filename, Output_DB_filename, Index_filename_reduced,
  ID_to_taxon_map_filename, Multi_fasta_filename;
bool Allow_extra_kmers = false;
bool Operate_in_RAM = false;

bool One_FASTA_file = false;
map<uint32_t, uint32_t> Parent_map;
map<string, uint32_t> ID_to_taxon_map;
KrakenDB Database;
uint64_t del_cnt = 0;
uint32_t *kmer_weight;
uint8_t *kmer_perc;
uint8_t perc_score = 0;
uint8_t rank_level = 7;
uint8_t Bin_key_nt = 15;
map<uint32_t, uint32_t> Pruned_parent_map;
map<uint32_t, uint8_t> Rank_map;
multimap<uint32_t, uint32_t> Children_multimap;
std::set<uint64_t> kmer_check_set;

int main(int argc, char **argv) {
  parse_command_line(argc, argv);
  Parent_map = build_parent_map(Nodes_filename);
  Rank_map = build_rank_map(Nodes_filename);

  QuickFile db_file(DB_filename, "r");
  Database = KrakenDB(db_file.ptr());
  KmerScanner::set_k(Database.get_k());

  char *temp_ptr = NULL;
  size_t db_file_size = db_file.size();
  if (Operate_in_RAM) {
    db_file.close_file();
    temp_ptr = new char[ db_file_size ];
    ifstream ifs(DB_filename.c_str(), ifstream::binary);
    ifs.read(temp_ptr, db_file_size);
    ifs.close();
    Database = KrakenDB(temp_ptr);
  }

  QuickFile idx_file(Index_filename);
  KrakenDBIndex db_index(idx_file.ptr());
  Database.set_index(&db_index);

  uint64_t key_ct = Database.get_key_ct();
  kmer_weight = new uint32_t [key_ct];

  for (uint64_t i = 0; i < key_ct; i++)
    kmer_weight[i] = 0;
  kmer_perc = new uint8_t [key_ct];

  if (One_FASTA_file)
    process_single_file();
  else
    process_files();
  if (Operate_in_RAM) {
    ofstream ofs(DB_filename.c_str(), ofstream::binary);
    ofs.write(temp_ptr, db_file_size);
    ofs.close();
    delete temp_ptr;
  }

  Children_multimap = build_children_multimap (Pruned_parent_map);
  count_children (1);
  std::cout << "Calculating Score..." << std::endl;
  calculate_percentage();

  return 0;
}

void process_single_file() {
  ifstream map_file(ID_to_taxon_map_filename.c_str());
  if (map_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", ID_to_taxon_map_filename.c_str());
  }
  string line;
  while (map_file.good()) {
    getline(map_file, line);
    if (line.empty())
      break;
    string seq_id;
    uint32_t taxid;
    istringstream iss(line);
    iss >> seq_id;
    iss >> taxid;
    ID_to_taxon_map[seq_id] = taxid;
  }

  FastaReader reader(Multi_fasta_filename);
  DNASequence dna;
  uint32_t seqs_processed = 0;
  uint32_t lvl_taxid = 0;
  uint32_t temp_taxid = 0;
  uint64_t taxid_count = 0;
  bool start = true;

  while (reader.is_valid()) {
    dna = reader.next_sequence();
    if (! reader.is_valid())
      break;
    uint32_t taxid = ID_to_taxon_map[dna.id];

    lvl_taxid = at_rank(taxid, rank_level);
    if (start) {
      temp_taxid = lvl_taxid;
      taxid_count++;
    }
    if (lvl_taxid != temp_taxid) {
      temp_taxid = lvl_taxid;
      taxid_count++;
      kmer_check_set.clear();
    }

    if (taxid) {
      for (size_t i = 0; i < dna.seq.size(); i += SKIP_LEN) {
        set_lcas(taxid, dna.seq, i, i + SKIP_LEN + Database.get_k() - 1);
      }
    }
    cerr << "\rProcessed " << ++seqs_processed << " sequences";
    start = false;
  }

  cerr << "\r                                                       ";
  cerr << "\rFinished processing " << seqs_processed << " sequences" << endl;
}

void process_files() {
  ifstream map_file(File_to_taxon_map_filename.c_str());
  if (map_file.rdstate() & ifstream::failbit) {
    err(EX_NOINPUT, "can't open %s", File_to_taxon_map_filename.c_str());
  }
  string line;
  uint32_t seqs_processed = 0;

  while (map_file.good()) {
    getline(map_file, line);
    if (line.empty())
      break;
    string filename;
    uint32_t taxid;
    istringstream iss(line);
    iss >> filename;
    iss >> taxid;
    process_file(filename, taxid);
    cerr << "\rProcessed " << ++seqs_processed << " sequences";
  }
  cerr << "\r                                                       ";
  cerr << "\rFinished processing " << seqs_processed << " sequences" << endl;
}

void process_file(string filename, uint32_t taxid) {
  FastaReader reader(filename);
  DNASequence dna;
  
  // For the purposes of this program, we assume these files are
  // single-fasta files.
  dna = reader.next_sequence();

  for (size_t i = 0; i < dna.seq.size(); i += SKIP_LEN)
    set_lcas(taxid, dna.seq, i, i + SKIP_LEN + Database.get_k() - 1);
}

void set_lcas(uint32_t taxid, string &seq, size_t start, size_t finish) {
  KmerScanner scanner(seq, start, finish);
  uint64_t *kmer_ptr;
  uint32_t *val_ptr;
  int64_t pos_kmer;
  uint64_t canonical_kmer;
  uint32_t parent_taxid;
  std::pair<uint32_t *, int64_t> ret_pair;
  std::set<uint64_t>::iterator set_it;

  while ((kmer_ptr = scanner.next_kmer()) != NULL) {
    if (scanner.ambig_kmer())
      continue;
    canonical_kmer = Database.canonical_representation(*kmer_ptr);
    ret_pair = Database.kmer_query_pair(canonical_kmer);
    val_ptr = ret_pair.first;
    if (val_ptr == NULL) {
      if (! Allow_extra_kmers)
        errx(EX_DATAERR, "kmer found in sequence that is not in database");
      else
        continue;
    }
    pos_kmer = ret_pair.second;
    set_it = kmer_check_set.find(canonical_kmer);
    if (set_it == kmer_check_set.end()) {
      pos_kmer = ret_pair.second;
      kmer_weight[pos_kmer]++;
      kmer_check_set.insert(canonical_kmer);
    }
    
    taxid = at_rank(taxid, rank_level);
    if ( Pruned_parent_map.find(taxid) == Pruned_parent_map.end() ) {
      parent_taxid = Parent_map[taxid];
      Pruned_parent_map[taxid] = parent_taxid;
      while (parent_taxid != 0) {
        Pruned_parent_map.insert ( 
          std::pair<uint32_t,uint32_t>(parent_taxid, Parent_map[parent_taxid] ));
        parent_taxid = Parent_map[parent_taxid];
      }
    }
  }
}

uint32_t count_children (uint32_t taxid) {
  //Leaf Node
  std::multimap<uint32_t,uint32_t>::iterator it = Children_multimap.find(taxid);
  if (it == Children_multimap.end()) {
    Pruned_parent_map[taxid] = 1;
    return 1;
  }
  std::pair <std::multimap<uint32_t,uint32_t>::iterator, 
             std::multimap<uint32_t,uint32_t>::iterator> ret;
  ret = Children_multimap.equal_range(taxid);
  uint32_t count = 0;
  for (it=ret.first; it!=ret.second; ++it)    
    count += count_children(it->second);
  Pruned_parent_map[taxid] = count;
  return count;
}

void calculate_percentage() {
  uint64_t key_len = Database.get_key_len();
  uint64_t val_len = Database.get_val_len();
  uint64_t key_ct = Database.get_key_ct();
  uint64_t pair_size = key_len + val_len;
  char pair[pair_size];
  char *data = new char[ key_ct * (key_len + val_len) ];

  float kmer_perc;

  ifstream input_file(DB_filename.c_str(), std::ifstream::in | std::ifstream::binary);
  input_file.seekg(Database.header_size(), ios_base::beg);
  QuickFile db_file(DB_filename, "r");

  for (uint64_t i = 0; i < key_ct; i++) {
    input_file.read(pair, pair_size);
    uint32_t val = 0;
    memcpy(&val, pair + key_len, val_len);
    kmer_perc = round((float)kmer_weight[i] * 255 / (float)Pruned_parent_map[val]);
    if (kmer_perc > 255)
      kmer_perc = 255;

    if (kmer_perc >= perc_score) {
      char *pair_pos = data + pair_size * (i - del_cnt);
      memcpy(pair_pos, pair, pair_size);
    }
    else 
      del_cnt++;
  }

  if (del_cnt > 0) {
    std::cout << "Filtering DB..." << std::endl;
    ofstream output_file(Output_DB_filename.c_str(), std::ofstream::binary);
    size_t db_header_size = 72 + 2 * (4 + Database.get_key_bits() * 8);
    char *header = new char[ db_header_size ];
    memcpy(header, db_file.ptr(), db_header_size);
    output_file.write(header, db_header_size);
    output_file.write(data, (key_ct - del_cnt) * (key_len + val_len));
    output_file.close();

    uint64_t disc_cnt = key_ct - del_cnt;
    ofstream mod_output_file_size(Output_DB_filename.c_str(), std::ofstream::binary | std::ios::out | std::ios::in);
    mod_output_file_size.seekp(48, ios::beg);
    mod_output_file_size.write((char*)&disc_cnt, sizeof(del_cnt));
    mod_output_file_size.close();

    QuickFile input_reduced_db_file(Output_DB_filename.c_str());
    KrakenDB *input_reduced_db = new KrakenDB(input_reduced_db_file.ptr());
    input_reduced_db->make_index(Index_filename_reduced.c_str(), Bin_key_nt);

    std::cout << "K-mers eliminated: " << del_cnt << ", " << (float)del_cnt * 100 / (float)key_ct << "% of full DB."<< std::endl;
  }
  else
    std::cout << "No filtering needed." << std::endl;
  input_file.close();
}

uint32_t at_rank(uint32_t taxid, uint8_t rank) {
  int last_rank = Rank_map[taxid];
  int last_taxid = taxid;

  while (Rank_map[taxid] > rank) {
    taxid = Parent_map[taxid];
    if (Rank_map[taxid] < last_rank) {
      last_rank = Rank_map[taxid];
      last_taxid = taxid;
    }
  }
  return last_taxid;
}

void parse_command_line(int argc, char **argv) {
  int opt;
  float sfg;
  long long sig;

  if (argc > 1 && strcmp(argv[1], "-h") == 0)
    usage(0);
  while ((opt = getopt(argc, argv, "f:d:i:n:o:m:r:l:F:S:xM")) != -1) {
    switch (opt) {
      case 'f' :
        File_to_taxon_map_filename = optarg;
        break;
      case 'd' :
        DB_filename = optarg;
        break;
      case 'i' :
        Index_filename = optarg;
        break;
      case 'F' :
        Multi_fasta_filename = optarg;
        break;
      case 'm' :
        ID_to_taxon_map_filename = optarg;
        break;
      case 'n' :
        Nodes_filename = optarg;
        break;
      case 'o' :
        Output_DB_filename = optarg;
        break;
      case 'r' :
        Index_filename_reduced = optarg;
        break;
      case 'l' :
        sig = atoll(optarg);
        if (sig < 1 || sig > 31)
          errx(EX_USAGE, "bin key length out of range");
        Bin_key_nt = (uint8_t) sig;
        break;
      case 'x' :
        Allow_extra_kmers = true;
        break;
      case 'M' :
        Operate_in_RAM = true;
        break;
      case 'S' :
        sfg = atof(optarg);
        if (sfg < 0 || sfg > 100)
          errx(EX_USAGE, "Score must be in the [0, 100] interval");
        perc_score = round((sfg)*2.55);
        break;
      default:
        usage();
        break;
    }
  }

  if (DB_filename.empty() || Index_filename.empty() ||
      Nodes_filename.empty())
    usage();
  if (File_to_taxon_map_filename.empty() &&
      (Multi_fasta_filename.empty() || ID_to_taxon_map_filename.empty()))
    usage();

  if (! File_to_taxon_map_filename.empty())
    One_FASTA_file = false;
  else
    One_FASTA_file = true;
}

void usage(int exit_code) {
  cerr << "Usage: set_lcas [options]" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -d filename      Kraken DB filename" << endl
       << "* -i filename      Kraken DB index filename" << endl
       << "* -n filename      NCBI Taxonomy nodes file" << endl
       << "* -o filename      Output Kraken DB filename" << endl
       << "* -r filename      Output Kraken DB index filename" << endl
       << "* -l NUM           Minimizer length" << endl
       << "* -S NUM           Score filter parameter" << endl
       << "  -x               K-mers not found in DB do not cause errors" << endl
       << "  -f filename      File to taxon map" << endl
       << "  -F filename      Multi-FASTA file with sequence data" << endl
       << "  -m filename      Sequence ID to taxon map" << endl
       << "  -h               Print this message" << endl
       << endl
       << "-F and -m must be specified together.  If -f is given, "
       << "-F/-m are ignored." << endl;
  exit(exit_code);
}
