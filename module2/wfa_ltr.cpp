#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include <unistd.h>   // getopt
#include <iomanip>
#include <cmath>
#include <cctype>
#include <algorithm>
#include "bindings/cpp/WFAligner.hpp"
#include <limits>

// g++ -O3 -std=c++11 module2/wfa_ltr.cpp ./WFA2-lib/bindings/cpp/WFAligner.cpp    -I./WFA2-lib -I./WFA2-lib/bindings/cpp     ./WFA2-lib/build/libwfa2.a -o module2/wfa_ltr

using namespace std;
using namespace wfa;

struct Candidate {
    // SCN fields we care about (1-based, inclusive)
    string chr;
    long s_lLTR;
    long e_lLTR;
    long s_rLTR;
    long e_rLTR;
};

// -------------------------- Utilities --------------------------
static inline string trim(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b-a+1);
}

void split_ws(const string& line, vector<string>& out) {
    out.clear();
    string tok;
    istringstream iss(line);
    while (iss >> tok) out.push_back(tok);
}

void print_usage(const char* prog) {
    cerr << "Usage: " << prog << " [options] <SCN> <genome.fa>\n"
         << "Options:\n"
         << "  -x <int>        mismatch penalty (default 6)\n"
         << "  -O <int>        gap opening1 penalty (default 12)\n"
         << "  -E <int>        gap extension1 penalty (default 4)\n"
         << "  -o <int>        gap opening2 penalty (default 80)\n"
         << "  -e <int>        gap extension2 penalty (default 2)\n"
         << "  -c [sam|full]   CIGAR format printed in table (default full)\n"
         << "  -f <int>        flank length around s/e positions (default 70)\n"
         << "  -A              also print basepair alignment under each row (off)\n"
         << "  -w <int>        wrap width for basepair alignment (default 200)\n"
         << "  --match-thresh <int>   threshold for long '=' cluster (default 14)\n"
         << "  --allowed-gaps  <int>  how many 'I'/'D' events can interrupt (default 0)\n"
         << "  --allowed-subs  <int>  how many 'X' bases can interrupt     (default 2)\n"
         << "  --vic-in <int>         bases to include inside boundary  (default 5)\n"
         << "  --vic-out <int>        bases to include outside boundary (default 10)\n"
         << "  --win-pairs <int>      window size (aligned pairs) for boundary QC (default 20)\n"
         << "  --thresh-high <dbl>    high identity threshold; checks '>' (default 0.70)\n"
         << "  --thresh-low  <dbl>    low  identity threshold; checks '<' (default 0.70)\n"
         << "  -h              show this help message\n\n"
         << "Output columns:\n"
         << "  header1, header2, CIGAR, [identity if full],\n"
         << "  total_length, total_substitutions, transitions, transversions,\n"
         << "  raw_d, JC69_d, K2P_d\n"
         << "Plus, for each candidate, one line:\n"
         << "  # TRUE_BOUNDARY <chr> <true lLTR start> <true lLTR end> <true rLTR start> <true rLTR end> <TSD> <motif>\n";
}

// Read a multi-FASTA into map<id, seq> (ID is up to first whitespace)
bool read_fasta(const string& path, unordered_map<string,string>& seqs, unordered_map<string,size_t>& lens) {
    ifstream in(path);
    if (!in) return false;
    string line, id;
    ostringstream seq;
    auto flush = [&](bool force=false){
        if (!id.empty() && (force || !seq.str().empty())) {
            string s = seq.str();
            for (char& c : s) c = toupper(c);
            seqs[id] = move(s);
            lens[id] = seqs[id].size();
            seq.str(""); seq.clear();
        }
    };
    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            flush();
            string hdr = line.substr(1);
            size_t sp = hdr.find_first_of(" \t");
            id = (sp==string::npos)? hdr : hdr.substr(0, sp);
        } else {
            seq << line;
        }
    }
    flush(true);
    return !seqs.empty();
}

string make_header(const string& chr, long start, long end) {
    ostringstream oss;
    oss << chr << ":" << start << "-" << end;
    return oss.str();
}

// Extract [start,end] 1-based, inclusive, clamped; return empty on invalid
string extract_region(const string& chrom_seq, long start1, long end1) {
    if (start1 > end1) return string();
    long n = static_cast<long>(chrom_seq.size());
    long s = max(1L, min(start1, n));
    long e = max(1L, min(end1,   n));
    if (s > e) return string();
    long off = s - 1;
    long len = e - s + 1;
    return chrom_seq.substr(static_cast<size_t>(off), static_cast<size_t>(len));
}

double compute_gc_identity(const string& cigar) {
    long matches = 0, mismatches = 0, gaps = 0;
    int num = 0;
    for (char c : cigar) {
        if (isdigit(c)) {
            num = num * 10 + (c - '0');
        } else {
            switch (c) {
                case '=': matches    += num; break;
                case 'X': mismatches += num; break;
                case 'D': case 'I': gaps += 1; break;
                default: break;
            }
            num = 0;
        }
    }
    long denom = matches + mismatches + gaps;
    return (denom>0)? double(matches)/double(denom) : 0.0;
}

void count_substitutions(
    const string& cigar,
    const string& s1,
    const string& s2,
    long &total_subs,
    long &transitions,
    long &transversions
) {
    size_t i = 0, j = 0, pos = 0;
    total_subs = transitions = transversions = 0;
    while (pos < cigar.size()) {
        int num = 0;
        while (pos < cigar.size() && isdigit(cigar[pos])) {
            num = num * 10 + (cigar[pos++] - '0');
        }
        if (pos >= cigar.size()) break;
        char op = cigar[pos++];
        if (op == '=') {
            i += num; j += num;
        } else if (op == 'X') {
            for (int k = 0; k < num; ++k) {
                if (i >= s1.size() || j >= s2.size()) break;
                char a = toupper(s1[i++]);
                char b = toupper(s2[j++]);
                ++total_subs;
                if ((a=='A' && b=='G') || (a=='G' && b=='A') ||
                    (a=='C' && b=='T') || (a=='T' && b=='C')) {
                    ++transitions;
                } else {
                    ++transversions;
                }
            }
        } else if (op == 'D') {
            i += num;
        } else if (op == 'I') {
            j += num;
        }
    }
}

void build_pretty_alignment(
    const string& cigar_full,
    const string& s1,
    const string& s2,
    string& out_top,
    string& out_mid,
    string& out_bot
) {
    out_top.clear(); out_mid.clear(); out_bot.clear();
    size_t i = 0, j = 0, pos = 0;
    while (pos < cigar_full.size()) {
        int num = 0;
        while (pos < cigar_full.size() && isdigit(cigar_full[pos])) {
            num = num * 10 + (cigar_full[pos++] - '0');
        }
        if (pos >= cigar_full.size()) break;
        char op = cigar_full[pos++];
        if (op == '=') {
            for (int k = 0; k < num; ++k) {
                char a = (i<s1.size()) ? s1[i++] : '-';
                char b = (j<s2.size()) ? s2[j++] : '-';
                out_top.push_back(a);
                out_mid.push_back('|');
                out_bot.push_back(b);
            }
        } else if (op == 'X') {
            for (int k = 0; k < num; ++k) {
                char a = (i<s1.size()) ? s1[i++] : '-';
                char b = (j<s2.size()) ? s2[j++] : '-';
                out_top.push_back(a);
                out_mid.push_back(' ');
                out_bot.push_back(b);
            }
        } else if (op == 'D') {
            for (int k = 0; k < num; ++k) {
                char a = (i<s1.size()) ? s1[i++] : '-';
                out_top.push_back(a);
                out_mid.push_back(' ');
                out_bot.push_back('-');
            }
        } else if (op == 'I') {
            for (int k = 0; k < num; ++k) {
                char b = (j<s2.size()) ? s2[j++] : '-';
                out_top.push_back('-');
                out_mid.push_back(' ');
                out_bot.push_back(b);
            }
        }
    }
}

void print_pretty_alignment_block(
    const string& h1, const string& h2,
    const string& top, const string& mid, const string& bot,
    int width
) {
    cout << "# ALIGNMENT " << h1 << "  vs  " << h2 << "\n";
    for (size_t p = 0; p < top.size(); p += width) {
        size_t n = min<size_t>(width, top.size() - p);
        cout << top.substr(p, n) << "\n";
        cout << mid.substr(p, n) << "\n";
        cout << bot.substr(p, n) << "\n\n";
    }
}

// ------- CIGAR helpers to find boundary & convert to sequence offsets -------
struct CigarOp { char op; int len; };

// parse full CIGAR "=\0/IXD..."
static inline vector<CigarOp> parse_cigar(const string& cigar_full) {
    vector<CigarOp> ops;
    int num = 0;
    for (char c : cigar_full) {
        if (isdigit(c)) {
            num = num*10 + (c-'0');
        } else {
            ops.push_back({c, num});
            num = 0;
        }
    }
    return ops;
}

// Compute prefix offsets (i,j) at op boundaries.
// start_off[i] = offsets BEFORE op i; end_off[i] = offsets AFTER op i.
static inline void compute_prefix_offsets(
    const vector<CigarOp>& ops,
    vector<pair<size_t,size_t>>& start_off,
    vector<pair<size_t,size_t>>& end_off
) {
    start_off.resize(ops.size());
    end_off.resize(ops.size());
    size_t i = 0, j = 0;
    for (size_t t = 0; t < ops.size(); ++t) {
        start_off[t] = {i,j};
        char op = ops[t].op;
        int  L  = ops[t].len;
        if (op=='=' || op=='X') { i += (size_t)L; j += (size_t)L; }
        else if (op=='D') { i += (size_t)L; }
        else if (op=='I') { j += (size_t)L; }
        end_off[t] = {i,j};
    }
}

// Build a cluster of '=' blocks anchored at index 'anchor' (which must be '='),
// extending in 'dir' (+1 for right, -1 for left). Disruptions (I/D/X) are allowed
// ONLY when followed by another '=' in the same direction; thus they appear
// strictly BETWEEN '=' blocks. Returns the [L,R] op-index span of the cluster,
// plus aggregated stats.
struct Cluster {
    int L=-1, R=-1;          // op indices inclusive
    long matches=0;          // sum of '=' lengths
    int gap_events=0;        // count of I/D ops (length ignored)
    long subs_events=0;      // sum of X lengths
    bool valid=false;
};

static inline Cluster build_cluster_from_anchor(
    const vector<CigarOp>& ops, int anchor, int dir,
    int allowed_gaps, int allowed_subs
) {
    Cluster c;
    if (anchor < 0 || anchor >= (int)ops.size() || ops[anchor].op!='=') return c;

    c.L = c.R = anchor;
    c.matches = ops[anchor].len;
    c.gap_events = 0;
    c.subs_events = 0;
    c.valid = true;

    // Helper to attempt to include an adjacent '=' (free) or
    // a (disruption + next '=') pair (cost), repeatedly.
    auto step_once = [&](int cur_edge)->bool{
        int next = cur_edge + dir;
        if (next < 0 || next >= (int)ops.size()) return false;
        if (ops[next].op=='=') {
            // add contiguous '=' without cost
            if (dir>0) c.R = next; else c.L = next;
            c.matches += ops[next].len;
            return true;
        }
        // disruption must be followed by another '=' in same direction
        if (ops[next].op=='X' || ops[next].op=='I' || ops[next].op=='D') {
            int next2 = next + dir;
            if (next2 < 0 || next2 >= (int)ops.size()) return false;
            if (ops[next2].op!='=') return false;
            // tentative new tallies
            int add_gap = (ops[next].op=='I' || ops[next].op=='D') ? 1 : 0;
            long add_sub = (ops[next].op=='X') ? ops[next].len : 0;
            int new_gaps = c.gap_events + add_gap;
            long new_subs = c.subs_events + add_sub;
            if (new_gaps > allowed_gaps || new_subs > allowed_subs) return false;
            // include disruption + following '='
            c.gap_events = new_gaps;
            c.subs_events = new_subs;
            if (dir>0) { c.R = next2; } else { c.L = next2; }
            c.matches += ops[next2].len;
            return true;
        }
        return false;
    };

    // Step outward until we can no longer include pieces under allowances
    while (step_once(dir>0 ? c.R : c.L)) {}
    return c;
}

// Sliding search for relaxed cluster meeting thresholds.
// For start boundary (right->left): pick FIRST cluster from right whose matches >= thresh.
// Return LEFT edge offsets (0-based) in (s1,s2) at the start of the cluster.
// For end boundary (left->right): pick LAST cluster from left; return RIGHT edge offsets.
static inline bool relaxed_boundary_for_start_pair(
    const string& cigar_full, int match_thresh,
    int allowed_gaps, int allowed_subs,
    size_t& i_off, size_t& j_off
) {
    auto ops = parse_cigar(cigar_full);
    if (ops.empty()) return false;

    vector<pair<size_t,size_t>> start_off, end_off;
    compute_prefix_offsets(ops, start_off, end_off);

    for (int t = (int)ops.size()-1; t >= 0; --t) {
        if (ops[t].op != '=') continue;
        Cluster cl = build_cluster_from_anchor(ops, t, -1, allowed_gaps, allowed_subs);
        if (cl.valid && cl.matches >= match_thresh) {
            // left edge at start of op cl.L
            i_off = start_off[cl.L].first;
            j_off = start_off[cl.L].second;
            return true;
        }
    }
    return false;
}

static inline bool relaxed_boundary_for_end_pair(
    const string& cigar_full, int match_thresh,
    int allowed_gaps, int allowed_subs,
    size_t& i_off_after, size_t& j_off_after
) {
    auto ops = parse_cigar(cigar_full);
    if (ops.empty()) return false;

    vector<pair<size_t,size_t>> start_off, end_off;
    compute_prefix_offsets(ops, start_off, end_off);

    bool found = false;
    size_t bi=0, bj=0;
    for (int t = 0; t < (int)ops.size(); ++t) {
        if (ops[t].op != '=') continue;
        Cluster cl = build_cluster_from_anchor(ops, t, +1, allowed_gaps, allowed_subs);
        if (cl.valid && cl.matches >= match_thresh) {
            // right edge AFTER op cl.R
            bi = end_off[cl.R].first;
            bj = end_off[cl.R].second;
            found = true; // keep the last one
            // continue to seek later clusters; this ensures "LAST"
        }
    }
    if (found) {
        i_off_after = bi;
        j_off_after = bj;
    }
    return found;
}

// Clamp helper (1-based)
static inline long clamp1(long x, long lo, long hi) {
    return max(lo, min(x, hi));
}

// Extract 1-based inclusive [L,R] from chromosome (empty if invalid)
static inline string chr_slice_1based(const string& chrseq, long L, long R) {
    if (L>R) return string();
    long n = (long)chrseq.size();
    L = clamp1(L, 1, n);
    R = clamp1(R, 1, n);
    if (L>R) return string();
    return chrseq.substr((size_t)(L-1), (size_t)(R-L+1));
}

// Try to find a TSD given the two boundary coordinates and vic windows.
static inline string detect_tsd(
    const string& chrseq,
    long left_start,   // true lLTR start
    long right_end,    // true rLTR end
    int vic_in,
    int vic_out
) {
    if (left_start<=0 || right_end<=0) return "NA";
    long n = (long)chrseq.size();

    // TSD length limits (explicit)
    const int MIN_TSD = 5;
    const int MAX_TSD = 8;

    // windows around both boundaries
    long L1 = clamp1(left_start - vic_out, 1, n);
    long R1 = clamp1(left_start + vic_in - 1, 1, n);
    string win5 = chr_slice_1based(chrseq, L1, R1);

    long L2 = clamp1(right_end - (vic_in - 1), 1, n);
    long R2 = clamp1(right_end + vic_out, 1, n);
    string win3 = chr_slice_1based(chrseq, L2, R2);

    if (win5.empty() || win3.empty()) return "NA";

    int len5 = (int)win5.size();
    int len3 = (int)win3.size();

    int out5 = min(vic_out, len5);
    int out3 = min(vic_out, len3);
    int in5  = min(vic_in,  len5);
    int in3  = min(vic_in,  len3);

    for (int k = min({MAX_TSD, out5, out3}); k >= MIN_TSD; --k) {
        if (out5 < k || out3 < k) break;
        string s5 = win5.substr(out5 - k, k);
        string s3 = win3.substr(in3, k);
        if (s5 == s3) return s5;
    }

    int maxL = min(MAX_TSD, min(out5 + in5, len5));
    int maxR = min(MAX_TSD, min(out3 + in3, len3));
    for (int k = min(maxL, maxR); k >= MIN_TSD; --k) {
        int max_take5_out = min(out5, k-1);
        for (int take5_out = max_take5_out; take5_out >= 1; --take5_out) {
            int take5_in = k - take5_out;
            if (take5_in > in5) continue;
            int start5 = out5 - take5_out;
            string s5 = win5.substr(start5, k);

            int max_take3_out = min(out3, k-1);
            for (int take3_out = max_take3_out; take3_out >= 1; --take3_out) {
                int take3_in = k - take3_out;
                if (take3_in > in3) continue;
                int start3 = in3 - take3_in;
                string s3 = win3.substr(start3, k);
                if (s5 == s3) return s5;
            }
        }
    }

    return "NA";
}

// Extract observed terminal dinucleotides at the called boundaries.
static inline string extract_terminal_motif(
    const string& chrseq,
    long true_l_start,  // first base of element (1-based)
    long true_r_end     // last  base of element (1-based)
) {
    if (true_l_start <= 0 || true_r_end <= 0) return "NA";
    long n = (long)chrseq.size();
    if (true_l_start + 1 > n || true_r_end - 1 < 1) return "NA";
    string left2  = chr_slice_1based(chrseq, true_l_start, true_l_start + 1);
    string right2 = chr_slice_1based(chrseq, true_r_end - 1, true_r_end);
    for (char& c : left2)  c = toupper(c);
    for (char& c : right2) c = toupper(c);
    return left2 + "…" + right2;
}

// Compute identity over first/last K aligned pairs (only '='/'X' count).
// from_left=true -> first K; false -> last K.
// Returns NaN if there are 0 aligned pairs in that window.
static inline double window_identity_aligned_pairs(const string& cigar_full, int K, bool from_left) {
    if (K <= 0) return numeric_limits<double>::quiet_NaN();
    vector<CigarOp> ops = parse_cigar(cigar_full);
    long aligned = 0;
    long matches = 0;

    if (from_left) {
        for (const auto& op : ops) {
            if (op.op=='=' || op.op=='X') {
                int take = min<int>(op.len, (int)(K - aligned));
                if (take > 0) {
                    if (op.op=='=') matches += take;
                    aligned += take;
                    if (aligned >= K) break;
                }
            }
            // skip I/D columns for the window definition
        }
    } else {
        for (int t=(int)ops.size()-1; t>=0 && aligned<K; --t) {
            const auto& op = ops[t];
            if (op.op=='=' || op.op=='X') {
                int take = min<int>(op.len, (int)(K - aligned));
                if (take > 0) {
                    if (op.op=='=') matches += take;
                    aligned += take;
                }
            }
        }
    }

    if (aligned == 0) return numeric_limits<double>::quiet_NaN();
    return double(matches) / double(aligned);
}



// Align two sequences and print one output row (+ optional pretty alignment)
void do_alignment_row(
    WFAlignerGapAffine2Pieces& aligner,
    const string& s1, const string& s2,
    const string& h1, const string& h2,
    const string& cigar_print_format,  // "full" or "sam" for table
    bool print_alignment_block,
    int wrap_width
) {
    int status = aligner.alignEnd2End(s1, s2);
    if (status < 0) {
        cerr << "Alignment failed (" << h1 << " vs " << h2 << "), status=" << status << "\n";
        return;
    }

    const string cigar_full  = aligner.getCIGAR(true);
    const string cigar_table = (cigar_print_format=="full") ? cigar_full : aligner.getCIGAR(false);

    long total_subs=0, trans_count=0, transv_count=0;
    count_substitutions(cigar_full, s1, s2, total_subs, trans_count, transv_count);

    double tot_len = double(s1.size());
    double raw_d  = (tot_len>0.0)? (trans_count + transv_count) / tot_len : 0.0;
    double P      = (tot_len>0.0)? trans_count / tot_len : 0.0;
    double Q      = (tot_len>0.0)? transv_count / tot_len : 0.0;

    double JC69_d = (1.0 - (4.0 * raw_d / 3.0) > 0.0) ? (-0.75 * log(1.0 - (4.0*raw_d/3.0))) : numeric_limits<double>::quiet_NaN();
    double inner  = (1.0 - 2.0*P - Q);
    double inner2 = (1.0 - 2.0*Q);
    double K2P_d  = (inner>0.0 && inner2>0.0) ? (-0.5 * log(inner * sqrt(inner2))) : numeric_limits<double>::quiet_NaN();

    cout << h1 << "\t" << h2 << "\t" << cigar_table;
    if (cigar_print_format=="full") {
        double id = compute_gc_identity(cigar_full);
        cout << "\t" << fixed << setprecision(6) << id;
    }
    cout << "\t"
         << long(tot_len) << "\t"
         << total_subs << "\t"
         << trans_count << "\t"
         << transv_count << "\t"
         << fixed << setprecision(6)
         << raw_d    << "\t"
         << JC69_d   << "\t"
         << K2P_d
         << "\n";

    if (print_alignment_block) {
        string top, mid, bot;
        build_pretty_alignment(cigar_full, s1, s2, top, mid, bot);
        print_pretty_alignment_block(h1, h2, top, mid, bot, wrap_width);
    }
}

// ---------------------------- main ----------------------------
int main(int argc, char** argv) {
    // defaults (updated)
    int mismatch  = 6;
    int go1       = 12;
    int ge1       = 4;
    int go2       = 80;
    int ge2       = 2;
    string format = "full";
    long flank    = 70;
    bool print_alignment_block = false;
    int wrap_width = 200;

    int match_thresh = 14;
    int vic_in = 5;
    int vic_out = 10;

    int allowed_gaps = 0;
    int allowed_subs = 2;

    // NEW: boundary QC window & thresholds
    int    win_pairs    = 20;   // replaces WIN
    double thresh_high  = 0.70; // replaces hard-coded >0.70 checks
    double thresh_low   = 0.70; // replaces hard-coded <0.70 checks

    // getopt long-ish handling
    std::vector<std::string> argv_vec;
    argv_vec.reserve(argc);
    for (int i = 0; i < argc; ++i) argv_vec.emplace_back(argv[i]);

    // Build a new argv with long options stripped, but KEEP STORAGE ALIVE
    std::vector<std::string> keep;
    keep.reserve(argv_vec.size());
    keep.push_back(argv_vec[0]);

    for (size_t i = 1; i < argv_vec.size(); ++i) {
        if (argv_vec[i] == "--match-thresh" && i + 1 < argv_vec.size()) {
            match_thresh = std::stoi(argv_vec[i + 1]); ++i;
        } else if (argv_vec[i] == "--allowed-gaps" && i + 1 < argv_vec.size()) {
            allowed_gaps = std::stoi(argv_vec[i + 1]); ++i;
        } else if (argv_vec[i] == "--allowed-subs" && i + 1 < argv_vec.size()) {
            allowed_subs = std::stoi(argv_vec[i + 1]); ++i;
        } else if (argv_vec[i] == "--vic-in" && i + 1 < argv_vec.size()) {
            vic_in = std::stoi(argv_vec[i + 1]); ++i;
        } else if (argv_vec[i] == "--vic-out" && i + 1 < argv_vec.size()) {
           vic_out = std::stoi(argv_vec[i + 1]); ++i;
        } else if (argv_vec[i] == "--win-pairs" && i + 1 < argv_vec.size()) {
            win_pairs = std::stoi(argv_vec[i + 1]); ++i;
        } else if (argv_vec[i] == "--thresh-high" && i + 1 < argv_vec.size()) {
            thresh_high = std::stod(argv_vec[i + 1]); ++i;
        } else if (argv_vec[i] == "--thresh-low" && i + 1 < argv_vec.size()) {
            thresh_low = std::stod(argv_vec[i + 1]); ++i;
        } else {
            keep.push_back(argv_vec[i]);
        }
    }

    // Create a stable argv array backed by `keep`
    std::vector<char*> argv_stable;
    argv_stable.reserve(keep.size());
    for (auto &s : keep) argv_stable.push_back(const_cast<char*>(s.c_str()));
    // No need for trailing nullptr for getopt in C++; we’ll pass argc/argv explicitly.

    argc = static_cast<int>(argv_stable.size());
    argv = argv_stable.data();

    // now run getopt as usual
    int opt;
    while ((opt = getopt(argc, argv, "x:O:E:o:e:c:f:Aw:h")) != -1) {

        switch (opt) {
            case 'x': mismatch = stoi(optarg);   break;
            case 'O': go1      = stoi(optarg);   break;
            case 'E': ge1      = stoi(optarg);   break;
            case 'o': go2      = stoi(optarg);   break;
            case 'e': ge2      = stoi(optarg);   break;
            case 'c': format   = optarg;         break;
            case 'f': flank    = stol(optarg);   break;
            case 'A': print_alignment_block = true; break;
            case 'w': wrap_width = max(1, stoi(optarg)); break;
            case 'h':
            default:
                print_usage(argv[0]);
                return (opt=='h'? 0 : 1);
        }
    }

    if (format!="full" && format!="sam") {
        cerr << "Error: unknown CIGAR format '" << format << "'\n";
        print_usage(argv[0]);
        return 1;
    }
    if (optind + 2 != argc) {
        print_usage(argv[0]);
        return 1;
    }
    const string scn_path   = argv[optind];
    const string genome_path= argv[optind+1];

    // load genome
    unordered_map<string,string> genome;
    unordered_map<string,size_t> chrom_lens;
    if (!read_fasta(genome_path, genome, chrom_lens)) {
        cerr << "Error: failed to read genome FASTA '" << genome_path << "'\n";
        return 1;
    }

    // read SCN
    auto read_scn = [](const string& path, vector<Candidate>& out)->bool{
        ifstream in(path);
        if (!in) return false;
        string line;
        while (getline(in, line)) {
            line = trim(line);
            if (line.empty() || line[0]=='#') continue;
            vector<string> f; split_ws(line, f);
            if (f.size() < 12) continue;
            Candidate c;
            c.s_lLTR = stol(f[3]);
            c.e_lLTR = stol(f[4]);
            c.s_rLTR = stol(f[6]);
            c.e_rLTR = stol(f[7]);
            c.chr    = f[11];
            out.push_back(move(c));
        }
        return true;
    };

    vector<Candidate> cand;
    if (!read_scn(scn_path, cand) || cand.empty()) {
        cerr << "Error: no candidates parsed from SCN '" << scn_path << "'\n";
        return 1;
    }

    // prepare aligner
    WFAlignerGapAffine2Pieces aligner(
        mismatch, go1, ge1, go2, ge2,
        WFAligner::Alignment, WFAligner::MemoryHigh
    );

    // process each candidate -> two alignments + true-boundary call
    for (const auto& c : cand) {
        auto it = genome.find(c.chr);
        if (it == genome.end()) {
            cerr << "Warning: chromosome '" << c.chr << "' not found; skipping\n";
            continue;
        }
        const string& chrom_seq = it->second;
        long chrN = (long)chrom_seq.size();

        // starts
        long a1 = c.s_lLTR - flank;
        long a2 = c.s_lLTR + flank;
        long b1 = c.s_rLTR - flank;
        long b2 = c.s_rLTR + flank;

        string s_left  = extract_region(chrom_seq, a1, a2);
        string s_right = extract_region(chrom_seq, b1, b2);

        string cigar_start_pair;
        size_t s_i = 0, s_j = 0; // boundary offsets for starts

        if (!s_left.empty() && !s_right.empty()) {
            string h1 = make_header(c.chr, max(1L,a1), min<long>(a2, chrN));
            string h2 = make_header(c.chr, max(1L,b1), min<long>(b2, chrN));

            int status = aligner.alignEnd2End(s_left, s_right);
            if (status < 0) {
                cerr << "Alignment failed (" << h1 << " vs " << h2 << "), status=" << status << "\n";
            } else {
                const string cigar_full  = aligner.getCIGAR(true);
                cigar_start_pair = cigar_full;
                // print the normal row
                {
                    const string cigar_table = (format=="full") ? cigar_full : aligner.getCIGAR(false);
                    long total_subs=0, trans_count=0, transv_count=0;
                    count_substitutions(cigar_full, s_left, s_right, total_subs, trans_count, transv_count);
                    double tot_len = double(s_left.size());
                    double raw_d  = (tot_len>0.0)? (trans_count + transv_count) / tot_len : 0.0;
                    double P      = (tot_len>0.0)? trans_count / tot_len : 0.0;
                    double Q      = (tot_len>0.0)? transv_count / tot_len : 0.0;
                    double JC69_d = (1.0 - (4.0 * raw_d / 3.0) > 0.0) ? (-0.75 * log(1.0 - (4.0*raw_d/3.0))) : numeric_limits<double>::quiet_NaN();
                    double inner  = (1.0 - 2.0*P - Q);
                    double inner2 = (1.0 - 2.0*Q);
                    double K2P_d  = (inner>0.0 && inner2>0.0) ? (-0.5 * log(inner * sqrt(inner2))) : numeric_limits<double>::quiet_NaN();

                    cout << h1 << "\t" << h2 << "\t" << cigar_table;
                    if (format=="full") {
                        double id = compute_gc_identity(cigar_full);
                        cout << "\t" << fixed << setprecision(6) << id;
                    }
                    cout << "\t" << long(tot_len) << "\t" << total_subs << "\t" << trans_count << "\t" << transv_count
                         << "\t" << fixed << setprecision(6) << raw_d << "\t" << JC69_d << "\t" << K2P_d << "\n";
                    if (print_alignment_block) {
                        string top, mid, bot;
                        build_pretty_alignment(cigar_full, s_left, s_right, top, mid, bot);
                        print_pretty_alignment_block(h1, h2, top, mid, bot, wrap_width);
                    }
                }
                // relaxed boundary for start pair
                if (!relaxed_boundary_for_start_pair(cigar_full, match_thresh, allowed_gaps, allowed_subs, s_i, s_j)) {
                    cerr << "Warning: start-pair boundary not found at " << h1 << " vs " << h2
                         << " (no '=' cluster >= " << match_thresh << " within allowances: gaps="
                         << allowed_gaps << ", subs=" << allowed_subs << ")\n";
                }
            }
        } else {
            cerr << "Warning: empty extraction for starts at " << c.chr << " s(lLTR)="
                 << c.s_lLTR << " s(rLTR)=" << c.s_rLTR << "\n";
        }

        // ends
        long a3 = c.e_lLTR - flank;
        long a4 = c.e_lLTR + flank;
        long b3 = c.e_rLTR - flank;
        long b4 = c.e_rLTR + flank;

        string e_left  = extract_region(chrom_seq, a3, a4);
        string e_right = extract_region(chrom_seq, b3, b4);

        string cigar_end_pair;
        size_t e_i_after = 0, e_j_after = 0; // boundary offsets (AFTER block) for ends

        if (!e_left.empty() && !e_right.empty()) {
            string h1 = make_header(c.chr, max(1L,a3), min<long>(a4, chrN));
            string h2 = make_header(c.chr, max(1L,b3), min<long>(b4, chrN));

            int status = aligner.alignEnd2End(e_left, e_right);
            if (status < 0) {
                cerr << "Alignment failed (" << h1 << " vs " << h2 << "), status=" << status << "\n";
            } else {
                const string cigar_full  = aligner.getCIGAR(true);
                cigar_end_pair = cigar_full;
                // print the normal row
                {
                    const string cigar_table = (format=="full") ? cigar_full : aligner.getCIGAR(false);
                    long total_subs=0, trans_count=0, transv_count=0;
                    count_substitutions(cigar_full, e_left, e_right, total_subs, trans_count, transv_count);
                    double tot_len = double(e_left.size());
                    double raw_d  = (tot_len>0.0)? (trans_count + transv_count) / tot_len : 0.0;
                    double P      = (tot_len>0.0)? trans_count / tot_len : 0.0;
                    double Q      = (tot_len>0.0)? transv_count / tot_len : 0.0;
                    double JC69_d = (1.0 - (4.0 * raw_d / 3.0) > 0.0) ? (-0.75 * log(1.0 - (4.0*raw_d/3.0))) : numeric_limits<double>::quiet_NaN();
                    double inner  = (1.0 - 2.0*P - Q);
                    double inner2 = (1.0 - 2.0*Q);
                    double K2P_d  = (inner>0.0 && inner2>0.0) ? (-0.5 * log(inner * sqrt(inner2))) : numeric_limits<double>::quiet_NaN();

                    cout << h1 << "\t" << h2 << "\t" << cigar_table;
                    if (format=="full") {
                        double id = compute_gc_identity(cigar_full);
                        cout << "\t" << fixed << setprecision(6) << id;
                    }
                    cout << "\t" << long(tot_len) << "\t" << total_subs << "\t" << trans_count << "\t" << transv_count
                         << "\t" << fixed << setprecision(6) << raw_d << "\t" << JC69_d << "\t" << K2P_d << "\n";
                    if (print_alignment_block) {
                        string top, mid, bot;
                        build_pretty_alignment(cigar_full, e_left, e_right, top, mid, bot);
                        print_pretty_alignment_block(h1, h2, top, mid, bot, wrap_width);
                    }
                }
                // relaxed boundary for end pair
                if (!relaxed_boundary_for_end_pair(cigar_full, match_thresh, allowed_gaps, allowed_subs, e_i_after, e_j_after)) {
                    cerr << "Warning: end-pair boundary not found at " << h1 << " vs " << h2
                         << " (no '=' cluster >= " << match_thresh << " within allowances: gaps="
                         << allowed_gaps << ", subs=" << allowed_subs << ")\n";
                }
            }
        } else {
            cerr << "Warning: empty extraction for ends at " << c.chr << " e(lLTR)="
                 << c.e_lLTR << " e(rLTR)=" << c.e_rLTR << "\n";
        }

        // If both boundaries exist, convert to chromosome coords & TSD
        long true_l_start = 0, true_r_start = 0;
        if (!s_left.empty() && !s_right.empty() && !cigar_start_pair.empty()) {
            if (s_i || s_j || relaxed_boundary_for_start_pair(cigar_start_pair, match_thresh, allowed_gaps, allowed_subs, s_i, s_j)) {
                long winL_start = max(1L, a1);
                long winR_start = max(1L, b1);
                true_l_start = winL_start + (long)s_i;
                true_r_start = winR_start + (long)s_j;
            }
        }

        long true_l_end = 0, true_r_end = 0;
        if (!e_left.empty() && !e_right.empty() && !cigar_end_pair.empty()) {
            if (e_i_after || e_j_after || relaxed_boundary_for_end_pair(cigar_end_pair, match_thresh, allowed_gaps, allowed_subs, e_i_after, e_j_after)) {
                long winL_end = max(1L, a3);
                long winR_end = max(1L, b3);
                true_l_end = winL_end + (long)e_i_after - 1;
                true_r_end = winR_end + (long)e_j_after - 1;
            }
        }

        string TSD = "NA";
        if (true_l_start>0 && true_r_end>0) {
            TSD = detect_tsd(chrom_seq, true_l_start, true_r_end, vic_in, vic_out);
        }

        string motif = "NA";
        if (true_l_start>0 && true_r_end>0) {
            motif = extract_terminal_motif(chrom_seq, true_l_start, true_r_end);
        }

        // ---------- PASS/FAIL based on parameterized window identities ----------
        auto ok = [](double x){ return !std::isnan(x); };

        double s_first = numeric_limits<double>::quiet_NaN();
        double s_last  = numeric_limits<double>::quiet_NaN();
        double e_first = numeric_limits<double>::quiet_NaN();
        double e_last  = numeric_limits<double>::quiet_NaN();

        if (!cigar_start_pair.empty()) {
            s_first = window_identity_aligned_pairs(cigar_start_pair, win_pairs, /*from_left=*/true);
            s_last  = window_identity_aligned_pairs(cigar_start_pair, win_pairs, /*from_left=*/false);
        }
        if (!cigar_end_pair.empty()) {
            e_first = window_identity_aligned_pairs(cigar_end_pair, win_pairs, /*from_left=*/true);
            e_last  = window_identity_aligned_pairs(cigar_end_pair, win_pairs, /*from_left=*/false);
        }

        bool pass_flag =
            ok(s_first) && ok(s_last) && ok(e_first) && ok(e_last) &&
            (s_first <  thresh_low)  &&
            (s_last  >  thresh_high) &&
            (e_first >  thresh_high) &&
            (e_last  <  thresh_low);

        string pass_fail = pass_flag ? "PASS" : "FAIL";

        // ---------- Output TRUE_BOUNDARY with PASS/FAIL ----------
        if (true_l_start>0 && true_l_end>0 && true_r_start>0 && true_r_end>0) {
            cout << "# TRUE_BOUNDARY " << c.chr << " "
                 << true_l_start << " " << true_l_end << " "
                 << true_r_start << " " << true_r_end << " "
                 << TSD << "\t" << motif << "\t" << pass_fail << "\n";
        } else {
            cout << "# TRUE_BOUNDARY " << c.chr << " "
                 << (true_l_start>0? to_string(true_l_start):"NA") << " "
                 << (true_l_end>0  ? to_string(true_l_end):"NA")   << " "
                 << (true_r_start>0? to_string(true_r_start):"NA") << " "
                 << (true_r_end>0  ? to_string(true_r_end):"NA")   << " "
                 << TSD << "\t" << motif << "\t" << pass_fail << "\n";
        }
    }

    return 0;
}
