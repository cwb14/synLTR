cat /home/chris/data/synLTR2/systematic_test/piecewise_regression/v3/lineage1/get_divK_Kmer2LTR_benchmarking.pl
#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname(__FILE__); # assume Distributions.pm is in the same directory
use Distributions;
use Getopt::Long;

# ---------------------------
# Defaults and usage
# ---------------------------
my $window  = 0.001;  # window size for scanning K
my $bootnum = 100;    # number of bootstraps (renamed from $b to avoid clash with sort's $b)
my $mode    = 1;      # bootstrapping mode (1 = with replacement; 0 = without replacement, 90%)
my $help    = 0;      # help flag

# NEW: scaling multipliers (old behavior; still supported)
my $smul    = 1.0;    # shared multiplier (-sm)
my $umul    = 1.0;    # unique multiplier (-um)

# NEW: target counts (take precedence over multipliers if provided)
my $scount;           # target count for shared/full (value==2)  (-sc)
my $ucount;           # target count for unique/flanking (value==1) (-uc)

my $usage = "\nFitting the simple piecewise regression (aka. fishing regression) to find the optimum K
Usage: perl get_divK_Kmer2LTR2.pl [options] <input.tsv>

Input format (TSV/XLS-style):
  - Column 1: ID string. If it begins with 'shared_', treat as Site_type:Full; otherwise Site_type:Flanking
  - Column 11 (1-based): divergence (if the value looks like identity (>0.6), it will be converted to divergence as 1 - value)

Options:
  -b [int]     Number of bootstraps (default: 100)
  -m [1,0]     Bootstrapping approach:
               1 = sample all datapoints with replacement
               0 = sample 90% of datapoints without replacement (default: 1)
  -w [float]   Window size for scanning K (default: 0.001)
  -sm [float]  Shared multiplication factor (value==2; default: 1.0)
  -um [float]  Unique multiplication factor (value==1; default: 1.0)
  -sc [int]    TARGET COUNT for shared/full (value==2). If set, overrides -sm.
  -uc [int]    TARGET COUNT for unique/flanking (value==1). If set, overrides -um.
  -h           Display this help message

Count vs Multiplier:
  - If -sc/-uc are provided, they take precedence over -sm/-um respectively.
  - Examples (with 1000 shared and 2500 unique originally):
      -sc 1500 -uc 2400  => shared: include all 1000 once, then duplicate 500 randomly to reach 1500;
                             unique: sample 2400 without replacement from the 2500.

Output:
  First line (added): SUMMARY\\torig_S=..\\torig_U=..\\torig_S/U=..\\trescaled_S=..\\trescaled_U=..\\trescaled_S/U=..
  Then per bootstrap:  OptimumK\\t[lowerCI,upperCI]

Notes:
  - If there are zero 'Flanking' points (value==1) in a bootstrap, prints 0.00000 [0.00000,0.00000]
  - If fewer than 50 'Full' points (value==2) lie above K, prints NA [NA,NA]\n\n";

GetOptions(
    'b=i'  => \$bootnum,
    'm=i'  => \$mode,
    'w=f'  => \$window,
    'sm=f' => \$smul,
    'um=f' => \$umul,
    'sc=i' => \$scount,
    'uc=i' => \$ucount,
    'h'    => \$help,
) or die $usage;

die $usage if $help;
die $usage if @ARGV != 1;

# Validate bootstrapping mode
if ($mode != 0 && $mode != 1) {
    die "Error: Invalid value for -m. Must be 1 (with replacement) or 0 (without replacement).\n$usage";
}

# Ensure that the number of bootstraps is at least 1
if ($bootnum < 1) {
    die "Error: Number of bootstraps (-b) must be at least 1.\n$usage";
}

# Validate multipliers
if ($smul < 0 || $umul < 0) {
    die "Error: -sm and -um must be >= 0.\n$usage";
}

# Validate counts if provided
if (defined $scount && $scount < 0) {
    die "Error: -sc must be >= 0.\n$usage";
}
if (defined $ucount && $ucount < 0) {
    die "Error: -uc must be >= 0.\n$usage";
}

# Parameters for scanning K and CI
my $mindist    = 0.0;   # minimum distance
my $maxdist    = 0.1;   # maximum searching distance
my $percentile = 0.05;  # two-tailed CI percentile (5%)

# ---------------------------
# Helpers
# ---------------------------

sub fisher_yates_shuffle_inplace {
    my ($aref) = @_;
    for (my $i = @$aref - 1; $i > 0; $i--) {
        my $j = int(rand($i + 1));
        @$aref[$i, $j] = @$aref[$j, $i];
    }
}

# Return a NEW arrayref containing items sampled according to factor
# - factor < 1: downsample without replacement to round(factor * N)
# - factor >=1: include all once, then:
#       add floor(factor-1) full extra copies of the entire set,
#       then add round(frac * N) extra single copies of randomly chosen items (without replacement)
sub sample_by_factor {
    my ($items_aref, $factor) = @_;
    my @items = @$items_aref;
    my $N = scalar @items;

    return [] if $N == 0 || $factor == 0;

    my @out;

    if ($factor < 1.0) {
        my $keep = int($factor * $N + 0.5);  # round to nearest
        $keep = 0 if $keep < 0;
        $keep = $N if $keep > $N;
        my @idx = (0 .. $N-1);
        fisher_yates_shuffle_inplace(\@idx);
        @idx = @idx[0 .. $keep-1] if $keep < @idx;
        foreach my $k (@idx) {
            push @out, $items[$k];
        }
        return \@out;
    }

    # factor >= 1
    push @out, @items;  # one base copy of everything

    my $extra_full = int($factor - 1.0);   # number of full extra passes
    my $frac       = ($factor - 1.0) - $extra_full;

    # add full extra copies
    for (my $t = 0; $t < $extra_full; $t++) {
        push @out, @items;
    }

    # add a fractional second pass: choose round(frac * N) unique items, each gets one extra copy
    if ($frac > 0) {
        my $extra_count = int($frac * $N + 0.5);
        if ($extra_count > 0) {
            my @idx = (0 .. $N-1);
            fisher_yates_shuffle_inplace(\@idx);
            @idx = @idx[0 .. $extra_count-1] if $extra_count < @idx;
            foreach my $k (@idx) {
                push @out, $items[$k];
            }
        }
    }

    return \@out;
}

# NEW: Return a NEW arrayref containing exactly 'target' items
# - if target <= N: sample without replacement to exactly target
# - if target > N : include floor(target/N) full copies; then add (target % N) unique extra items (one extra copy each)
sub sample_to_count {
    my ($items_aref, $target) = @_;
    my @items = @$items_aref;
    my $N = scalar @items;

    return [] if $N == 0 || !defined $target || $target <= 0;

    my @out;

    if ($target <= $N) {
        my @idx = (0 .. $N-1);
        fisher_yates_shuffle_inplace(\@idx);
        @idx = @idx[0 .. $target-1] if $target < @idx;
        foreach my $k (@idx) { push @out, $items[$k]; }
        return \@out;
    }

    # target > N
    my $full = int($target / $N);
    my $rem  = $target % $N;

    # add full copies
    for (my $t = 0; $t < $full; $t++) {
        push @out, @items;
    }

    # add remainder as unique extras (one extra copy each)
    if ($rem > 0) {
        my @idx = (0 .. $N-1);
        fisher_yates_shuffle_inplace(\@idx);
        @idx = @idx[0 .. $rem-1] if $rem < @idx;
        foreach my $k (@idx) { push @out, $items[$k]; }
    }

    return \@out;
}

# Safe ratio printer (S/U), returns string with 4 decimals or NA if denom==0
sub ratio_str {
    my ($s, $u) = @_;
    return "NA" if !defined $u || $u == 0;
    return sprintf("%.4f", $s / $u);
}

# ---------------------------
# Read input and collect data
# ---------------------------
my $file = $ARGV[0];
open my $fh, '<', $file or die "Cannot open file '$file': $!\n$usage";

my @all_flanking;  # value==1 (Unique/Flanking)
my @all_full;      # value==2 (Shared/Full)

while (<$fh>) {
    chomp;
    next if /^$/;           # skip empty
    next if /^\s*#/;        # skip comments

    my @fields = split(/\t/, $_);

    # Need at least 11 columns (1-based col 11 => 0-based index 10)
    if (@fields < 11) {
        warn "Skipping malformed line (needs >= 11 columns): $_\n";
        next;
    }

    my $id  = $fields[0];   # column 1
    my $raw = $fields[10];  # 1-based Column 11 => divergence or identity-to-convert

    # Site type from column 1
    my $value = (defined $id && $id =~ /^shared_/) ? 2 : 1;

    # Parse divergence (allow numeric; if identity-like (>=0.6), convert to divergence)
    next if !defined $raw || $raw eq '' || $raw =~ /NA/i;

    if ($raw !~ /^-?\d*\.?\d+(?:[eE][+-]?\d+)?$/) {
        warn "Skipping non-numeric divergence in column 11: $raw\tLine: $_\n";
        next;
    }

    my $div = $raw + 0;

    # If value looks like identity (>= 0.6), convert to divergence
    if ($div >= 0.6) {
        $div = 1.0 - $div;
    }

    # Sanity: clamp to [0,1]
    if ($div < 0) { $div = 0; }
    if ($div > 1) { $div = 1; }

    if ($value == 1) {
        push @all_flanking, [1, $div];
    } else {
        push @all_full,     [2, $div];
    }
}
close $fh;

my $orig_S = scalar @all_full;
my $orig_U = scalar @all_flanking;
my $orig_total = $orig_S + $orig_U;
die "Error: No valid data points found in the input file.\n" if $orig_total == 0;

# ---------------------------
# Apply scaling BEFORE bootstrapping
# ---------------------------
my $scaled_flanking_aref; # Unique/Flanking (value==1)
my $scaled_full_aref;     # Shared/Full (value==2)

if (defined $ucount) {
    $scaled_flanking_aref = sample_to_count(\@all_flanking, $ucount);
} else {
    $scaled_flanking_aref = sample_by_factor(\@all_flanking, $umul);
}

if (defined $scount) {
    $scaled_full_aref = sample_to_count(\@all_full, $scount);
} else {
    $scaled_full_aref = sample_by_factor(\@all_full, $smul);
}

my $rescaled_S = scalar @$scaled_full_aref;
my $rescaled_U = scalar @$scaled_flanking_aref;

my @all_input = (@$scaled_flanking_aref, @$scaled_full_aref);

# Precompute totals over the ADJUSTED data
my $total_n          = scalar @all_input;
my $total_count_ones = $rescaled_U;
my $total_count_twos = $rescaled_S;

# total sum of divergence (adjusted)
my $total_sum_div = 0;
$total_sum_div += $_->[1] for @all_input;

# ---------------------------
# Print the requested SUMMARY line
# ---------------------------
my $orig_ratio     = ratio_str($orig_S,     $orig_U);
my $rescaled_ratio = ratio_str($rescaled_S, $rescaled_U);

print "SUMMARY\t";
print "orig_S=$orig_S\torig_U=$orig_U\torig_S/U=$orig_ratio\t";
print "rescaled_S=$rescaled_S\trescaled_U=$rescaled_U\trescaled_S/U=$rescaled_ratio\n";

# If there are zero 'Flanking' (value==1) points overall after rescaling, keep original behavior
if ($total_count_ones == 0) {
    for (my $i = 0; $i < $bootnum; $i++) {
        print "0.00000\t[0.00000,0.00000]\n";
    }
    exit;
}

# ---------------------------
# Bootstrapping
# ---------------------------
for (my $bootstrap = 0; $bootstrap < $bootnum; $bootstrap++) {
    my @input;
    my $sum_div    = 0;
    my $n          = 0;
    my $count_ones = 0;
    my $count_twos = 0;

    if ($bootstrap == $bootnum - 1) {
        # Last bootstrap: sample 100% without replacement from ADJUSTED data
        @input         = @all_input;
        $n             = $total_n;
        $sum_div       = $total_sum_div;
        $count_ones    = $total_count_ones;
        $count_twos    = $total_count_twos;
    } else {
        if ($mode == 1) {
            # Sample all datapoints with replacement
            for (my $i = 0; $i < $total_n; $i++) {
                my $rand_index = int(rand($total_n));
                my $pair       = $all_input[$rand_index];
                push @input, [$pair->[0], $pair->[1]];
                $sum_div += $pair->[1];
                $n++;
                $count_ones++ if $pair->[0] == 1;
                $count_twos++ if $pair->[0] == 2;
            }
        } else {
            # Sample 90% without replacement
            my $sample_size = int(0.9 * $total_n);
            $sample_size = 1 if $sample_size < 1;

            my @shuffled = @all_input;
            fisher_yates_shuffle_inplace(\@shuffled);

            for (my $i = 0; $i < $sample_size; $i++) {
                my $pair = $shuffled[$i];
                push @input, [$pair->[0], $pair->[1]];
                $sum_div += $pair->[1];
                $n++;
                $count_ones++ if $pair->[0] == 1;
                $count_twos++ if $pair->[0] == 2;
            }
        }
    }

    # If no 'Flanking' points in this bootstrap, report zeros
    if ($count_ones == 0) {
        print "0.00000\t[0.00000,0.00000]\n";
        next;
    }

    # ---------------------------
    # SSE scan over K
    # ---------------------------
    my %sse;
    my ($x_min, $x_max) = ("NA", "NA");
    my $SSE = "NA";

    for (my $x = $mindist; $x <= $maxdist + 1e-12; $x += $window) {
        my $current_x = sprintf("%.5f", $x);
        my $sse = 0;

        foreach my $pair (@input) {
            my ($val, $div) = @{$pair}[0, 1];
            if ($div < $current_x) {
                $sse += abs($val - 1);
            } else {
                $sse += abs($val - 2);
            }
        }

        $sse{$current_x} = $sse;

        if ($SSE eq "NA") {
            ($x_min, $x_max, $SSE) = ($current_x, $current_x, $sse);
        }

        if ($sse <= $SSE) {
            $x_min = $current_x if $sse < $SSE;
            $x_max = $current_x;
            $SSE   = $sse;
        }
    }

    my $distance = sprintf("%.5f", ($x_min + $x_max) / 2);
    my $s2 = ($n > 2) ? $SSE / ($n - 2) : 0;

    # Count 'Full' (value==2) points with divergence > K
    my $count_twos_above_k = 0;
    foreach my $pair (@input) {
        my ($val, $div) = @{$pair}[0, 1];
        if ($val == 2 && $div > $distance) {
            $count_twos_above_k++;
        }
    }

    # If fewer than 50 'Full' above K, report NA
    if ($count_twos_above_k < 50) {
        print "NA\t[NA,NA]\n";
        next;
    }

    # ---------------------------
    # Confidence interval via F
    # ---------------------------
    my $f = Statistics::Distributions::fdistr(1, $n - 1, $percentile);

    my ($lower, $upper) = ("NA", "NA");
    if ($s2 > 0) {
        foreach my $x (sort { $a <=> $b } keys %sse) {
            my $F_x = ($sse{$x} - $SSE) / $s2; # Toms & Lesperance (2003)
            if ($F_x <= $f) {
                $lower = $x if $lower eq "NA";
                $upper = $x;
            }
        }
    }

    # If no CI region found, fall back to 2nd-highest divergence as a rough estimate
    if ($lower eq "NA") {
        my @div_values = map { $_->[1] } @input;
        @div_values = sort { $a <=> $b } @div_values;
        if (@div_values >= 2) {
            $distance = sprintf("%.5f", $div_values[-2]);
        } elsif (@div_values == 1) {
            $distance = sprintf("%.5f", $div_values[0]);
        } else {
            $distance = "NA";
        }
        $lower = "NA";
        $upper = "NA";
    }

    print "$distance\t[$lower,$upper]\n";
}
# END
