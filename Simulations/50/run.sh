# define sequence length
if [ -z "$SL" ]; then
export SL=100
fi

# create 100 beast XML files
applauncher CoverageTest -xml input.out.xml -w . -log truth.log -tree truth.trees -o xml
cd xml
mkdir run-sl$SL
cd run-sl$SL


# run 100 beast XML files
for i in {0..4};   do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {5..9};   do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {10..14}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {15..19}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {20..24}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {25..29}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {30..34}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {35..39}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {40..44}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {45..49}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {50..54}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {55..59}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {60..64}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {65..69}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {70..74}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {75..79}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {80..84}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {85..89}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {90..94}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &
for i in {95..99}; do beast -D sl=$SL -seed 127$i -overwrite ../analysis-out$i.xml > out$i 2>&1 ; done &

# create summary of trace logs	
loganalyser -threads 20 -oneline a*.log > summary

# calculate coverage of true values by estimates
applauncher beastvalidation.experimenter.CoverageCalculator -log ../../../truth.log -logA summary

applauncher beastvalidation.experimenter.CladeCoverageCalculator -tr ../../../truth.trees -pre analysis-out


