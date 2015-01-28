Code to reproduce results shown in arXiv:1404.7353
--------------------------------------------------

This code is intended to allow easy reproduction of the results given
in (arXiv:1404.7353)[http://arxiv.org/abs/arXiv:1407.0408].

Usage
-----

You will need to carry out the following steps:

1. Make sure you have FastJet (fastjet-config should be in your path),
   at least version 3.1.0, fastjet contrib version > 1.014 (should be
   installed to the same location as fastjet) and the Boost C++ library.

2. Get the lightweight event-analysis framework from the (Pileup 2014
   workshop)[https://indico.cern.ch/event/306155/]: in some directory,
   do the following

        git clone https://github.com/PileupWorkshop/2014PileupWorkshop.git
        cd 2014PileupWorkshop/Framework
        ./mkmk
        make 

   (If you already got a copy from the pileup workshop in May 2014,
   update it, because we use some new features added since then)

3. Then edit the file mkmk in this directory (the same as this
   README.md file) so that PU14CODE points to the location of the
   pileup workshop top-level directory. Then build things with

        ./mkmk
        make

   On a Mac, if you got Boost from some source other than macports,
   you may need to edit the parts that point to /opt/local.

4. Get some event files. These are at
   http://puws2014.web.cern.ch/puws2014/events/ and also on afs at
   /afs/cern.ch/user/p/puws2014/public/events/

   You will need at least:

   - http://puws2014.web.cern.ch/puws2014/events/lhc14-pythia8-4C-Zprime500-nobdecay-noUE-nev1e5.pu14.gz
   - http://puws2014.web.cern.ch/puws2014/events/lhc14-pythia8-4C-minbias-nev1e6.pu14.gz

Then you're ready to run. You can get most of the CHS results with a
command such as the following

        ./baseline  -chs -massless \
          -hard events/lhc14-pythia8-4C-Zprime500-nobdecay-noUE-nev1e5.pu14.gz \
          -pileup events/lhc14-pythia8-4C-minbias-nev1e6.pu14.gz \
          -npu 100 -nev 100 -out results.txt

where events/ is the directory that contains the event files (or a
symbolic link to it).

100 events (-nev 100) won't give you very reliable answers but should
run in 15-30s. The -npu 100 flag adds exactly 100 pileup events to
each hard event.

To start getting numerically stable results you should increase the
number of events to at least 10000. The results should then agree with
those in

  file:reference-results-npu100.txt

(if you need more than 10k events, download a larger minbias file).


----------------------------------------------------------------------
Reproducing results in arXiv:1407.0408v2
========================================

Create a directory results/ and then try out the following commands in
order to get a subset of the results shown in the paper.

The code used to generate the actual plots in the paper was mostly
different from the code in this directory here. We believe that the
results are consistent across codes.

Also some minor details of the running were different between what's given
here and what is given below. In particular, fig.4 (correlation
coefficients v. NPU) in the paper is given for bins of NPU ranging
e.g. from 1-10, 11-20, etc., whereas the program runs just for a
single fixed pileup.

Finally, to run the examples below, you'll need a wider range of event
files (all should be in the same directory linked to above).

## Figure 4 (correlation coefficient for Z' decay versus NPU)

Run commands such as the following (with and without B's in the Z' decay)

    ./baseline -nev 100000 -chs -massless \
               -pileup events/lhc14-pythia8-4C-minbias-nev1e7.pu14.gz -npu 105 \
               -hard events/lhc14-pythia8-4C-Zprime500-nobdecay-noUE-nev1e5.pu14.gz \
               -out results/Zp500-npu105.res

    ./baseline -nev 100000 -chs -massless \
               -pileup events/lhc14-pythia8-4C-minbias-nev1e7.pu14.gz -npu 105 \
               -hard events/lhc14-pythia8-4C-Zprime2uds500-nobdecay-noUE-nev1e5.pu14.gz \
               -out results/Zp500uds-npu105.res

You'll need to scan over the "-npu 105" argument to get the full range
of number of pileup events.

The reference results for the first command are in file:ref-results/Zp500-npu105.res

## Figure 5 ("fingerprint" plots for Z' events)

The same, but replacing "-npu 100" with "-upu 140" (and changing the
output file name). Just a single run is needed since the -upu argument
gives pileup distributed uniformly over the range 1..140.

Note that the program only outputs dispersions, not the 90% width
information. For that you'd need to extend the program to write out
histograms or ntuples and deduce the 90% width from those.

## Figure 6 (50GeV dijets with R=0.4 jets)

Run the following:

    ./baseline -nev 100000 -chs -massless \
               -pileup events/lhc14-pythia8-4C-minbias-nev1e7.pu14.gz -upu 140 \
               -hard events/lhc14-pythia8-4C-dijetsel50-noUE-nevsel1e5.pu14.gz \
               -R 0.4 -Rsub 0.2 -jet-ptmin 50 \
               -out results/dijets50-R04-Rsub02-upu140.res

## Figure 7 (WW with f=0.05 trimmed R=1 jets)

Run

    ./baseline -nev 100000 -chs -massless \
               -pileup events/lhc14-pythia8-4C-minbias-nev1e7.pu14.gz -upu 140 \
               -hard events/lhc14-pythia8-4C-WW500-noUE-nev1e5.pu14.gz \
               -jet-ptmin 500 -fcut 0.05 \
               -out results/WW500-fcut05-upu140.res

Note that this run (which recycles a hard event file taken from the 2014
Pileup Workshop) differs slightly from that in the paper, where the
generation cut for the events is 450 GeV (here it is 500 GeV). Both
apply a 500 GeV jet selection cut.

