# MyUGens
Some Supercollider plugins:

AdachiAyers for trombone trumpet physical model

DWGClarinet

DWGFlute

DWGReverb

Karplus 

PitchTracker: when Tartini is not enough.

IIRf for infinite impulse response filters

KLJunction for Kelly-Lochbaum filter implementation

SonLPC for linear prediction analysis and synthesis.

PluckSynth emulation of a plucked string from Fredrik Eckerholm, Gianpaolo Evangelista paper

##compilation

On a sibling folder to MyUGens folder:
<code>
cmake -G"Your generator or skip G to let cmake choose" -DCMAKE_BUILD_TYPE=Release  -DSUPERNOVA=ON -DCMAKE_INSTALL_PREFIX="./scinstall or were you want to install" -DSC_PATH="../supercollider/ or where to find SC repo folder"  ../MyUGens
</code>

wait command to finish and
<code>
cmake --build . --target install
</code>
