DWGClarinet3 : UGen
{
	*ar { arg freq=440, pm=1, pc=1, m=0.8, gate=1, release=0.01, c1=0.25, c3=7;
		^this.multiNew('audio',freq, pm,pc,m, gate,release,c1,c3);
	}
}
