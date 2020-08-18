DWGFlute : UGen
{
	*ar { arg freq=400, pm = 1, endr = 0.5, jetr = 0.25 , jetRa = 0.33, gate = 1, release = 0.1;
		^this.multiNew('audio',freq, pm,endr,jetr,jetRa, gate,release);
	}
}

DWGFlute2 : UGen
{
	*ar { arg freq=400, pm = 1, endr = 0.5, jetr = 0.25 , jetRa = 0.33, gate = 1, release = 0.1;
		^this.multiNew('audio',freq, pm,endr,jetr,jetRa, gate,release);
	}
}
