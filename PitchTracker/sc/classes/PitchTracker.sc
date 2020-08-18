PitchTracker : MultiOutUGen {

	*kr { arg in = 0.0, size = 1024, hop = 0.25, dfreq = 440, thresh = 0.85, kind = 0, useef = 0, bufnum = -1, t_clear = 0;
		^this.multiNew('control', in,size,hop,dfreq,thresh,kind,useef,bufnum,t_clear)
	}
	init { arg ... theInputs;
		inputs = theInputs;
		^this.initOutputs(3, rate);
	}
}