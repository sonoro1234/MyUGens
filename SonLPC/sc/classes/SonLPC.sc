SonLPC : UGen
{
	*ar { arg buff = -1, in=0, hop=0.5, poles=10;
		^this.multiNew('audio',buff,in,hop,poles);
	}
}

SonLPCSynth : UGen
{
	*ar { arg chain = -1;
		^this.multiNew('audio',chain);
	}
}

SonLPCSynthIn : UGen
{
	*ar { arg chain = -1, in=0;
		^this.multiNew('audio',chain, in);
	}
}


