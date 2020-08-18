MembraneCircleV : UGen {
	*ar { arg excitation, tension=0.05, loss = 0.99999,ewidth=0.5,epos=0,size=1,a1=0.5,doprint=0, mul = 1.0, add = 0.0;
		^this.multiNew('audio', excitation, tension, loss,ewidth,epos,size,a1,doprint).madd(mul, add)
	}
	checkInputs { ^this.checkSameRateAsFirstInput }
}

MembraneHexagonV : UGen {
	*ar { arg excitation, tension=0.05, loss = 0.99999,ewidth=0.5,epos=0,size=1,a1=0.5,doprint=0, mul = 1.0, add = 0.0;
		^this.multiNew('audio', excitation, tension, loss,ewidth,epos,size,a1,doprint).madd(mul, add)
	}
}
