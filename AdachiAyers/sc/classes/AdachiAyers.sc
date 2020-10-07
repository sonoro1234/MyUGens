AdachiAyers : UGen
{
	*ar { arg flip = 231 ,pO = 5000,radio = 0.0085,buffnum = -1,buffnum2 = -1,buffnum3 = -1,yequil = 0.0014,gate = 1,delay = 100;
		var delaycl = delay.min(1400).max(1);
		^this.multiNew('audio',flip,pO,radio,buffnum,buffnum2,buffnum3,yequil,gate,delay);
	}
}