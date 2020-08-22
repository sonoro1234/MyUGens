PluckSynth : UGen
{
	*ar { arg freq = 440, amp = 1, gate = 1, pos = 1/7, c1 = 0.25, c3 = 5, release = 0.1, f = 0, m = 0, k = 0, r = 0, l = 0.65, ra = 0.001, rho = 7850;
		^this.multiNew('audio',freq, amp,gate,pos,c1,c3,release,f,m,k,r,l,ra,rho);
	}
}
