require"victor.utils"
require"victor.ta"
require"victor.filter_design"
require"tubes.number2string"
require"tubes.matrices"
require"victor.gslpoly"
dir = [[C:\SupercolliderRepos\Mios\DWGClarinet\]]
Radiation = LevineScavone
DatosCopa={31e-3-29.5e-3	,.5*24.1e-3	,.5*23.4e-3,
					29.5e-3	-  27.5e-3	,.5*23.4e-3	,.5*22.2e-3,
					27.5e-3 -  25.5e-3	,.5*22.2e-3	,.5*20.6e-3,
					25.5e-3 -  23.5e-3	,.5*20.6e-3	,.5*18.9e-3,
					23.5e-3 -  21.5e-3	,.5*18.9e-3	,.5*16.8e-3,
					21.5e-3 - 19.5e-3	,.5*16.8e-3	,.5*14.3e-3,
					19.5e-3 - 17.5e-3	,.5*14.3e-3	,.5*11.4e-3,
					17.5e-3 - 15.5e-3	,.5*11.4e-3	,.5*8.4e-3,
					15.5e-3 - 13.5e-3	,.5*8.4e-3	,.5*7.e-3,
					13.5e-3 - 11.5e-3	,.5*7.e-3	,.5*6.8e-3,
					11.5e-3 - 11.e-3	,.5*6.8e-3	,.5*6.7e-3,
					11.e-3 - 0.e-3		,.5*6.7e-3 	,.5*6.7e-3,};
DatosBell={     5.7e-2		,1.225e-2	,1.225e-2,
					5.8e-2	-  5.7e-2	,1.225e-2	,1.265e-2,
					 5.9e-2	-  5.8e-2	,1.265e-2	,1.265e-2,
					 8.3e-2 -  5.9e-2	,1.265e-2	,1.3e-2,
					12.4e-2 -  8.3e-2	,1.3e-2		,1.4e-2,
					16.2e-2 - 12.4e-2	,1.4e-2		,1.5e-2,
					20.3e-2 - 16.2e-2	,1.5e-2		,1.65e-2,
					25.3e-2 - 20.3e-2	,1.65e-2	,1.85e-2,
					31.2e-2 - 25.3e-2	,1.85e-2	,2.15e-2,
					36e-2   - 31.2e-2	,2.15e-2	,2.45e-2,
					40e-2   - 36e-2		,2.45e-2	,2.8e-2,
					43.2e-2 - 40e-2		,2.8e-2 	,3.2e-2,
					45.7e-2 - 43.2e-2	,3.2e-2 	,3.7e-2,
					48.2e-2 - 45.7e-2	,3.7e-2 	,4.45e-2,
					50e-2   - 48.2e-2	,4.45e-2	,5.2e-2,
					51.5e-2 - 50e-2		,5.2e-2 	,6.2e-2,
					52.7e-2 - 51.5e-2	,6.2e-2 	,7.45e-2,
					53.6e-2 - 52.7e-2	,7.45e-2	,8.95e-2,
					54.4e-2 - 53.6e-2	,8.95e-2	,11.45e-2,};
--[[					
chokelength=4.e-2;
Datos[trozoscopa].l=chokelength;	 
Datos[trozoscopa].r1=0.3*radiotubo --Minimradiocopa;
Datos[trozoscopa].r2=radiotubo;
	
Datos[trozoscopa+1].l=170e-2;	 
Datos[trozoscopa+1].r1=radiotubo;
Datos[trozoscopa+1].r2=radiotubo;
--]]
radiotubo = DatosBell[2]
Tubo ={170e-2,radiotubo,radiotubo}
Choke ={4.e-2,DatosCopa[#DatosCopa],radiotubo}
function CalcVolumen(Seg)
	local sumV=0.;
	local RR,Rr,rr;
	for i,seg in ipairs(Seg) do
		RR=seg.r1*seg.r1;
		Rr=seg.r1*seg.r2;
		rr=seg.r2*seg.r2;
		sumV = sumV+math.pi*seg.l*(RR+rr+Rr)/3.;
	end
	return sumV;
end

function data2seg(Datos)
	local dCopa = TA{}
	for i=1,#Datos,3 do
		local seg = {}
		seg.l = Datos[i]
		seg.r1 = Datos[i+1]
		seg.r2 = Datos[i+2]
		dCopa[#dCopa + 1] = seg
	end
	return dCopa
end

function TroceaSeg(len,r1,r2,trozos)
	local stepr = (r2 - r1)/trozos
	local ltrozo = len/trozos
	local seg = {}
	local defecto = 0
	for i=1,trozos do
		seg[#seg + 1] = {l=ltrozo,r1=r1 + stepr*(i-1),r2=r1 + stepr*(i)}
		defecto = defecto + math.abs(seg[#seg].r2 - seg[#seg].r1)*seg[#seg].l
	end
	print("defecto:",defecto)
	return seg
end
function Segmenta(datos,faltamin)
	local datosseg = TA()
	for i,v in ipairs(datos) do
		local falta = (v.r2 - v.r1)*v.l
		--local tro = math.log(math.abs(falta/faltamin))/math.log(2)
		local tro = math.sqrt(math.abs(falta/faltamin))
		local trozos = math.max(1,math.floor(tro))
		print(falta,falta/faltamin,tro,trozos)
		datosseg = datosseg..TroceaSeg(v.l,v.r1,v.r2,trozos)
	end
	return datosseg
end


function MatrixSegmentsBAK(Segs,freqs)
	
	local Matrices = {}
	for i2,freq in ipairs(freqs) do
		Matrices[i2] =  matrix.unit(2)
	end
	local tc = os.clock()
	for i,seg in ipairs(Segs) do
		
		print("viisco "..i)
		tc = os.clock()
		local viscos,Z0 = KViscosity(freqs,0.5*(seg.r1+seg.r2))
		print(os.clock() -tc)
		local invZ0 = 1/Z0
		print("freqw "..i)
		tc = os.clock()
		--for i2,_ in ipairs(freqs) do
		for i2,MAT in ipairs(Matrices) do
			--local K = viscos[i2]
			--local L = seg.l
			local KL  = viscos[i2] * seg.l
			local coskl = complex.cos(KL)
			local sinkl = -complex.i*complex.sin(KL)
			--Matrices[i2] =  Matrices[i2] * matrix.cdef{{coskl,Z0*sinkl},{invZ0*sinkl,coskl}}
			Matrices[i2] =  MAT * matrix.cdef{{coskl,Z0*sinkl},{invZ0*sinkl,coskl}}
		end
		print(os.clock() -tc)
	end
	return Matrices
end
---version conos
function MatrixSegmentsF(Segs,freqs,viscofunc)
	local viscofunc =  viscofunc or KViscosity
	local Matrices = {}
	for i2,freq in ipairs(freqs) do
		Matrices[i2] =  matrix.unit(2)
	end
	local prealloc = matrix.cdef{{1,1},{1,1}}
	local tc = os.clock()
	for i,seg in ipairs(Segs) do
		
		print("viisco "..i)
		tc = os.clock()
		local viscos,Z0 = viscofunc(freqs,0.5*(seg.r1+seg.r2))
		print(os.clock() -tc)
		--local invZ0 = 1/Z0
		print("freqw "..i)
		tc = os.clock()
		--for i2,_ in ipairs(freqs) do
		local MatTemp = {}
		for i2,MAT in ipairs(Matrices) do
			
			--MatTemp[i2] = MAT * TuboMatrix(seg.l,0.5*(seg.r1+seg.r2),freqs[i2],viscos[i2],Z0)
			MatTemp[i2] = MAT * ConoMatrix2(seg.l,seg.r1,seg.r2,viscos[i2],Z0,freqs[i2],prealloc)
		end
		Matrices = MatTemp
		print(os.clock() -tc)
	end
	return Matrices
end
function MatrixSegmentsScavone(Segs,freqs)
	local viscofunc =  viscofunc or ViscoScavone
	local Matrices = {}
	for i2,freq in ipairs(freqs) do
		Matrices[i2] =  matrix.unit(2)
	end
	local prealloc = matrix.cdef{{1,1},{1,1}}
	local tc = os.clock()
	for i,seg in ipairs(Segs) do
		
		local tc = os.clock()
		--for i2,_ in ipairs(freqs) do
		local MatTemp = {}
		for i2,MAT in ipairs(Matrices) do
			local K,Z0 = ViscoScavone(freqs[i2],seg.r1,seg.r2)
			--MatTemp[i2] = MAT * TuboMatrix(seg.l,0.5*(seg.r1+seg.r2),freqs[i2],viscos[i2],Z0)
			MatTemp[i2] = MAT * ConoMatrixScavoneVisco(seg.l,seg.r1,seg.r2,K,Z0,freqs[i2])
		end
		Matrices = MatTemp
		print(os.clock() -tc)
	end
	return Matrices
end
--version optimizada
function MatrixSegments(Segs,freqs,viscofunc)
	local viscofunc =  viscofunc or KViscosity
	local Matrices = {}
	for i2,freq in ipairs(freqs) do
		Matrices[i2] =  {1,0,0,1}
	end
	local tc = os.clock()
	for i,seg in ipairs(Segs) do
		
		print("viisco "..i)
		tc = os.clock()
		local viscos,Z0 = viscofunc(freqs,0.5*(seg.r1+seg.r2))
		print(os.clock() -tc)
		local invZ0 = 1/Z0
		print("freqw "..i)
		tc = os.clock()
		--for i2,_ in ipairs(freqs) do
		local MatTemp = {}
		for i2,MAT in ipairs(Matrices) do
			--local K = viscos[i2]
			--local L = seg.l
			local KL  = viscos[i2] * seg.l
			local coskl = complex.cos(KL)
			local sinkl = -complex.i*complex.sin(KL)
			local Z0sinkl = Z0*sinkl
			local invZ0sinkl = invZ0*sinkl
			--Matrices[i2] =  Matrices[i2] * matrix.cdef{{coskl,Z0*sinkl},{invZ0*sinkl,coskl}}
			--Matrices[i2] =  MAT * matrix.cdef{{coskl,Z0*sinkl},{invZ0*sinkl,coskl}}

			local A = MAT[1]*coskl + MAT[2]*invZ0sinkl
			local B = MAT[1]*Z0sinkl + MAT[2]*coskl
			local C = MAT[3]*coskl + MAT[4]*invZ0sinkl
			local D = MAT[3]*Z0sinkl + MAT[4]*coskl
			MatTemp[i2] = {A,B,C,D}
		end
		Matrices = MatTemp
		print(os.clock() -tc)
	end
	for i,m in ipairs(Matrices) do
		Matrices[i] = matrix.cdef{{m[1],m[2]},{m[3],m[4]}}
	end
	return Matrices
end
function GetFreqs(Nf,fmin,fmax)
	local fmin = fmin or 0
	local fmax = fmax or FRECMAX

	local fstep = (fmax - fmin)/Nf
	local freqs = {}
	for i=0,Nf do
		freqs[#freqs + 1] = fstep * i + fmin
	end
	return freqs
end
function Mat2POPI(mats,Zls)
	local popis = {}
	for i,Zl in ipairs(Zls) do
		local res = mats[i] * matrix.vec{Zl,1}
		if i == 1 then
			print("POPI res[1],res[2]",res[1],res[2],Zl)
			popis[i] = Zl/(res[1]) 
		else
			popis[i] = Zl/(res[1])
		end
	end
	return popis
end
function Mat2POPMas(mats,Zls,radio)
	local Zc = RO_AIRE*C_AIRE/(math.pi*radio^2)
	local matConv = matrix.cdef{{1/2,Zc/2},{1/2,-Zc/2}}
	local popis = {}
	for i,Zl in ipairs(Zls) do
		local res = matConv * mats[i] * matrix.vec{Zl,1}
		if i == 1 then
			print("POPI res[1],res[2]",res[1],res[2],Zl)
			popis[i] = Zl/(res[1]) 
		else
			popis[i] = Zl/(res[1])
		end
	end
	return popis
end
function Mat2Zin(mats,Zls)
	local popis = {}
	local res
	for i,Zl in ipairs(Zls) do
		if Zl ~= math.huge then
		--if true then
			res = mats[i] * matrix.vec{Zl,1}
		else
			res = mats[i] * matrix.vec{1,1/Zl}
		end
		--if res[2] == 0 then 
		if i == 1 then
			print("res[1],res[2]",tostring(res[1]),tostring(res[2]),tostring(Zl))
			--print("Zin bin ",i,"res[2] == 0")
			popis[i] = Zl 
		else
			popis[i] = res[1]/res[2]
		end
	end
	return popis
end
function Mat2POUI(mats,Zls)
	local popis = {}
	for i,Zl in ipairs(Zls) do
		local res = mats[i] * matrix.vec{1,1/Zl}
		if i==1 then 
			popis[i] = 0 
		else
			popis[i] = 1/(res[2])
		end
	end
	return popis
end
function Mat2PIUO(mats,Zls)
	local popis = {}
	for i,Zl in ipairs(Zls) do
		local res = mats[i] * matrix.vec{Zl,1}
		if i==1 then 
			popis[i] = 0 
		else
			popis[i] = res[1]
		end
	end
	return popis
end


function ToCausal(ceps)
	--if true then return ceps end
	--fold to causal in time domain
	local nfft = #ceps
	print("#ceps",#ceps,nfft)
	local Nik=nfft/2 +1
	local ceps2 = {} --matrix.new(nfft, 1)
	ceps2[1] = 2*ceps[1]
	ceps2[Nik] = ceps[Nik]
	for k=2,Nik-1 do
		ceps2[k] = ceps[k] + ceps[nfft - k + 2]
	end
	for k=Nik+1,nfft do
		ceps2[k] = 0
	end
	return ceps2
end
function linearmap(s,e,ds,de,v)
	return ((de-ds)*(v-s)/(e-s)) + ds
end

function Impulse(popis,notcausal)
	local sq = matrix.new((#popis-1)*2, 1, |i| i == 1 and 1 or 0)
	print("#popis-1)*2",(#popis-1)*2)
	print("#popis-",#popis)
	local ft = num.fft(sq)
	ft[0] = complex.real(popis[1]) 
	for i=2,#popis-1  do
		--local Zi = complex.conj(popis[i])
		local Zi = popis[i]
		ft[i-1] = Zi 
	end
	ft[#popis-1] = complex.abs(popis[#popis]) 
	--[[
	local Nf = #popis-1
	local bincorte = math.floor(Nf*0.7)
	for i=bincorte,Nf do
		ft[i] = ft[i] * linearmap(bincorte,Nf,1,0,i)
	end
	--]]
	local res = num.fftinv(ft)
	if notcausal then
		return res,ft
	else
		return ToCausal(res),ft
	end
end

function Reflection(popis,radio)
	local area = radio*radio*math.pi
	local Z0 = RO_AIRE*C_AIRE/area
	local refle = {}
	for i,Zi in ipairs(popis) do
		if Zi == math.huge then 
			refle[i] = 1
		else
			refle[i] = (Zi-Z0)/(Zi+Z0);
		end
	end
	return refle
end
function Transmision(popis,radio)
	local area = radio*radio*math.pi
	local Z0 = RO_AIRE*C_AIRE/area
	local refle = {}
	for i,Zi in ipairs(popis) do
		if Zi == math.huge then 
			refle[i] = 2
		else
			refle[i] = 2*Zi/(Zi+Z0);
		end
	end
	return refle
end
function plotXY(X,Y,N)
	if N then
		X = TA(X)(1,N)
		Y = TA(Y)(1,N)
	end
	local Xm = matrix.vec(X)
	local Ym = matrix.vec(Y)
	local p = graph.plot('')
	local line = graph.xyline(Xm,Ym)
	p:add(line,"red",{{"stroke"}}) --,{"marker",size=10}
	p:show()
end
function Normalize(buf)
	local maxi = - math.huge
	for i=1,#buf do
		if math.abs(buf[i]) > maxi then maxi = math.abs(buf[i]) end
	end
	local invmaxi = 1/maxi
	for i=1,#buf do
		buf[i] = buf[i] * invmaxi
	end
	return buf
end
function writefile(buf,name)
	local file = io.open(name,"wb")
	for i=1,#buf do
		file:write(float2str(buf[i]):reverse())
	end
	file:close()
end
function writefileComplex(buf,name)
	local file = io.open(name,"wb")
	for i=1,#buf do
		file:write(float2str(complex.real(buf[i])):reverse())
		file:write(float2str(complex.imag(buf[i])):reverse())
	end
	file:close()
end
function Dibuja(segs)
	local function prepareplot2(data)
		local x ={}
		local y = {}
		local oldx = 0
		for i,v in ipairs(data) do
			x[i*2-1] = oldx
			y[i*2-1] = 0.5*(v.r1+v.r2)
			x[i*2] = oldx + v.l
			y[i*2] = 0.5*(v.r1+v.r2)
			oldx = oldx + v.l
		end
		--x[#x + 1]= oldx
		--y[#y + 1]= data[#data].r2
		return x,y
	end
	local X,Y = prepareplot2(segs)
	plotXY(X,Y)
end
function InfiniteZl(N,val)
	local val = val or math.huge
	local res = {}
	for i=1,N do
		res[#res +1] = val
	end
	return res
end
function radiofinal(ss)
	return 0.5*(ss[#ss].r1 + ss[#ss].r2)
end
function AyersDm(mats,Ris,radio)
	local area=radio*radio*math.pi;
	local  Zo=RO_AIRE*C_AIRE/area;
	local res = {}
	for i,Ri in ipairs(Ris) do
		local A = mats[i][1][1]
		local B = mats[i][1][2]
		--res[i]=A*(1+Ri) --+B*(1-Ri)/Zo;
		res[i]=(1+Ri)/A
		--res[i]=A*(1+Ri) -B*(1+Ri)/Zo;
	end
	return res
end
----------------------------draw
DatosTodo = TA(DatosCopa)..TA(Choke)..TA(Tubo)..TA(DatosBell)
--DatosTodo = TA(Tubo)..TA(DatosBell)
DatosTodo = data2seg(DatosTodo)
DatosTodo1 = TA(DatosCopa)..TA(Choke)
DatosTodo1 = data2seg(DatosTodo1)
DatosTodo2 = TA(Tubo)..TA(DatosBell)
DatosTodo2 = data2seg(DatosTodo2)


---------------------------------------------------------------------
function isnan(x) return x ~= x end

function TODO()
---[[
-- zin todo
sss = data2seg(TA(DatosCopa)..TA(Choke)..TA(Tubo)..TA(DatosBell))
sss=Segmenta(sss,1e-7)
Dibuja(sss)
freqs = GetFreqs(2^13)
mat = MatrixSegmentsScavone(sss,freqs)
Zrad = Radiation(freqs, sss[#sss].r2)
--Zrad = TA(freqs):Do(function(v) return LevineScavone(v,sss[#sss].r2) end)
Zin = Mat2Zin(mat,Zrad)
--Zin = Reflection(Zin,0.5*(sss[1].r1+sss[1].r2))
--
print("radio reflection",0.5*(sss[1].r1+sss[1].r2))
buf1 = Impulse(Zin,true)
graph.fxplot(|k| complex.abs(Zin[math.floor(k)]), 1, #Zin,"blue",1000)
--writefile(Normalize(buf1),dir.."zintodo.raw")
writefile(Normalize(buf1),dir.."Zin.raw")
Zin = Mat2POUI(mat,Zrad)
buf1 = Impulse(Zin,true)
writefile(Normalize(buf1),dir.."POUI.raw")
--]]
end
-------------------------------------------------------------------
function GI()
---[[
-- gi
sss=Segmenta(DatosTodo1,1e-8)
Dibuja(sss)
freqs = GetFreqs(2^13)
--mat = MatrixSegmentsF(sss,freqs)
mat = MatrixSegmentsScavone(sss,freqs)
Zrad = InfinitePiston(freqs, radiofinal(sss)) --sss[#sss].r2)
Zin = Mat2Zin(mat,Zrad)
buf1 = Impulse(Zin,true)
--writefile(Normalize(buf1),dir.."gi.raw")
graph.fxplot(|k| complex.abs(Zin[math.floor(k)]), 1, #Zin,"blue",1000)
writefile(buf1,dir.."gi.raw")
--]]
end
--------------------------------------------------------
function DM()
---[=[
--  d+
--sss=DatosTodo1 
sss = Segmenta(DatosTodo1,1e-7)
Dibuja(sss)
freqs = GetFreqs(2^13)
mat = MatrixSegmentsScavone(sss,freqs)
Zrad = InfinitePiston(freqs, radiofinal(sss))
Zin = Mat2POPI(mat,Zrad)
buf1,ft = Impulse(Zin,true)
writefile(buf1,dir.."dMas.raw")
--graph.fibars(|k| complex.abs(ft1[k]), 0, #ft1/2)
graph.fxplot(|k| complex.abs(ft[math.floor(k)]), 0, #ft/2,"blue",1000)
graph.fxplot(|k| Phase(ft[math.floor(k)]), 0, #ft/2,"red",1000)
--[[
Zin = Mat2Zin(mat,Zrad)
Zin = Reflection(Zin,0.5*(sss[1].r1 + sss[1].r2))
Zin2 = AyersDm(mat,Zin,0.5*(sss[1].r1 + sss[1].r2))
graph.fxplot(|k| complex.abs(Zin2[math.floor(k)]), 1, #Zin2/20,"red",1000)
buf1 = Impulse(Zin2)
writefile(buf1,dir.."dMas2.raw")
--]]
--]=]
end
------------------------------------------------------------------
function RB()
---[[
-- rb

sss=Segmenta(data2seg(DatosBell),1e-7)
radioini = 0.5*(sss[1].r1 + sss[1].r2)
Dibuja(sss)
freqs = GetFreqs(2^13)
mat = MatrixSegmentsScavone(sss,freqs)
Zrad = Radiation(freqs,radiofinal(sss))
Zin = Mat2Zin(mat,Zrad)
Zin = Reflection(Zin,radioini)
buf1 = Impulse(Zin,true)
--writefile(Normalize(buf1),dir.."rb.raw")
writefile(buf1,dir.."rb.raw")
graph.fxplot(|k| complex.abs(Zin[math.floor(k)]), 1, #Zin,"blue",1000)

Zin = Mat2POPI(mat,Zrad)
print("Zin[[1]]",Zin[1])
Zin[1] = 0
writefileComplex(Zin,dir.."popifinComplex.raw")
buf1 = Impulse(Zin,true)
writefile(Normalize(buf1),dir.."popifin.raw")
graph.fxplot(|k| complex.abs(Zin[math.floor(k)]), 1, #Zin,"blue",1000)

Zin = Mat2POPMas(mat,Zrad,radioini)
print("Zin[[1]]",Zin[1])
Zin[1] = 0
writefileComplex(Zin,dir.."popifinpmasComplex.raw")
buf1 = Impulse(Zin,true)
writefile(Normalize(buf1),dir.."popifinpmas.raw")
graph.fxplot(|k| complex.abs(Zin[math.floor(k)]), 1, #Zin,"blue",1000)
--]]
end
----------------------------------------------------------------
function RI()
---[[
-- ri
sss = TA(deepcopy(DatosTodo1)):reverse()
sss = sss:Do(function(v) v.r1,v.r2=v.r2,v.r1; return v end)
sss=Segmenta(sss,1e-7)
radioopen = 0.5*(sss[#sss].r1 + sss[#sss].r2)
Dibuja(sss)
freqs = GetFreqs(2^13)
mat = MatrixSegmentsScavone(sss,freqs)
print("mat[1]",mat[1])
Zrad = InfiniteZl(#freqs)
--Zrad = InfinitePiston(freqs,0.5*(sss[#sss].r1 + sss[#sss].r2))
Zin = Mat2Zin(mat,Zrad)
Zin = Reflection(Zin,0.5*(sss[1].r1 + sss[1].r2))
buf1 = Impulse(Zin,true)
--writefile(Normalize(buf1),dir.."ri.raw")
writefile(buf1,dir.."ri.raw")
graph.fxplot(|k| complex.abs(Zin[math.floor(k)]), 1, #Zin/20,"blue",1000)
--res=Mat[a]*(1.f+Ri)+Mat[b]*(1.f-Ri)/Zo;
Zin2 = AyersDm(mat,Zin,0.5*(sss[1].r1 + sss[1].r2))
graph.fxplot(|k| complex.abs(Zin2[math.floor(k)]), 1, #Zin2/20,"red",1000)
buf1 = Impulse(Zin2,true)
writefile(buf1,dir.."dm.raw")
--]]
end
-----------------------------------------------------------
function GetBin(x,zin,f)
	local bin = 1 + (#zin-1)*x/math.pi
	bin = math.floor(bin)
	return f(zin[bin])
end
function FluteRadiation()
local freqs = GetFreqs(2^13)
local Zrad,zin = Radiation(freqs,7.5*1e-3)

--weight = TA():Fill(#zin,function(i) if i < #zin/20 then return 1 else return clip(line(i,1*#zin/10,1,3*#zin/10,0),0,1) end end)
--local weight = TA():Fill(#zin,function(i) if i < #zin/10 then return ((#zin-1)/math.pi)^3 else return ((#zin-1)/((i-1)*math.pi))^3 end end)
weight = 1 / TA():series(#zin,0,math.pi/(#zin-1))^2
weight[1] = weight[2] 
maxpoles = 2
maxerror = 1e-8
--for np=2,maxpoles do
for np=maxpoles,maxpoles do
	HzminPH,chi = FilterFromFFT(np,np,1e-9,zin,weight)
	--HzminPH,chi = FilterFromFFTAmp(np,np,1e-9,zin,weight)
	--HzminPH,chi = FilterFromAmp((#zin - 1)*2,np,np,1e-9,model)
	HzminPH = HzminPH:Simplify()
	print(np," Chi ",chi)
	print(tostring(HzminPH))
	if chi < maxerror then break end
end

roots = HzminPH.Den:roots()
magroots = TA(roots):Do(complex.abs)
print("magroots",magroots)

impulse = HzminPH:filter({1},150)
--graph.fiplot(|x| impulse[x],#impulse -150,#impulse).title ="recons"
graph.fiplot(|x| impulse[x],150).title ="recons"

HzPDU=PhaseDelayUnwrap(HzminPH)
HzGD = HzminPH:GroupDelayFunction()

w = graph.window('v..')
p1 = graph.plot('Magnitude')
p1:addline(graph.fxline(|x| GetBin(x,zin,complex.abs), 0, math.pi),'blue')
p1:addline(graph.fxline(|x| HzminPH:magnitude(x), 0, math.pi),'red')
w:attach(p1, '1')
p1 = graph.plot('phases')
p1:addline(graph.fxline(|x| GetBin(x,zin,Phase), 0, math.pi),'blue')
p1:addline(graph.fxline(|x| HzminPH:phase(x) , 0, math.pi),'red')
--p1:addline(graph.fxline(|x| x>0 and (HzPDU(x)*x+math.pi)/x or 0, 0, math.pi),'yellow')
--p1:addline(graph.fxline(|x| HzPDU(x), 0, math.pi),'green')
--p1:addline(graph.fxline(|x| HzGD(x), 0, math.pi),'green')
w:attach(p1, '2')
lafreq = 65.3674
print("phasedel",HzPDU(2*math.pi*lafreq/44100))
end

FluteRadiation()
-- TODO()
-- DM()
--RB()
--RI()
--GI()
-------------------------------------------------------------

print("fin")