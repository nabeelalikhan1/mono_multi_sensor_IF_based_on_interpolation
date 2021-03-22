
num=2;
tic
for i=1:100
[findex] = FAST_IF(Sig,win_length, num, delta,L,0,0)*2*SampFreq;
end
toc

tic
for i=1:100
findex =FASTEST_IF(Sig,win_length, num, delta,L,0,0,4,length(Sig))*2*SampFreq;
end
toc

tic
for i=1:100
findex =FASTEST_IF(Sig,win_length, num, delta,L,0,0,8,length(Sig))*2*SampFreq;
end
toc

tic
for i=1:100
findex =FASTEST_IF(Sig,win_length, num, delta,L,0,0,16,length(Sig))*2*SampFreq;
end
toc