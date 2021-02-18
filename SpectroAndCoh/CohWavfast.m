function out=CohWavfast(pow,wav)

[y0 f0] =mtchd(wav,2^9,1250,2^7);
out{1}=y0;
out{2}=f0;