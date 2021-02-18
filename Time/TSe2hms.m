function time_string=TSe2hms(time_in_eeg_samples)

time=time_in_eeg_samples/1250;


time_string=sec2hms(time);