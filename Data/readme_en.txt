The file registros_EIS.m contains a MATLAB-formatted matrix with the following column order:
ID | device | electrode | sample | freq_range | excitation | v_i | origin | eis

Significado de cada uno de las etiquetas de cada columna:

ID: Sequential number of the measurement, assigned randomly based on the order in which the data was loaded.
device: Model and brand of the device used to perform the bioimpedance measurement.
electrode: Electrode configuration used during the bioimpedance measurement.
sample: Biological sample on which the bioimpedance measurement was taken.
freq_range: Frequency range over which the frequency sweep was conducted.
excitation: Amplitude of the excitation signal.
v_i: Indicates whether the excitation was performed using voltage (v) or current (i).
origin: Author or experiment from which the dataset was sourced.
eis: Spectroscopic bioimpedance measurement values.

Values and/or format of each column label.
ID: Unique and unrepeatable numeric value.

device: Numeric value that identifies a device and model.
	0: No data
	1: Solartron
	2: HF2IS Zurich Instruments
	3: Xitron 4001B
	4: SFB7 Impedimed

electrode:
	0: No data
	2: bipolar
	3: tripolar
	4: tetrapolar

sample:
	0: No data
	1: human
	2: vegetable
	3: meat
	4: animal
	5: AgarAgarPhantom

freq_range:
	0: No data
	1: 1Hz-65kHz
	2: 400Hz-1MHz
	3: 1kHz-1.2MHz
	4: 5kHz-1MHz
	5: 3kHz-1MHz

excitation:
	
v_i:

origin:
	0: desconocido
	1: Dell'Osa et al, 2019
	2: Gonzalez Correa, 2018
	3: DeLorenzo et al, 1998
	4: Ayllon et al, 
	5: Bertemes Filho (no Publish measurements)

eis: It is a numeric matrix. The number of rows depends on the data origin. Each row corresponds to a different measurement point:
	Column 1: Frequency
	Column 2: Real part of the bioimpedance measurement
	Column 3: Imaginary part of the bioimpedance measurement
	Column 4: Magnitude or modulus of the bioimpedance measurement
	Column 5: Phase of the bioimpedance measurement

Each row contains a different measurement.