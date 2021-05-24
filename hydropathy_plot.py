from Bio import SeqIO
import matplotlib.pyplot as plt


record = SeqIO.read() #read fasta-file

# https://doi.org/10.1016/0022-2836(82)90515-0 values above 1.6 indicate membranespanning sequences (using a window of about 20 aminoacids)
Kyte_Doolittle = {
                "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, 
                "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
                "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
                "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3
                }

length = len(record.seq)
x = []
y = []

for idx, aa in enumerate(record.seq):
    s = 0.0
    if idx < length-19:
        x.append(idx)
        for i in record.seq[idx:idx+20]:
            s = s + Kyte_Doolittle[i]
        y.append(s)

    elif idx == length-19:
        x.append(idx)
        for i in record.seq[idx:idx+20]:
            s = s + Kyte_Doolittle[i]
        y.append(s)
        break

plt.figure()
ax = plt.axes()
plt.rcParams["font.family"] = "Arial"
plt.plot(x,y)
plt.plot([0,length],[1.6, 1.6], "r")
ax.set(xlim=(0, length-19))
plt.ylabel(r"$\sum_{i=0}^{20}\ hydropathy\ of\ AA_i$", fontsize=16, math_fontfamily='stixsans')
plt.xlabel(r"$AA_i$", fontsize=16, math_fontfamily='stixsans')
plt.show()