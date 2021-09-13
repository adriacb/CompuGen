import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
#!sudo apt-get install emboss


"""**2. Configure main files and directories**"""

# Generating the main directories files

def create_dirs(dirName):
  if not os.path.exists(dirName):
      os.mkdir(dirName)
      print("Directory " , dirName ,  " Created ")
  else:    
      print("Directory " , dirName ,  " already exists")

main = os.getcwd()

directories = {'bin': os.path.join(main, 'bin'),
               'data': os.path.join(main, 'data'),
               'docs': os.path.join(main, 'docs'),
               'images': os.path.join(main, 'images'),
               'stats': os.path.join(main, 'stats')}

for d in directories:
  print("Creating ",directories[d])
  create_dirs(directories[d])

# Generating the main report files

# Downloading the data
ucsc_file_to_download = 'mrna.fa.gz'
GenomeID = 'cb3'

def download_ucsc(GenomeID: str, file: str):
  try:
    print(f"Downloading from: http:hgdownload.soe.ucsc.edu/goldenPath/{GenomeID}/bigZips/{file} ")

    os.system(f"wget https://hgdownload.soe.ucsc.edu/goldenPath/{GenomeID}/bigZips/{file} -O {directories['data']}/{file}")

  except:
    print(f"Failed to download http:hgdownload.soe.ucsc.edu/goldenPath/{GenomeID}/bigZips/{file}")

  #os.system(f"ls {directories['data']}/{file}")

def download_ucsc_DB(GenomeID: str, file: str):
  try:
    print(f"Downloading from: http:hgdownload.soe.ucsc.edu/goldenPath/{GenomeID}/database/{file} ")

    os.system(f"wget https://hgdownload.soe.ucsc.edu/goldenPath/{GenomeID}/database/{file} -O {directories['data']}/{file}")

  except:
    print(f"Failed to download http:hgdownload.soe.ucsc.edu/goldenPath/{GenomeID}/database/{file}")


download_ucsc(GenomeID,ucsc_file_to_download)

"""**3. Simple analysis of mRNA sequences**"""

# Display basic information about sequences:

#Display Fasta head
os.system(f"gunzip -c {directories['data']}/{ucsc_file_to_download} | head -10")

#How many
os.system(f"gunzip -c {directories['data']}/{ucsc_file_to_download} | egrep -c '^>' ")

#Display information and save it into file
outfname = GenomeID + ucsc_file_to_download.split('.fa.gz')[0] + '.tbl'
outpath = os.path.join(directories['data'], outfname)
os.system(f"gunzip -c {directories['data']}/{ucsc_file_to_download} | infoseq -sequence fasta::stdin -outfile {directories['data']}/{outfname} -noheading -only -name -length -pgc")


"""**4. Visualizing the analysis of mRNA sequences**"""

df = pd.read_csv(str(outpath), delimiter=r'[ \t\n]+', names=['Sequence','Length', 'GC content (%)'])
print(df.head())
print(df.describe())


hist = df.hist(bins=len(df))
plt.show()


df.plot.scatter(x='GC content (%)',y='Length')
plt.show()



# Analyzing differences among chromosomes
ucsc_file_to_download = 'all_mrna.txt.gz'

download_ucsc_DB(GenomeID, ucsc_file_to_download)
os.system(f"gunzip -c {directories['data']}/{ucsc_file_to_download} | " + "gawk '{ print $10, $11, $15, $19 }' >"+ f"{directories['data']}/all_mrnabychr.tbl")
os.system("awk 'NR==FNR {ID[$1] = $0; next} {print ID[$2], $1" "$3" "$4;}' "+ f"{directories['data']}/{outfname}  {directories['data']}/all_mrnabychr.tbl > {directories['data']}/mrna_bychr_bylengc.tbl")

import statsmodels.api as sm

def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)


    m, b = np.polyfit(x, y, 1)
    
    ax.plot(x, m*x + b, '-', color='red')

    # now determine nice limits by hand:
    ax_histx.hist(x, bins=len(x))
    ax_histy.hist(y, bins=len(y), orientation='horizontal')

# start with a square Figure
fig = plt.figure(figsize=(8, 8))

# Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
# the size of the marginal axes and the main axes in both directions.
# Also adjust the subplot parameters for a square plot.
gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), 
  height_ratios=(2, 7),
  left=0.1, right=0.9, bottom=0.1, top=0.9,
  wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

# use the previously defined function
scatter_hist(df['GC content (%)'], df['Length'], ax, ax_histx, ax_histy)

plt.show()


# read the mRNA with the added extra information about chr

df_chr = pd.read_csv(f"{directories['data']}/mrna_bychr_bylengc.tbl", sep='\s+', names=['ID','Length','GC %','Chromosome']) 


#split chromosome column: i.e. +chrI2 into "+" "chrI2"
df_chr['Strand'] = df_chr['Chromosome'].str.extract('([+-])', expand=True)
df_chr['Chr'] = df_chr['Chromosome'].str.extract('([a-zA-z]+)', expand=True)

print(df_chr.head())


#HISTOGRAMS: GC and length

g = sns.displot(data=df_chr, x="GC %", hue="Chr", multiple="stack",
col="Strand",
    kind="hist", height=4, aspect=.7)
plt.show()


g2 = sns.displot(data=df_chr, x="Length", hue="Chr", multiple="stack",
col="Strand",
    kind="hist", height=4, aspect=.7, bins=len(df_chr))
g2.set(yscale="log")
plt.show()



#BOXPLOTS: GC and length

g3 = sns.catplot(data=df_chr, y="GC %", x='Chr', hue="Strand", 
    kind="box")
plt.show()

g4 = sns.catplot(data=df_chr, y="Length", x='Chr', hue="Strand", 
    kind="box")
plt.show()



###########################FALTA CANVIAR COLOR PER CHR


# start with a square Figure
fig = plt.figure(figsize=(8, 8))

# Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
# the size of the marginal axes and the main axes in both directions.
# Also adjust the subplot parameters for a square plot.
gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), 
  height_ratios=(2, 7),
  left=0.1, right=0.9, bottom=0.1, top=0.9,
  wspace=0.05, hspace=0.05)

ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

# use the previously defined function
scatter_hist(df_chr['GC %'], df['Length'], ax, ax_histx, ax_histy)

plt.show()