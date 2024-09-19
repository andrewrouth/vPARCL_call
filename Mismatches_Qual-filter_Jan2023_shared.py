import re
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("Pileup", help="Description")
parser.add_argument("Output", help="Description")
parser.add_argument("Quals", help="Min PHRED score")
parser.add_argument("--MinCov", help= "Description")
args = parser.parse_args()

QualLim = int(args.Quals) + 33
if args.MinCov:
    MinCov = int(args.MinCov)
else:
    MinCov = 1
Out = open(str(args.Output), 'w')
Out.write("Ref\tPosition\tRefNuc\tCoverage\tA\tT\tG\tC\tErrorrate\t\n")

All_Data = {}
Nucs = ['A','T','G','C','a','t','g','c']
Bases = ['.',',','A','T','G','C','N']
InDels = ['-', '+']
with gzip.open(str(args.Pileup), 'rt') as In:
        Data = In.readline().split()
        while Data:
            Ref = Data[0]
            if Ref not in All_Data:
                All_Data[Ref] = []
                All_Data[Ref].append(["Ref", "Coverage", "A", "T", "G", "C"])
            All_Data[Ref].append([0,0,0,0,0,0])
            Position = int(Data[1])
            RefNuc = Data[2].upper()
            Seqs = Data[4].upper()
            Quals = Data[5]
            while len(All_Data[Ref]) <= Position:
                    All_Data[Ref].append([0,0,0,0,0,0])
            Seqs_filt = []
            if '^' in Seqs:
                    x,y = Seqs.split('^',1)
                    y = y.split('^')
                    y = [i[1:] for i in y]
                    y = ''.join(y)
                    Seqs = x+y
            else:
                    pass
            for j in InDels:
                    if j in Seqs:
                            Seqs,y = Seqs.split(j,1)
                            y = y.split(j)
                            for k in y:
                                Ks = re.split('(\d*)',k,1)
                                Seqs += Ks[2][int(Ks[1]):]
                    else:
                            pass
            Seqs = Seqs.replace('$','')
            #Seqs = Seqs.replace('^','')
            Seqs = Seqs.replace('~','')
            if len(Seqs) != len(Quals):
                print(Position, Data)
            n=0 
            Seqs_filt = []
            for i in range(len(Seqs)):
                Qual = Quals[i]
                if ord(Qual) > QualLim:
                    if Seqs[i] != '>' and Seqs[i] != '<':
                        Seqs_filt += Seqs[i]
                    else:
                        pass
                else:
                    pass
            Coverage = float(len(Seqs_filt))
            All_Data[Ref][Position][0:2] = [RefNuc, Coverage]
            if RefNuc == "A":
                    for i in range(len(Seqs_filt)):
                            if Seqs_filt[i] != '.' and Seqs_filt[i] != ',':    
                                         if Seqs_filt[i] == "T":
                                                All_Data[Ref][int(Position)][3] += 1
                                         elif Seqs_filt[i] == "G":
                                                All_Data[Ref][int(Position)][4] += 1
                                         elif Seqs_filt[i] == "C":
                                                All_Data[Ref][int(Position)][5] += 1
                                         else:
                                                 pass
                            else:
                                 pass
                            n+=1
            elif RefNuc == "T":
                        for i in range(len(Seqs_filt)):
                            if Seqs_filt[i] != '.' and Seqs_filt[i] != ',':
                                         if Seqs_filt[i] == "A":
                                                All_Data[Ref][int(Position)][2] += 1
                                         elif Seqs_filt[i] == "G":
                                                All_Data[Ref][int(Position)][4] += 1
                                         elif Seqs_filt[i] == "C":
                                                All_Data[Ref][int(Position)][5] += 1
                                         else:
                                                 pass
                            else:
                                 pass
                            n+=1
            elif RefNuc == 'G':
                  for i in range(len(Seqs_filt)):
                            if Seqs_filt[i] != '.' and Seqs_filt[i] != ',':
                                         if Seqs_filt[i] == "A":
                                                All_Data[Ref][int(Position)][2] += 1
                                         elif Seqs_filt[i] == "T":
                                                All_Data[Ref][int(Position)][3] += 1
                                         elif Seqs_filt[i] == "C":
                                                All_Data[Ref][int(Position)][5] += 1
                                         else:
                                                 pass
                            else:
                                 pass
                            n+=1
            elif RefNuc == "C":
                        for i in range(len(Seqs_filt)):
                            if Seqs_filt[i] != '.' and Seqs_filt[i] != ',':
                                         if Seqs_filt[i] == "A":
                                                All_Data[Ref][int(Position)][2] += 1
                                         elif Seqs_filt[i] == "T":
                                                All_Data[Ref][int(Position)][3] += 1
                                         elif Seqs_filt[i] == "G":
                                                All_Data[Ref][int(Position)][4] += 1
                                         else:
                                                 pass
                            else:
                                 pass
                            n+=1
            else:
                    pass
            Mismatches = 0
            for i in range(2,6):
                  Mismatches += All_Data[Ref][int(Position)][i]
            if Coverage >= MinCov:
                ErrorRate = Mismatches/Coverage
            else:
                ErrorRate = 0 #'NA'
            Out.write(str(Ref) + '\t'
              + str(Position) + '\t'
              + str(RefNuc) + "\t"
              + str(int(Coverage)) + "\t"
              + str(All_Data[Ref][int(Position)][2]) + "\t"
              + str(All_Data[Ref][int(Position)][3]) + "\t"
              + str(All_Data[Ref][int(Position)][4]) + "\t"
              + str(All_Data[Ref][int(Position)][5]) + "\t"
              + str(ErrorRate) + '\t\n')
            Data = In.readline().split()
Out.close()


   



        

